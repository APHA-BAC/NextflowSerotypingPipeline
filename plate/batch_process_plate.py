import subprocess
import os
import glob
import argparse
import logging
import signal
from textwrap import dedent


DEFAULT_READS_DIRECTORY = os.path.expanduser('~/wgs-reads')
DEFAULT_KMER_URI = "s3://s3-ranch-046/KmerID_Ref_Genomes"


def run(cmd, record_output=False):
    """
        Run a command and assert that the process exits with a non-zero
        exit code. If record_output=True, then the stdout of the
        subcommand is logged

        Parameters:
            cmd (list): List of strings defining the command, see
            (subprocess.run in python docs)
    """
    ps = subprocess.run(cmd, capture_output=True)
    returncode = ps.returncode
    if record_output:
        logging.info(format_subprocess_output(ps.stdout))
    if returncode:
        raise Exception(dedent(f"""
                                   *****
                                   cmd '{(" ").join(cmd)}' failed with exit code {returncode}
                                   {format_subprocess_output(ps.stderr)}
                                   *****
                                """))


def format_subprocess_output(output):
    """
        Removes any output which proceeds '\r' without a trailing '\n'
        character, i.e. removes any output with carriage return and
        without newline. Lots of AWS CLI commands produce output which
        overwrites itself by using '\r'. This output looks ugly in a
        logfile so we remove it.
    """
    return "\n".join([split_line.split("\r")[-1] for split_line in
                      output.decode().split("\n")])


def run_pipeline(plate_name, **kwargs):
    """ Run the Salmonella pipeline """
    run(["/root/nextflow/nextflow", "SCE3_pipeline_update.nf", "--runID",
         plate_name], **kwargs)


def download_s3(s3_uri, destination, **kwargs):
    """
        Recursively download a S3 Object
    """
    run(["aws", "s3", "cp", "--recursive", s3_uri, destination], **kwargs)


def upload_s3(file_path, s3_destination, **kwargs):
    """
        Uploads a file to S3
    """
    run(["aws", "s3", "cp", file_path, s3_destination, "--acl",
         "bucket-owner-full-control"], **kwargs)


def s3_object_release_date(s3_uri):
    """
        Date s3 object was published. Returns a 3 element list with
        format [year, month, day]
    """
    # Retrieve metadata from S3
    ls_cmd = f"aws s3 ls {os.path.join(s3_uri, '')}"
    contents = [x.decode("utf-8") for x in
                subprocess.check_output(ls_cmd, shell=True).splitlines()]

    # Extract date
    return contents[0].split()[0].split("-")


def s3_uri_to_plate_name(s3_uri):
    """
        Convert a S3 URI from CSU to a plate name with consistent naming
        convention
    """
    # Remove trailing slash
    s3_uri = s3_uri.rstrip("/")

    # Format
    year, month, day = s3_object_release_date(s3_uri)
    run_name = s3_uri.split("/")[-1]

    return f"{day}{month}{year[-2:]}_APHA_{run_name}"


def rename_fastq_file(filepath):
    """
        Rename a fastq file from CSU's convention to BGE
    """
    # Parse
    directory = os.path.dirname(filepath)
    filename = os.path.join(os.path.basename(filepath), "")
    sample_name = filename.split("_")[0]
    # Determine Read Number
    if "_R1" in filename:
        read_number = 1
    elif "_R2" in filename:
        read_number = 2
    else:
        raise Exception("Unable to determine read pair number: ", filename)
    # Rename
    renamed = f"{directory}/{sample_name}_R{read_number}.fastq.gz"
    os.rename(filepath, renamed)


def run_plate(reads_uri, reads_dir, results_uri, kmer_uri):

    """
        Download, process and store a plate of raw Salmonella data
    """
    # Set paths
    plate_name = s3_uri_to_plate_name(reads_uri)
    plate_reads_dir = os.path.join(reads_dir, plate_name)

    # Download reference genomes from s3
    logging.info(f"Downloading KmerID reference genomes: {kmer_uri}\n")
    download_s3(kmer_uri, "/root/KmerID_Ref_Genomes")

    # Download reads
    logging.info(f"Downloading reads: {reads_uri}\n")
    download_s3(reads_uri, plate_reads_dir, record_output=True)

    # Rename fastq files
    logging.info(f"Renaming fastq files: {reads_dir}\n")
    for filepath in glob.glob(plate_reads_dir + '/*.fastq.gz'):
        rename_fastq_file(filepath)

    logging.info("Running Nextflow pipeline\n")
    run_pipeline(plate_name, record_output=True)

    # Upload results to s3
    TableFile_name = plate_name + "_SummaryTable_plusLIMS.csv"
    summaryTable_path = os.path.join("~/wgs-results/", plate_name,
                                     TableFile_name)
    summaryTable_path = os.path.expanduser(summaryTable_path)
    logging.info(f"Uploading results: {results_uri}\n")
    upload_s3(summaryTable_path, os.path.join(results_uri, TableFile_name),
              record_output=True)


def upload_logfile(results_uri):
    log_uri = os.path.join(results_uri, "batch_process_plate.log")
    logging.info(f"Uploading log file: {log_uri}")
    upload_s3(log_file_path, log_uri)


if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(
        description="run pipeline on a routine Salmonella Plate")
    parser.add_argument("reads_uri",
                        help="s3 uri corresponding to the fastq plate to run")
    parser.add_argument("results_uri",
                        help="s3 uri where results are to be uploaded")
    parser.add_argument("--reads-dir", default=DEFAULT_READS_DIRECTORY,
                        help="local directory for storing reads")
    parser.add_argument("--kmer_uri", default=DEFAULT_KMER_URI,
                        help="s3 uri of KmerID reference genomes")

    args = parser.parse_args()

    # setup logging
    log_file_path = os.path.expanduser("~/batch_process_plate.log")
    logging.basicConfig(level=logging.INFO, format="%(message)s",
                        handlers=[logging.StreamHandler(),
                                  logging.FileHandler(log_file_path)])
    # setup handler to upload log if batch job times out
    signal.signal(signal.SIGTERM,
                  upload_logfile(
                      f"{args.results_uri.rstrip('/')}_timeout_error"))
    # Run
    try:
        run_plate(args.reads_uri, args.reads_dir, args.results_uri,
                  args.kmer_uri)
        results_uri = args.results_uri
    except Exception as e:
        # if the run fails, append "_failed" to the results_uri
        results_uri = f"{args.results_uri.rstrip('/')}_failed"
        logging.exception(e)
        # re-raise the caught exception
        raise e
    # the finally block runs before re-raising 'e'.
    finally:
        upload_logfile(results_uri)
