import subprocess
import os
import glob
import argparse
import logging
import signal
import tempfile
import pandas as pd
from textwrap import dedent
import update_master_sum
import time


DEFAULT_READS_DIRECTORY = os.path.expanduser('~/wgs-reads')
DEFAULT_KMER_URI = "s3://s3-ranch-046/KmerID_Ref_Genomes"
DEFAULT_MASTER_SUMMARY_URI = "s3://s3-ranch-050/master_summary.csv"


class TimeoutHandler:
    def __init__(self, results_uri):
        self.results_uri = results_uri
        signal.signal(signal.SIGTERM, self.handler)

    def handler(self, *_):
        raise TimeoutError


def run(cmd, record_output=False):
    """
        Run a command and assert that the process exits with a non-zero
        exit code. If record_output=True, then the stdout of the
        subcommand is logged

        Parameters:
            cmd (list): List of strings defining the command, see
            (subprocess.run in python docs)
    """
    tf = tempfile.NamedTemporaryFile(delete=False)
    try:
        with open(tf.name, mode="wb") as stdout_fo:
            ps = subprocess.Popen(cmd, stdout=stdout_fo, stderr=subprocess.PIPE)
            ps.wait()
        with open(tf.name, mode="rb") as stdout_fo:
            if record_output:
                logging.info(format_subprocess_output(stdout_fo.read()))
            if ps.returncode:
                raise Exception(dedent(f"""
                                        *****
                                        cmd '{(" ").join(cmd)}' failed with exit code {ps.returncode}
                                        {format_subprocess_output(ps.stderr.read())}
                                        *****
                                        """))
    except TimeoutError as e:
        with open(tf.name, mode="rb") as stdout_fo:
            ps.kill()
            logging.info(format_subprocess_output(stdout_fo.read()))
            raise e
    finally:
        tf.close()
        os.unlink(tf.name)


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
         plate_name, "--batchRun", "True"], **kwargs)


def download_s3(s3_uri, destination, recursive=False, **kwargs):
    """
        Download a S3 Object
    """
    if recursive:
        run(["aws", "s3", "cp", "--recursive", s3_uri, destination], **kwargs)
    else:
        run(["aws", "s3", "cp", s3_uri, destination], **kwargs)


def upload_s3(file_path, s3_destination, **kwargs):
    """
        Uploads a file to S3
    """
    run(["aws", "s3", "cp", file_path, s3_destination, "--acl",
         "bucket-owner-full-control"], **kwargs)

def upload_fasta(assemblies_dir):

    logging.info(f"Uploading fasta")
    run(["aws", "s3", "cp", "--recursive", assemblies_dir, "s3://s3-ranch-055/assemblies", "--acl", "bucket-owner-full-control"])

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


def update_master_summary(TableFile_name,
                          master_sum_uri=DEFAULT_MASTER_SUMMARY_URI):
    """
        Updates the master summary table stored in s3-ranch-050,
        appending the latest results to a master CSV containing all
        historical results
    """
    df_new_sum = pd.read_csv(TableFile_name)
    download_s3(master_sum_uri, "master_sum.csv", record_output=True)
    df_new_sum.to_csv("master_sum.csv", mode="a", header=False, index=False)
    upload_s3("master_sum.csv", master_sum_uri, record_output=True)


def run_plate(reads_uri, reads_dir, results_uri, kmer_uri):

    """
        Download, process and store a plate of raw Salmonella data
    """
    # Set paths
    plate_name = s3_uri_to_plate_name(reads_uri)
    plate_reads_dir = os.path.join(reads_dir, plate_name)
    assemblies_dir = os.path.join("/root/wgs-results/{}/assemblies/".format(plate_name))


    # Download reference genomes from s3
    logging.info(f"Downloading KmerID reference genomes: {kmer_uri}\n")
    download_s3(kmer_uri, "/root/KmerID_Ref_Genomes", recursive=True)

    # Download reads
    logging.info(f"Downloading reads: {reads_uri}\n")
    download_s3(reads_uri, plate_reads_dir, recursive=True, record_output=True)

    # Rename fastq files
    logging.info(f"Renaming fastq files: {reads_dir}\n")
    for filepath in glob.glob(plate_reads_dir + '/*.fastq.gz'):
        rename_fastq_file(filepath)

    logging.info("Running Nextflow pipeline\n")

    try:
        run_pipeline(plate_name, record_output=True)
        with open("{}_finished.scemail".format(plate_name), 'w') as text_file:
            text_file.write("The following plate has completed:\n{}".format(plate_name))
        run(['aws', 's3', 'cp', './{}_finished.scemail'.format(plate_name), 's3://s3-scemail-3cfh-salmonella-serotyping-pipeline-1-0-0/salmonella-serotyping-pipeline/request/', 
            '--acl', 'bucket-owner-full-control'])
        # Upload results to s3
        TableFile_name = plate_name + "_SummaryTable_plusLIMS.csv"
        summaryTable_path = os.path.join("~/wgs-results/", plate_name,
                                         TableFile_name)
        summaryTable_path = os.path.expanduser(summaryTable_path)
        logging.info(f"Uploading results: {results_uri}\n")
        upload_s3(summaryTable_path, os.path.join(results_uri, TableFile_name),
                  record_output=True)

        # Update master summary table
        # logging.info(F"Updating master summary table: {DEFAULT_MASTER_SUMMARY_URI}")
        run(["aws", "s3", "cp", "s3://s3-ranch-055/master_summary.csv", "/root/"])
        update_master_sum.update_summary("/root/master_summary.csv",summaryTable_path)
        run(["aws", "s3", "cp", "/root/master_summary.csv", "s3://s3-ranch-055/", "--acl", "bucket-owner-full-control"])
        logging.info(f"New sum path: {summaryTable_path}\n")
        
        # upload_fasta(assemblies_dir)


    except:
        with open("{}_failed.scemail".format(plate_name), 'w') as text_file:
            text_file.write("The following plate has failed: \n{} \nPlease see log file for more details".format(plate_name))
        time.sleep(10)
        run(['aws', 's3', 'cp', './{}_failed.scemail'.format(plate_name), 's3://s3-scemail-3cfh-salmonella-serotyping-pipeline-1-0-0/salmonella-serotyping-pipeline/request/', 
            '--acl', 'bucket-owner-full-control'])

    


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
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s",
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        handlers=[logging.StreamHandler(),
                                  logging.FileHandler(log_file_path)])
    # setup handler to upload log if batch job times out
    timeout_handler = \
        TimeoutHandler(f"{args.results_uri.rstrip('/')}_timeout_error")
    # Run
    try:
        run_plate(args.reads_uri, args.reads_dir, args.results_uri,
                  args.kmer_uri)
        results_uri = args.results_uri
    except Exception as e:
        if isinstance(e, TimeoutError):
            # append "_timeout" to the results_uri
            results_uri = f"{args.results_uri.rstrip('/')}_timeout_error"
            logging.info("\nAWS batch job timed out\n")
        else:
            # append "_failed" to the results_uri
            results_uri = f"{args.results_uri.rstrip('/')}_failed"
            logging.exception(e)
        # re-raise the caught exception
        raise e
    # the finally block runs before re-raising 'e'.
    finally:
        upload_logfile(results_uri)
