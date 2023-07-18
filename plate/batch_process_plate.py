import subprocess
import os
import glob
import argparse
import logging
from textwrap import dedent
import boto3
import os 

def downloadDirectoryFroms3(bucketName, remoteDirectoryName, dest):
    s3_resource = boto3.resource('s3')
    bucket = s3_resource.Bucket(bucketName) 
    for obj in bucket.objects.filter(Prefix = remoteDirectoryName):
        if not os.path.exists(os.path.dirname(obj.key)):
            os.makedirs(os.path.dirname(obj.key))
        bucket.download_file(obj.key, os.path.join(dest, os.path.basename(obj.key))) # save to same path

DEFAULT_READS_DIRECTORY = os.path.expanduser('~/wgs-reads')
DEFAULT_KMER_URI = "s3://s3-ranch-046/KmerID_Ref_Genomes"


def run(cmd, *args, **kwargs):
    """
        Run a command and assert that the process exits with a non-zero
        exit code. See python's subprocess.run command for args/kwargs.
        If capture_output=True, then the stdout of the subcommand is
        logged

        Parameters:
            cmd (list): List of strings defining the command, see
            (subprocess.run in python docs)
    """
    ps = subprocess.run(cmd, *args, **kwargs)
    returncode = ps.returncode
    if "capture_output" in kwargs and kwargs["capture_output"]:
        logging.info(ps.stdout.decode().strip('\n'))
#    if returncode:
#        raise Exception(dedent(f"""
#                                   *****
#                                   cmd '{(" ").join(cmd)}' failed with exit \
#                                   code {returncode}
#                                   *****
#                                """))


def run_pipeline(plate_name, **kwargs):
    """ Run the Salmonella pipeline using docker """
    run(["/root/nextflow/nextflow", "SCE3_pipeline_update.nf", "--local",
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


def s3_object_release_date(s3_key):
    """
        Date s3 object was published. Returns a 3 element list with
        format [year, month, day]
    """
    # Retrieve metadata from S3
    ls_cmd = f"aws s3 ls {s3_key}/"
    contents = [x.decode("utf-8") for x in
                subprocess.check_output(ls_cmd, shell=True).splitlines()]

    # Extract date
    return contents[0].split()[0].split("-")


def s3_uri_to_plate_name(s3_key):
    """
        Convert a S3 URI from CSU to a plate name with consistent naming
        convention
    """
    # Remove trailing slash
    s3_key = s3_key.strip('/')

    # Format
    year, month, day = s3_object_release_date(s3_key)
    run_name = s3_key.split("/")[-1]

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
    # Download reference genomes from s3
    logging.info(f"Downloading KmerID reference genomes: {kmer_uri}")
    download_s3(kmer_uri, "/root/KmerID_Ref_Genomes")
    run(["df", "-h"], capture_output=True)

    # Download reads
    logging.info(f"Downloading reads: {reads_uri}")
    #downloadDirectoryFroms3("s3-csu-001", "FZ2000/M01765_0638/", reads_dir)
    download_s3(reads_uri, reads_dir, capture_output=True)

    # Rename fastq files
    logging.info(f"Renaming fastq files: {reads_dir}")
    for filepath in glob.glob(reads_dir + '/*.fastq.gz'):
        rename_fastq_file(filepath)

    logging.info("Running Nextflow pipeline")
    run_pipeline(reads_dir, capture_output=True)

    # Upload results to s3
    plate_name = s3_uri_to_plate_name(reads_uri)
    TableFile_name = plate_name + "_SummaryTable_plusLIMS.csv"
    summaryTable_path = os.path.join("~/wgs-results/", plate_name,
                                     TableFile_name)
    summaryTable_path = os.path.expanduser(summaryTable_path)
    logging.info(f"Uploading results: {results_uri}")
    upload_s3(summaryTable_path, os.path.join(results_uri, TableFile_name),
              capture_output=True)


if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(
        description="run pipeline on a routine Salmonella Plate")
    parser.add_argument("-i", "--reads_uri",
                        help="s3 uri corresponding to the fastq plate to run")
    parser.add_argument("-o", "--results_uri",
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

    # Run
    try:
        run_plate(args.reads_uri, args.reads_dir, args.results_uri,
                  args.kmer_uri)
    except Exception as e:
        # if the run fails, append "_failed" to the results_uri
        results_uri = \
            f"{args.results_uri.rstrip('/')}_failed"
        logging.exception(e)
        # re-raise the caught exception
        raise e
    # the finally block runs before re-raising 'e'.
    finally:
        # upload log file
        log_uri = os.path.join(results_uri, "batch_process_plate.log")
        logging.info(f"Uploading log file: {log_uri}")
        upload_s3(log_file_path, log_uri)
