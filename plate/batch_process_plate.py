import subprocess
import os
import glob
import argparse
import logging

DEFAULT_READS_DIRECTORY = os.path.expanduser('~/root/wgs-reads')


def run(cmd):
    """ Run a command and assert that the process exits with a non-zero
        exit code.

        Parameters:
            cmd (list): List of strings defining the command, see
            (subprocess.run in python docs)
    """
    # TODO: store stdout to a file
    returncode = subprocess.run(cmd).returncode

    if returncode:
        raise Exception("""*****
            %s
            cmd failed with exit code %i
        *****""" % (cmd, returncode))


def run_pipeline(plate_name):
    """ Run the Salmonella pipeline using docker """
    run(["/root/nextflow/nextflow", "SCE3_pipeline_update.nf", "--local",
         plate_name])


def download_s3(s3_uri, destination):
    """ Recursively download a S3 Object """
    run(["aws", "s3", "cp", "--recursive", s3_uri, destination])


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


def upload_s3(file_path, s3_destination):
    """
        Uploads a file to S3
    """
    run(["aws", "s3", "cp", file_path, s3_destination])


def download_kmerid():
    """
        Downloads reference genomes from s3
    """
    run(["aws", "s3", "cp", "--acl", "bucket-owner-full-control", "--recursive",
         "s3://s3-ranch-046/KmerID_Ref_Genomes", "/root/KmerID_Ref_Genomes/"])


def run_plate(reads_uri, reads_dir, results_uri):

    """
        Download, process and store a plate of raw Salmonella data
    """
    download_kmerid()

    # Download reads
    download_s3(reads_uri, reads_dir)

    # Rename fastq files
    for filepath in glob.glob(reads_dir + '/*.fastq.gz'):
        rename_fastq_file(filepath)

    run_pipeline(reads_dir)

    plate_name = s3_uri_to_plate_name(reads_uri)
    TableFile_name = plate_name + "_SummaryTable_plusLIMS.csv"
    summaryTable_path = os.path.join("~/root/wgs-results/", plate_name,
                                     TableFile_name)
    summaryTable_path = os.path.expanduser(summaryTable_path)
    upload_s3(summaryTable_path, os.path.join(results_uri, "."))


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

    args = parser.parse_args()

    log_file_path = os.path.expanduser("~/root/batch_process_plate.log")
    logging.basicConfig(level=logging.INFO,
                        handlers=[logging.FileHandler(log_file_path)])

    # Run
    try:
        run_plate(args.reads_uri, args.reads_dir, args.results_uri)
        upload_s3(log_file_path, os.path.join(args.results_uri, "."))
    except Exception as e:
        logging.exception(e)
        upload_s3(log_file_path, f"{args.results_uri.rstrip('/')}_failed/.")
