import subprocess
import os
import glob
import argparse

# TODO: Rename directories to BGE defaults
DEFAULT_READS_DIRECTORY = os.path.expanduser('~/wgs-reads/')
DEFAULT_RESULTS_DIRECTORY = os.path.expanduser('~/wgs-results/')
DEFAULT_IMAGE = "jguzinski/salmonella-seq:prod"

def run(cmd):
    """ Run a command and assert that the process exits with a non-zero exit code.
        
        Parameters:
            cmd (list): List of strings defining the command, see (subprocess.run in python docs)
    """
    # TODO: store stdout to a file
    returncode = subprocess.run(cmd).returncode

    if returncode:
        raise Exception("""*****
            %s
            cmd failed with exit code %i
        *****""" % (cmd, returncode))

def run_pipeline(reads, results, plate_name, image=DEFAULT_IMAGE):
    """ Run the Salmonella pipeline using docker """
    run(["sudo", "docker", "pull", image])
    run([
        "sudo", "docker", "run", "--rm", "-it", 
        "-v", f"{reads}:/root/WGS_Data/{plate_name}/",
        "-v", f"{results}:/root/WGS_Results/{plate_name}/", 
        image,
        "/root/nextflow/nextflow", "SCE3_pipeline_update.nf",
        "--runID", plate_name
    ])

def download_s3(s3_uri, destination):
    """ Recursively download a S3 Object """
    run(["aws", "s3", "cp", "--recursive",
        s3_uri,
        destination
    ])

def s3_object_release_date(s3_key):
    """ Date s3 object was published. Returns a 3 element list with format [year, month, day] """

    # Retrieve metadata from S3
    ls_cmd = f"aws s3 ls {s3_key}/"
    contents = [x.decode("utf-8") for x in subprocess.check_output(ls_cmd, shell=True).splitlines()]
    
    # Extract date
    return contents[0].split()[0].split("-")

def s3_uri_to_plate_name(s3_key):
    """ Convert a S3 URI from CSU to a plate name with consistent naming convention """
    
    # Remove trailing slash
    s3_key = s3_key.strip('/')

    # Format
    year, month, day = s3_object_release_date(s3_key)
    run_name = s3_key.split("/")[-1]

    return f"{day}{month}{year[-2:]}_APHA_{run_name}"

def rename_fastq_file(filepath):
    """ Rename a fastq file from CSU's convention to BGE """

    # Parse
    directory = os.path.dirname(filepath)
    filename = os.path.basename(filepath) + '/'
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

def run_plate(s3_uri, reads_dir, results_dir):
    """ Download, process and store a plate of raw Salmonella data """

    # Add trailing slash to directory names
    reads_dir = os.path.join(reads_dir, '')
    results_dir = os.path.join(results_dir, '')

    # Storage paths
    plate_name = s3_uri_to_plate_name(s3_uri)
    plate_reads_dir = reads_dir + plate_name + '/'
    plate_results_dir = results_dir + plate_name + '/'

    # Download
    download_s3(s3_uri, plate_reads_dir)

    for filepath in glob.glob(plate_reads_dir + '/*.fastq.gz'):
        rename_fastq_file(filepath)

    # Process
    run_pipeline(plate_reads_dir, plate_results_dir, plate_name, image=args.image)

    # TODO: Backup in fsx

if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="run pipeline on a routine Salmonella Plate")
    parser.add_argument("s3_uri", help="s3 uri that corresponds to the fastq plate to run")
    parser.add_argument("--reads-dir", default=DEFAULT_READS_DIRECTORY,  help="base directory that s3 objects are stored to")
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIRECTORY,  help="base directory where pipeline results are stored")
    parser.add_argument("--image", default=DEFAULT_IMAGE, help="docker image to use")

    args = parser.parse_args()

    # Run
    run_plate(args.s3_uri, args.reads_dir, args.results_dir)
