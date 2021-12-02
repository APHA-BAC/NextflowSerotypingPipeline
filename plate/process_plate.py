import subprocess
import os
import glob
import argparse

RESULTS_DIRECTORY = '/home/aaronfishman/temp/salmonella-results/'
READS_DIRECTORY = "/home/aaronfishman/temp/salmonella-reads/" 

# TODO: get plate name from Josh's cpde
def run_pipeline(reads, results, plate_name):
    run([
        "sudo", "docker", "run", "-it", 
        "-v", f"{reads}:/root/WGS_Data/{plate_name}/",
        "-v", f"{results}:/root/WGS_Results/{plate_name}/", 
        "jguzinski/salmonella-seq:prod",
        "/root/nextflow/nextflow", "SCE3_pipeline_update.nf",
        "--runID", plate_name
    ])

def run(cmd):
    """ Run a command and assert that the process exits with a non-zero exit code.
        See python's subprocess.run command for args/kwargs
        Parameters:
            cmd (list): List of strings defining the command, see (subprocess.run in python docs)
            cwd (str): Set surr
        Returns:
            None
    """
    # TODO: store stdout to a file
    returncode = subprocess.run(cmd).returncode

    if returncode:
        raise Exception("""*****
            %s
            cmd failed with exit code %i
        *****""" % (cmd, returncode))

def rename_WGS(readFiles, homeWGSDir):
    for readFile in readFiles:
        if "_R1" in readFile:
            newName = readFile.split("_")[0] + "_R1.fastq.gz"
        elif "_R2" in readFile:
            newName = readFile.split("_")[0] + "_R2.fastq.gz"
        print("Renaming readfile {} --> {}".format(readFile, newName))
        os.rename("{}/{}".format(homeWGSDir, readFile), "{}/{}".format(homeWGSDir, newName))
    print("\n\n")

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


def run_plate(s3_uri):
    """ Download, process and store a plate of raw Salmonella data """

    # 
    plate_name = s3_uri_to_plate_name(s3_uri.strip('/'))
    plate_reads_dir = READS_DIRECTORY + plate_name + '/'
    plate_results_dir = RESULTS_DIRECTORY + plate_name + '/'

    # Download
    download_s3(s3_uri, plate_reads_dir)

    for filepath in glob.glob(plate_reads_dir + '/*.fastq.gz'):
        rename_fastq_file(filepath)

    # Run
    run_pipeline(plate_reads_dir, plate_results_dir, plate_name)

    # TODO: Backup in fsx

if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Run pipeline on a routine Salmonella Plate")
    parser.add_argument("s3_uri", help="s3 uri to the plate you want to run")
    
    args = parser.parse_args()

    # Run
    run_plate(args.s3_uri)
