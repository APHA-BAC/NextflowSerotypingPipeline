import subprocess
import os
import glob
import argparse
import boto3
from archiver import *

# TODO: Rename directories to BGE defaults
DEFAULT_READS_DIRECTORY = os.path.expanduser('/root/WGS_Data')
DEFAULT_RESULTS_DIRECTORY = os.path.expanduser('/root/WGS_Results')
DEFAULT_IMAGE = "jguzinski/salmonella-seq:prod"
DEFAULT_KMERID_REF = os.path.expanduser('/root//KmerID_Ref_Genomes/ref/')
DEFAULT_KMERID_CONFIG = os.path.expanduser('/root/KmerID_Ref_Genomes/config/')
s3_destination = "s3://s3-staging-area/arslanhussaini/"

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

def run_pipeline(plate_name):
    """ Run the Salmonella pipeline using docker """
    
    run(["/root/nextflow/nextflow", "SCE3_pipeline_update.nf",
        "--runID", plate_name])

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


def upload_s3(summaryTable_path, s3_destination):
    """ Upload summary table to S3 bucket """
    print("****** " + summaryTable_path + " ******")
    print("****** " + s3_destination + " ******")
    try:
        run(["aws", "s3", "cp",summaryTable_path,s3_destination])
    except:
        print("Does the destination path exist?")

def download_kmerid():
    run(["aws", "s3", "cp", "--acl", "bucket-owner-full-control", "--recursive", "s3://s3-ranch-046/KmerID_Ref_Genomes", "/root/KmerID_Ref_Genomes/"])
     

def run_plate(s3_uri, reads_dir, results_dir, local, runID, upload, transfer):

    """ Download, process and store a plate of raw Salmonella data """
    download_kmerid()

    # Add trailing slash to directory names
    reads_dir = os.path.join(reads_dir, '')
    results_dir = os.path.join(results_dir, '')
    plate_reads_dir = ''
    plate_results_dir = ''
    plate_name = ''

    if local == 0:
        # Storage paths
        plate_name = s3_uri_to_plate_name(s3_uri)
        plate_reads_dir = reads_dir + plate_name + '/'
        plate_results_dir = results_dir + plate_name + '/'

        # Download
        download_s3(s3_uri, plate_reads_dir)

        for filepath in glob.glob(plate_reads_dir + '/*.fastq.gz'):
            rename_fastq_file(filepath)

    elif local == 1:

        plate_name = runID
        plate_reads_dir = reads_dir + runID + '/'
        plate_results_dir = results_dir + runID +'/'

        for filepath in glob.glob(plate_reads_dir + '/*.fastq.gz'):
            rename_fastq_file(filepath)
    # Process
    run_pipeline(plate_name)
    if upload == 1:
        for file in glob.glob(results_dir + runID + "/" + r'*plusLIMS.csv'):
            try:
                shutil.copy(file, '~/mnt/my_share/CR2009')
                print("Summary table copied")
            except:
                print("Copy of summary table failed. Is the drive mounted?")

    # If running plate from s3_uri, backup
    # try:
    #     if local == 0:
    #         new_s3_uri = s3_uri[16:-1]            
    #         check_mount()
    #         outDir, readFiles, readSizes = check_WGS(new_s3_uri)
    #         homeWGSDir = retrieve_from_bucket(new_s3_uri, outDir, readFiles, readSizes)
    #         archive_WGS(outDir, readFiles, homeWGSDir)
    #         shutil.rmtree(homeWGSDir)
    # except:
    #     print("Archive failed")

    if transfer:
        # Sets up the string that is the path to the summary table
        TableFile = plate_name + "_SummaryTable_plusLIMS.csv"
        summaryTable_path = os.path.join("~/root/WGS_Results/",plate_name,TableFile)
        summaryTable_path = os.path.expanduser(summaryTable_path)
        upload_s3(summaryTable_path,transfer)



if __name__ == '__main__':
    # Parse
    
    parser = argparse.ArgumentParser(description="run pipeline on a routine Salmonella Plate")
    parser.add_argument("-s","--s3_uri", help="s3 uri that corresponds to the fastq plate to run")
    parser.add_argument("--reads-dir", default=DEFAULT_READS_DIRECTORY,  help="base directory that s3 objects are stored to")
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIRECTORY,  help="base directory where pipeline results are stored")
    parser.add_argument("--image", default=DEFAULT_IMAGE, help="docker image to use")
    parser.add_argument("-l","--local", default=0, help="Set to 1 if your reads are in a local directory. Default is 0")
    parser.add_argument("-r","--runID", help="The name of the run which will also be the name of the directory for the results. Only needed if running locally")
    parser.add_argument("-u", "--upload", default=0, help="Set to 1 if you want to upload to SMB staging area")
    parser.add_argument("-t", "--transfer", default=False, help="Set to to 1 to transfer to S3 bucket")


    args = parser.parse_args()
    args.local = int(args.local)
    if args.runID == None and args.local == 1:
         parser.error("--local requires --runID")

    # Run

    run_plate(args.s3_uri, args.reads_dir, args.results_dir, args.local, args.runID, args.upload, args.transfer)

