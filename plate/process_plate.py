import subprocess
import os
import glob
import argparse
import boto3
import update_master_sum

# TODO: Rename directories to BGE defaults
DEFAULT_READS_DIRECTORY = os.path.expanduser('~/wgs-reads/')
DEFAULT_RESULTS_DIRECTORY = os.path.expanduser('~/wgs-results/')
# DEFAULT_IMAGE = "ahussaini96/serotypingpipeline:process_cleanup"
DEFAULT_IMAGE = "jguzinski/salmonella-seq:prod"
DEFAULT_KMERID = os.path.expanduser('s3://s3-ranch-046/KmerID_Ref_Genomes/')
DEFAULT_KMERID_REF = os.path.expanduser(str(DEFAULT_KMERID) + "ref/")
DEFAULT_KMERID_CONFIG = os.path.expanduser(str(DEFAULT_KMERID) + "config/")
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

def run_pipeline(reads, results, plate_name, assemblies_dir, image=DEFAULT_IMAGE, kmerid_ref=DEFAULT_KMERID_REF, kmerid_config=DEFAULT_KMERID_CONFIG):
    """ Run the Salmonella pipeline using docker """
    if not os.path.isdir("/opt/kmerid/"):
        run(['sudo', 'aws', 's3', 'cp', '--recursive', 's3://s3-ranch-046/KmerID_Ref_Genomes/ref', '/opt/kmerid/ref'])
        run(['sudo', 'aws', 's3', 'cp', '--recursive', 's3://s3-ranch-046/KmerID_Ref_Genomes/config', '/opt/kmerid/config'])
    run(["sudo", "docker", "pull", image])
    run([
        "sudo", "docker", "run", "--rm", "-it",
        "-v", f"{reads}:/root/wgs-reads/{plate_name}/",
        "-v", f"{results}:/root/wgs-results/{plate_name}/",
        "-v", f"{assemblies_dir}:/root/wgs-results/{plate_name}/assemblies/",
        "-v", f"/opt/kmerid/ref:/opt/kmerid/ref",
        "-v", f"/opt/kmerid/config:/opt/kmerid/config",
        image,
        "/root/nextflow/nextflow", "SCE3_pipeline_update.nf",
        "--runID", plate_name, "--plateRun", "False"
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


def upload_s3(summaryTable_path, s3_destination):
    """ Upload summary table to S3 bucket """
    try:
        run(["aws", "s3", "cp",summaryTable_path,s3_destination])
        print("Summary table uploaded to: {}".format(s3_destination))
    except:
        print("Does the destination path exist?")

def run_plate(s3_uri, reads_dir, results_dir, image, runID, transfer, updateSum, upload_fasta):

    """ Download, process and store a plate of raw Salmonella data """

    # Add trailing slash to directory names
    reads_dir = os.path.join(reads_dir, '')
    results_dir = os.path.join(results_dir, '')
    
    plate_reads_dir = ''
    plate_results_dir = ''
    plate_name = ''

    if not runID:
        # Storage paths
        plate_name = s3_uri_to_plate_name(s3_uri)
        plate_reads_dir = reads_dir + plate_name + '/'
        plate_results_dir = results_dir + plate_name + '/'

        # Download
        download_s3(s3_uri, plate_reads_dir)

        for filepath in glob.glob(plate_reads_dir + '/*.fastq.gz'):
            rename_fastq_file(filepath)

    elif runID:

        plate_name = runID
        plate_reads_dir = reads_dir + runID + '/'
        plate_results_dir = results_dir + runID +'/'

        for filepath in glob.glob(plate_reads_dir + '/*.fastq.gz'):
            rename_fastq_file(filepath)
    # run(["mkdir", assemblies_dir])
    assemblies_dir = str(plate_results_dir) + "assemblies/"
    run_pipeline(plate_reads_dir, plate_results_dir, plate_name, assemblies_dir, image=args.image)
    TableFile = plate_name + "_SummaryTable_plusLIMS.csv"
    summaryTable_path = os.path.join("~/wgs-results/",plate_name,TableFile)
    summaryTable_path = os.path.expanduser(summaryTable_path)

    if updateSum:
        try:
            update_master_sum.update_summary(updateSum, summaryTable_path)
            print("Master summary table updated")
        except:
            print("Unable to update master summary table")


    # Transfer to S3 if transfer flag set
    if transfer:
        # Sets up the string that is the path to the summary table
        TableFile = plate_name + "_SummaryTable_plusLIMS.csv"
        summaryTable_path = os.path.join("~/wgs-results/",plate_name,TableFile)
        summaryTable_path = os.path.expanduser(summaryTable_path)
        upload_s3(summaryTable_path,transfer)

    if upload_fasta:
        
        run(["aws", "s3", "cp", "--recursive", assemblies_dir, "s3://s3-ranch-050/assemblies", "--acl", "bucket-owner-full-control"])
        


if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="run pipeline on a routine Salmonella Plate")
    parser.add_argument("-s","--s3_uri", help="s3 uri that corresponds to the fastq plate to run")
    parser.add_argument("--reads-dir", default=DEFAULT_READS_DIRECTORY,  help="base directory that s3 objects are stored to")
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIRECTORY,  help="base directory where pipeline results are stored")
    parser.add_argument("--image", default=DEFAULT_IMAGE, help="docker image to use")
    parser.add_argument("-r","--runID", default=False, help="The name of the run which should be the name of the folder with your reads")
    parser.add_argument("-t", "--transfer", default=0, help="Set to the path for your desired s3 destination bucket")
    parser.add_argument("-u", "--updateSum", default=False, help="Set to path of master summary table if you want to update the master the summary table with the results from this run")
    parser.add_argument("-k", "--kmerID", default = DEFAULT_KMERID, help="Set path to KMER ID genome files")
    parser.add_argument("-f", "--fasta", default=False, help="Upload fasta files to S3 bucket. Default is False")

    args = parser.parse_args()
    
    # Run
    run_plate(args.s3_uri, args.reads_dir, args.results_dir, args.image, args.runID, args.transfer, args.updateSum, args.fasta)

