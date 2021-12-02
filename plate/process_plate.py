import subprocess
import os
import glob
import argparse

RESULTS_DIRECTORY = '/home/aaronfishman/temp/salmonella-results/'
READS_DIRECTORY = "/home/aaronfishman/temp/salmonella-reads/" 

# TODO: get plate name from Josh's cpde
def run_pipeline(reads, results, plate_name):
    subprocess.run([
        "sudo", "docker", "run", "-it", 
        "-v", f"{reads}:/root/WGS_Data/{plate_name}/",
        "-v", f"{results}:/root/WGS_Results/{plate_name}/", 
        "jguzinski/salmonella-seq:prod",
        "/root/nextflow/nextflow", "SCE3_pipeline_update.nf",
        "--runID", plate_name
    ])



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
    subprocess.run(["aws", "s3", "cp", "--recursive",
        s3_uri,
        destination
    ])

    # TODO: throw exception if it fails


def s3_uri_to_plate_name(s3_uri):
    # <date ddmmyy>_APHA_<plate_name>
    run_name = os.path.basename(s3_uri.strip('/'))

    return run_name

    # Retrieve and display information about the WGS data to be retrieved
def check_WGS(s3Key):
    s3Key = s3Key.strip('/')

    # TODO: Bit of a tidy, not that bad though really
    # print("This is your chosen WGS location: {}'\n".format(s3Key))
    lsCommand = "aws s3 ls {}/".format(s3Key)
    contents = [x.decode("utf-8") for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
    releaseDate = contents[0].split()[0]
    outfmtDate = releaseDate.split("-")[2] + releaseDate.split("-")[1] + releaseDate.split("-")[0][2:4]
    readFiles = [x.split()[3] for x in contents]
    # print("This location contains the following readfiles:\n")
    # print(", ".join(readFiles) + "\n")
    # print("Number of read-pairs (isolates):", len([x for x in readFiles if "_R1" in x]))
    # print("This data was released on: {}\n".format(releaseDate))
    runDir = s3Key.split("/")[-1]
    outDir = "{}_APHA_{}".format(outfmtDate, runDir)
    return outDir


def rename_WGS(readFiles, homeWGSDir):
    for readFile in [os.path.basename(x) for x in readFiles]:
        if "_R1" in readFile:
            newName = readFile.split("_")[0] + "_R1.fastq.gz"
        elif "_R2" in readFile:
            newName = readFile.split("_")[0] + "_R2.fastq.gz"
        print("Renaming readfile {} --> {}".format(readFile, newName))
        os.rename("{}/{}".format(homeWGSDir, readFile), "{}/{}".format(homeWGSDir, newName))
    print("\n\n")


def run_plate(s3_uri):
    plate_name = check_WGS(s3_uri)
    plate_reads_dir = READS_DIRECTORY + plate_name + '/'
    plate_results_dir = RESULTS_DIRECTORY + plate_name + '/'

    download_s3(s3_uri, plate_reads_dir)
    rename_WGS(glob.glob(plate_reads_dir+'/*.fastq.gz'), plate_reads_dir)

    run_pipeline(plate_reads_dir, plate_results_dir, plate_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run pipeline on a routine Salmonella Plate")
    parser.add_argument("s3_uri", help="s3 uri to the plate you want to run")
    
    args = parser.parse_args()

    run_plate(args.s3_uri)
