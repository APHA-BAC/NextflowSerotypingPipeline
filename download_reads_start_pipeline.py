#!/usr/bin/python3
# download_reads_start_pipeline.py by Josh 

import os
import shutil
import argparse
import subprocess
import sys

nextflowDir = os.path.expanduser("~/nextflow_new_ver_test")

# Get command line arguments.
def main():
    parser = argparse.ArgumentParser(description="Retrieve WGS data from the CSU S3 bucket and launch Salmonella typing pipeline.")
    parser.add_argument("bucket", help="The s3 key for the sequencing run, e.g. FZ2000/M02410_5275")
    args = parser.parse_args()
    return args

# Retrieve and display information about the WGS data to be retrieved
def check_WGS(s3Key):
    print("This is your chosen WGS location: {}'\n".format(s3Key))
    lsCommand = "aws s3 ls {}/".format(s3Key)
    contents = [x.decode("utf-8") for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
    releaseDate = contents[0].split()[0]
    outfmtDate = releaseDate.split("-")[2] + releaseDate.split("-")[1] + releaseDate.split("-")[0][2:4]
    readFiles = [x.split()[3] for x in contents]
    print("This location contains the following readfiles:\n")
    print(", ".join(readFiles) + "\n")
    print("Number of read-pairs (isolates):", len([x for x in readFiles if "_R1" in x]))
    print("This data was released on: {}\n".format(releaseDate))
    runDir = s3Key.split("/")[-1]
    outDir = "{}_APHA_{}".format(outfmtDate, runDir)
    return outDir, readFiles


def retrieve_from_bucket(s3Key, outDir):
    homeWGSDir = os.path.expanduser("~/WGS_Data/{}".format(outDir))
    if not os.path.isdir(homeWGSDir):
        os.makedirs(homeWGSDir)
    retrieveCommand = "aws s3 cp {} {} --recursive --include \"*.fastq.gz\"".format(s3Key, homeWGSDir)
    print(retrieveCommand)
    subprocess.call(retrieveCommand, shell=True)
    print("All readfiles downloaded.\n\n")
    return(homeWGSDir)


def rename_WGS(readFiles, homeWGSDir):
    for readFile in readFiles:
        if "_R1" in readFile:
            newName = readFile.split("_")[0] + "_R1.fastq.gz"
        elif "_R2" in readFile:
            newName = readFile.split("_")[0] + "_R2.fastq.gz"
        print("Renaming readfile {} --> {}".format(readFile, newName))
        os.rename("{}/{}".format(homeWGSDir, readFile), "{}/{}".format(homeWGSDir, newName))
    print("\n\n")


def start_pipeline_and_exit(outDir, nextflowDir):
    oldWork = os.path.join(nextflowDir, "work")
    if os.path.exists(oldWork):
        print("Removing old working directory:", oldWork)
        shutil.rmtree(oldWork)
    os.chdir(nextflowDir)
    nfCommand = "./nextflow SCE3_pipeline_update.nf --runID {}".format(outDir)
    print(nfCommand)
    subprocess.Popen(nfCommand, shell=True)
    sys.exit()

# Run
if __name__ == '__main__':
    args = main()
    print(args)
    s3Key = args.bucket
    s3Key = "s3://s3-csu-001/{}".format(s3Key)
    outDir, readFiles = check_WGS(s3Key)
    homeWGSDir = retrieve_from_bucket(s3Key, outDir)
    rename_WGS(readFiles, homeWGSDir)
    print("Setup complete. Please find downloaded and renamed readfiles in: {}\n\n".format(os.path.expanduser("~/WGS_Data/{}".format(outDir))))
    start_pipeline_and_exit(homeWGSDir, nextflowDir)


