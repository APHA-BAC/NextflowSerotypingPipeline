#!/usr/bin/python3
# download_wgs_from_csu_bucket.py by Josh 

import os
import sys
import shutil
import argparse
import subprocess
import hashlib

# Check presence of mounted FSx archive disk
# Needs to be updated to mount the disk it if not mounted, rather than just exit!
def check_mount():
    if not os.path.isdir(os.path.expanduser("~/mnt/Salmonella/BAC3_NGS_Archive/Salmonella")):
        sys.exit("Error, FSx-Salmonella drive not mounted! Please mount first and try again.")

# This is for the user to interactively choose the WGS directory when the --dir option is not specified
def choose_dir():
    lsCommand = "aws s3 ls s3://s3-csu-001/"
    budgetCodes = [x.decode("utf-8").split()[-1][:-1] for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
    budgetCodeStr = ", ".join(budgetCodes)
    print("Here are the project/budget code sublocations  in 's3://s3-csu-001/' :\n")
    print(budgetCodeStr + "\n")
    budgetCode = input("Please input your chosen budget code, e.g. FZ2000: ")
    while budgetCode not in budgetCodes:
        print("Sorry, budget code not recognised. Please try again.")
        budgetCode = input("Input budget code: ")
    print("\nThanks. Checking data releases in this location; please wait...\n")
    lsCommand = "aws s3 ls s3://s3-csu-001/{}/".format(budgetCode)
    runCodes = [x.decode("utf-8").split()[-1] for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
    runCodes = [x[:-1] for x in runCodes if x[-1] == "/"]
    releases = []
    for runCode in runCodes:
        lsCommand = "aws s3 ls s3://s3-csu-001/{}/{}/".format(budgetCode, runCode)
        contents = [x.decode("utf-8") for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
        readSizes = [int(x.split()[2]) for x in contents]
        volume = str(round(sum(readSizes) / 1000000000, 2)) + " GB"
        releaseDate = contents[0].split()[0]
        outfmtDate = "/".join([releaseDate.split("-")[2], releaseDate.split("-")[1], releaseDate.split("-")[0]])
        releaseDate = int(releaseDate.replace("-", ""))
        releases.append([runCode, releaseDate, outfmtDate, volume])
    releases.sort(key=lambda x:x[1])
    print("Here are the sequencing runs/batches in 's3://s3-csu-001/{}/' :\n".format(budgetCode))
    for release in releases:
        print("\t\t".join([release[0], release[2], release[3]]))
    print("\n")
    runCode = input("Please input your chosen sequencing batch: ")
    while runCode not in runCodes:
        print("Sorry, sequencing run/batch not recognised. Please try again.")
        runCode = input("Input sequencing batch: ")
    print("\nThanks!\n")
    s3Path = "{}/{}".format(budgetCode, runCode)
    return(s3Path)

# Determine that each R1 file has a corresponing R2 and vice versa
def check_readpairs(readFiles):
    counterParts = [x.replace("_R1", "_R2") if "_R1" in x else x.replace("_R2", "_R1") for x in readFiles]
    if sorted(readFiles) != sorted(counterParts):
        unmatched = [x for x in readFiles if x not in counterParts]
        print("Readfiles without counterparts:", unmatched)
        return unmatched
    else:
        return None

# Retrieve and display information about the WGS data to be retrieved
def check_WGS(s3Path):
    print("This is your chosen run: {}\n".format(s3Path.split("/")[1]))
    lsCommand = "aws s3 ls s3://s3-csu-001/{}/".format(s3Path)
    contents = [x.decode("utf-8") for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
    contents = [x for x in contents if len(x.split()) == 4 and "fastq.gz" in x.split()[-1]]
    if len(contents) == 0:
        print("No reads found in this location.\n")
        return
    releaseDate = contents[0].split()[0]
    outfmtDate = "/".join([releaseDate.split("-")[2], releaseDate.split("-")[1], releaseDate.split("-")[0]])
    readSizes = [x.split()[2] for x in contents]
    readFiles = [x.split()[3] for x in contents]
    print("This run contains the following readfiles:\n")
    print(", ".join(readFiles) + "\n")
    unmatched = check_readpairs(readFiles)
    if unmatched:
        readSizes = [x for idx, x in enumerate(readSizes) if readFiles[idx] not in unmatched]
        readFiles = [x for x in readFiles if x not in unmatched]
    volume = str(round(sum([int(x) for x in readSizes]) / 1000000000, 2)) + " GB"
    print("Number of read-pairs (isolates):", len([x for x in readFiles if "_R1" in x]))
    print("This data was released on: {}".format(outfmtDate))
    print("Approx data vol: {}\n".format(volume))
    print("S3 uri: 's3://s3-csu-001/{}/'\n".format(s3Path))

# Confirmation user wants to proceed; interactive mode only
def confirm_WGS():
    confirm = input("Would you like to examine a different sequencing batch, Y or N?: ")
    confirm = confirm.upper()
    while confirm not in ("Y", "N"):
        print("Sorry, please answer Y or N.")
        confirm = input("Would you like to examine a different sequencing batch, Y or N?: ")
        confirm = confirm.upper()
    return confirm

# Get command line arguments. NB: Default is NO arguments = interactive mode!
def main():
    parser = argparse.ArgumentParser(description="Retrieve info about sequencing data releases in CSU S3 bucket. Run without --uri for default interactive mode.")
    parser.add_argument("--uri", required=False, help="OPTIONAL! The location in s3-csu-001 corresponding to the sequencing run, e.g. FZ2000/NB552234_0097")
    args = parser.parse_args()
    return args


# Run
confirm = "Y"

if __name__ == '__main__':
    check_mount()
    args = main()
    print(args)
    if args.uri:
        s3Path = args.uri
        check_WGS(s3Path)
    else:
        while confirm == "Y":
            s3Path = choose_dir()
            check_WGS(s3Path)
            confirm = confirm_WGS()

quit()