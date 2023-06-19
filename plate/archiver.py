#!/usr/bin/python3
#!/usr/bin/python3
# download_wgs_from_csu_bucket.py by Josh 

import os
import sys
import shutil
import argparse
import subprocess
import hashlib
import shutil


# Get command line arguments.
# --dir is optional; there is user interactivity if dir is not specified, but the intention will be
# to always call the script using the --dir option for automated running of the pipeline
def main():
    parser = argparse.ArgumentParser(description="Retrieve WGS data from the CSU S3 bucket and perform setup for the Salmonella typing pipeline.")
    # parser.add_argument("-b", "--budget", required=True, help="The project/budget code (and thus parent directory) under which the WGS data has been released. e.g. if data has been released in s3://s3-csu-001/FZ2000/n1234/, BUDGET = FZ2000")
    # parser.add_argument("-r", "--run", required=True, help="The directory name for the specific run, e.g., if data has been released in s3://s3-csu-001/FZ2000/n1234/, RUN = n1234")
    parser.add_argument("--uri", required=False, help="OPTIONAL! The location in s3-csu-001 corresponding to the sequencing run, e.g. FZ2000/NB552234_0097")
    args = parser.parse_args()
    return args


# Check presence of mounted FSx archive disk
# Needs to be updated to mount the disk it if not mounted, rather than just exit!
def check_mount():
    if not os.path.isdir(os.path.expanduser("~/mnt/Salmonella/BAC3_NGS_Archive/Salmonella")):
        sys.exit("Error, FSx-Salmonella drive not mounted! Please mount first and try again.")


# This is for the user to interactively choose the WGS directory when the --dir option is not specified
def choose_dir():
    lsCommand = "aws s3 ls s3://s3-csu-001/"
    budgetDirs = [x.decode("utf-8").split()[-1][:-1] for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
    budgetDirStr = ", ".join(budgetDirs)
    print("Here are the project/budget code directories in 's3://s3-csu-001/' :\n")
    print(budgetDirStr + "\n")
    budgetDir = input("Please input your chosen budget code, e.g. FZ2000: ")
    while budgetDir not in budgetDirs:
        print("Sorry, budget code directory not recognised. Please try again.")
        budgetDir = input("Input budget code: ")
    print("\nThanks. Checking data releases in this folder; please wait...\n")
    lsCommand = "aws s3 ls s3://s3-csu-001/{}/".format(budgetDir)
    runDirs = [x.decode("utf-8").split()[-1][:-1] for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
    releaseDates = []
    for runDir in runDirs:
        lsCommand = "aws s3 ls s3://s3-csu-001/{}/{}/".format(budgetDir, runDir)
        contents = [x.decode("utf-8") for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
        releaseDate = contents[0].split()[0]
        outfmtDate = "/".join([releaseDate.split("-")[2], releaseDate.split("-")[1], releaseDate.split("-")[0]])
        releaseDate = int(releaseDate.replace("-", ""))
        releaseDates.append([releaseDate, outfmtDate])
    runDirs = [[runDirs[j]] + releaseDates[j] for j in range(len(runDirs))]
    runDirs = [e for e in sorted(runDirs, key=lambda x:x[1])]
    print("Here are the sequencing runs/batches in 's3://s3-csu-001/{}/' :\n".format(budgetDir))
    for runDir in runDirs:
        print("\t\t".join([runDir[0], runDir[2]]))
    print("\n")
    runDir = input("Please input your chosen sequencing batch: ")
    while runDir not in [x[0] for x in runDirs]:
        print("Sorry, sequencing run/batch not recognised. Please try again.")
        runDir = input("Input sequencing batch: ")
    print("\nThanks.\n")
    s3uri = "{}/{}".format(budgetDir, runDir)
    return(s3uri)


# Confirmation user wants to proceed; interactive mode only
def confirm_WGS():
    confirm = input("Would you like to download these data (Y) or choose a different batch (N)?: ")
    while confirm not in ("Y", "N"):
        print("Sorry, please answer Y or N.")
        confirm = input("Would you like to download these data (Y) or choose a different batch (N)?: ")
    return confirm


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
def check_WGS(s3uri):
    print("This is your chosen s3 location: 's3://s3-csu-001/{}/'\n".format(s3uri))
    lsCommand = "aws s3 ls s3://s3-csu-001/{}/".format(s3uri)
    contents = [x.decode("utf-8") for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
    releaseDate = contents[0].split()[0]
    outfmtDate = releaseDate.split("-")[2] + releaseDate.split("-")[1] + releaseDate.split("-")[0][2:4]
    readSizes = [x.split()[2] for x in contents]
    readFiles = [x.split()[3] for x in contents]
    print("This directory contains the following readfiles:\n")
    print(", ".join(readFiles) + "\n")
    unmatched = check_readpairs(readFiles)
    if unmatched:
        readSizes = [x for idx, x in enumerate(readSizes) if readFiles[idx] not in unmatched]
        readFiles = [x for x in readFiles if x not in unmatched]
    print("Number of read-pairs (isolates):", len([x for x in readFiles if "_R1" in x]))
    print("This data was released on: {}\n".format(releaseDate))
    runDir = s3uri.split("/")[1]
    outDir = "{}_APHA_{}".format(outfmtDate, runDir)
    return outDir, readFiles, readSizes


def check_retrieval(homeWGSDir, readFile, readSize):
    if os.path.isfile("{}/{}".format(homeWGSDir, readFile)) and str(os.path.getsize("{}/{}".format(homeWGSDir, readFile))) == readSize:
        return True
    else:
        return False

def retrieve_from_bucket(s3uri, outDir, readFiles, readSizes):
    homeWGSDir = os.path.expanduser("~/wgs-reads/{}".format(outDir))
    if not os.path.isdir(homeWGSDir):
        os.makedirs(homeWGSDir)
    retrieveCommand = "aws s3 cp s3://s3-csu-001/{} {} --recursive --include \"*.fastq.gz\"".format(s3uri, homeWGSDir)
    print(retrieveCommand)
    subprocess.call(retrieveCommand, shell=True)
    for idx, readFile in enumerate(readFiles):
        n = 0
        while not check_retrieval(homeWGSDir, readFile, readSizes[idx]):
            if n > 4:
                print(readFile)
                sys.exit("Error, readfile not retrieved correctly from CSU S3 bucket after 5 tries!")
                break
            subprocess.call("aws s3 cp s3://s3-csu-001/{}/{} {}".format(s3uri, readFile, homeWGSDir), shell=True)
            n += 1
    print("All readfiles successfully downloaded.")
    print("Archiving now...\n")
    return(homeWGSDir)

def get_md5(readFile):
    # The 'b' in the open clause is very important!
    # This allows open() to read the gzipped file correctly and generate the md5 of the file
    with open(readFile, "rb") as f:
        readData = f.read()
        readMD5 = hashlib.md5(readData).hexdigest()
    readData = None
    return(readMD5)

def check_archival(destFile, origMD5):
    if os.path.isfile(destFile) and get_md5(destFile) == origMD5:
        return True
    else:
        return False

# Copy data to FSx disk
def archive_WGS(outDir, readFiles, homeWGSDir):
    archiveDir = os.path.expanduser("~/mnt/Salmonella/BAC3_NGS_Archive/Salmonella")
    archiveDir = "{}/{}".format(archiveDir, outDir)
    if not os.path.isdir(archiveDir):
        os.makedirs(archiveDir)
    for idx, readFile in enumerate(readFiles):
        origMD5 = get_md5("{}/{}".format(homeWGSDir, readFile))
        destFile = "{}/{}".format(archiveDir, readFile)
        if check_archival(destFile, origMD5):
            print("Readfile",readFile,"already in archive.")
        else:
            n = 0
            while not check_archival(destFile, origMD5):
                if n > 4:
                    sys.exit("Error, readfile not copied correctly to archive after 5 tries!")
                    break
                print("Copying readfile", readFile, "to archive now...")
                shutil.copy("{}/{}".format(homeWGSDir, readFile), archiveDir)
                n += 1
            print("Done.")
    print("All readfiles copied to archive.\n\n")


def rename_WGS(readFiles, homeWGSDir):
    for readFile in readFiles:
        if "_R1" in readFile:
            newName = readFile.split("_")[0] + "_R1.fastq.gz"
        elif "_R2" in readFile:
            newName = readFile.split("_")[0] + "_R2.fastq.gz"
        else:
            print(readFile)
            sys.exit("Error, read file name not in recognised format!")
        print("Renaming readfile {} --> {}".format(readFile, newName))
        os.rename("{}/{}".format(homeWGSDir, readFile), "{}/{}".format(homeWGSDir, newName))
    print("\n\n")



confirm = "N"

# Run
# if __name__ == '__main__':
#     check_mount()
#     args = main()
#     print(args)
#     if args.uri:
#         s3uri = args.uri
#         outDir, readFiles, readSizes = check_WGS(s3uri)
#     else:
#         while confirm == "N":
#             s3uri = choose_dir()
#             outDir, readFiles, readSizes = check_WGS(s3uri)
#             confirm = confirm_WGS()
#     homeWGSDir = retrieve_from_bucket(s3uri, outDir, readFiles, readSizes)
#     archive_WGS(outDir, readFiles, homeWGSDir)
#     shutil.rmtree(homeWGSDir)
#     print("All finished.")

# quit()




# def check_s3uri(s3uri):
#     budgetDir = s3uri.split("/")[0]
#     runDir = s3uri.split("/")[1]
#     lsCommand = "aws s3 ls s3://s3-csu-001/"
#     budgetDirs = [x.decode("utf-8").split()[-1][:-1] for x in subprocess.check_output(lsCommand, shell=True).splitlines()]
