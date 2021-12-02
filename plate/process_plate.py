import subprocess


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
    return ""

if __name__ == '__main__':
    # TODO: Use Josh's interactive plate selection code

    s3_uri = 's3://s3-csu-001/temp/new-plate/'

    plate_name = <date ddmmyy>_APHA_<plate_name>

    quit()

    download_s3(s3_uri, )


    # Download raw reads from S3 --> EC2

   
    plate_name = "new-plate"

    # Run plate
    plate_reads_dir = READS_DIRECTORY + plate_name + '/'

    # TODO: Backup raw reads from EC2 --> FSX
    # Integrity checks should be done, Josh's script already does this

    plate_results_dir = RESULTS_DIRECTORY + plate_name + '/'

    # Run the pipeline

    run_pipeline(plate_reads_dir, plate_results_dir, plate_name)

    # TODO: Backup to fsx
    