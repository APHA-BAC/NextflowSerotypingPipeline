import subprocess

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

if __name__ == '__main__':
    results_directoy = '/home/aaronfishman/temp/salmonella-results/'
    reads_directory = "/home/aaronfishman/temp/salmonella-reads/"


    # Download raw reads from S3 --> EC2
    plate_name = "new-plate"
    plate_reads_dir = reads_directory + plate_name + '/'


    # TODO: Backup raw reads from EC2 --> FSX
    # Integrity checks should be done, Josh's script already does this

    plate_results_dir = results_directoy + plate_name + '/'

    print("plate_reads_dir", plate_reads_dir)
    print("plate_results_dir", plate_results_dir)

    run_pipeline(plate_reads_dir, plate_results_dir, plate_name)
