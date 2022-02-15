[![Build Status](https://apha.teamcity.com/app/rest/builds/buildType:(id:NextflowSerotypingPipeline_Pipeline)/statusIcon)](https://apha.teamcity.com/viewType.html?buildTypeId=NextflowSerotypingPipeline_Pipeline&guest=1)

# Introduction

Salmonella whole genome sequencing (WGS) serotyping pipeline developed by APHA and written in [Nextflow](https://www.nextflow.io/). The pipeline compares outputs from several publicly available serotyping tools to increase performance.

# Installation
To install the nextflow salmonella serotyping pipeline:
```
  $ cd install
  $ bash install.bash
```

This script installs the following dependancies:
- `conda`
- `fastp`
- `FastQC`
- `KmerID`
- `MOST`
- `Nextflow`
- `Quast`
- `SeqSero2`
- `Seqtk`
- `Shovill`
- `sistr`
- `SRST2`
- `sratoolkit`

# Running the pipeline

To run the pipeline on a batch a samples, the raw `.fastq.gz` files must be stored in `~/WGS_Data/<runID>`.  Each read-pair sample is represented by a pair of files named `*_R1.fastq.gz` and `*_R2.fastq.gz`. For example, to batch two samples named `salmonella_a` and `salmonella_b`, a `~/WGS_Data/<runID>` directory containing four files is required: `salmonella_a_R1.fastq.gz`, `salmonella_a_R2.fastq.gz`,  `salmonella_b_R1.fastq.gz` and `salmonella_b_R2.fastq.gz`, needs to be defined.

Then, to run the pipeline from the terminal call:
```
  $ nextflow run SCE3_pipeline_update.nf --runID <runID> 
```

Pipeline output is stored in  `~/WGS_Results/<runID>/` and contains:
- TODO

## Run from docker

**Note:** While running from the terminal is the easiest method for developers and data analysts, the pipeline can also be run from docker. This method has the benefit of working across platforms while guaranteeing consistency with automated tests (see below). 

A docker image containing all required dependancies is provided [here](https://hub.docker.com/r/jguzinski/salmonella-seq). 

This pull the latest image (if it's not already fetched) from dockerhub and run the container on data
```
sudo docker run jguzinksi/salmonella-seq:prod
sudo docker run --rm -it -v /ABS/PATH/TO/READS/:/reads/ -v /ABS/PATH/TO/RESULTS/:/results/ jguzinksi/salmonella-seq:prod /root/nextflow/nextflow SCE3_pipeline_update.nf --runID <runID>
```

# Pipeline Algorithm 

The pipelines processes data in three stages, as shown below. During the preprocessing stage; low quality bases, adapter sequences are removed from the fastq sample file and then subsampled to a maximum of 3M reads. Following this, the analysis stage runs multiple serotyping tools in parallel.
This strategy has been demonstrated to provide more accurate serovar detection than any tools running individually.
Bespoke typing of vaccine strains and particular servoars of interest to APHA is included as part of this analysis. 
Outputs from each tool are compared in the consensus call step. 
The final postprocessing stage assigns an `Outcome` to each sample by analysing data gathered during the analysis stage. The following `Outcome`s are used to signify subsequent lab processing steps:

- **Pass**: The sample contains Salmonella and has a serovar assigned
- **Check Required**: Further scrutiny of the data is required
- **Inconclusive**: The sample contains insufficient data volumes for analysis, is contaminated, or has low assembly quality. 

![image](https://user-images.githubusercontent.com/6979169/154100120-0199a72f-aec6-482f-9dc0-5ddd38c13c3c.png)

## Validation

This pipeline has been internally validated and approved by the APHA validation team against a large dataset in excess of 3,500 samples as well as non-Salmonella isolates to determine sensitivity and specificity.


## Automated Tests

The automated tests provided here ensure the software runs as expected. If you make changes to the algorithm, it is **strongly** reccomended that you run these tests to verify the pipeline is behaving as intended. The tests are also automatically run by [TEAMCITY](https://apha.teamcity.com/viewType.html?buildTypeId=NextflowSerotypingPipeline_Pipeline&guest=1) on each pull-request. 

The automated tests assert the correct serovar is assigned to known samples. These are called `inclusivity` tests. To run an inclusivity test, call:
```
$ bash -e test_isolates/jobs/inclusivity.bash 0
```

# Release Process

Release a new version of the software, the master branch needs only to be merged into `prod` branch. To perform this merge, a pull-request from the `master` branch into the prod branch needs to be made. Approval of pull-requests to `prod` is made by the CODEOWNER (Liljana Petrovska). The CODEOWNER is responsible for ensuring the code conforms to the reliability tests defined in BAC 0429 Section 4. A positive test result is required for approval.

To release a new version of the software:
1. A developer makes a pull-request from the `master` to `prod` branch. The CODEOWNER is automatically notified by e-mail.
1. The CODEOWNER runs the reliability tests defined in BAC 0429 Section 4 and reviews the code changes. 
1. The CODEOWNER approves the pull-request if they satisfied, or requests changes.
1. The dev tags the current of head of master as the next version. Versions are numbered incrementally with numbers, for example `v1`, `v2`, etc. This can be performed by navigating to github master branch and selecting `Create a release`
1. The dev merges the `master` branch into `prod`




## Utilization
1) Newly sequenced samples in paired FASTQ will be deposited by the APHA Sequencing Unit (SCU) into the SCU Amazon Simple Storage Service (Amazon S3) bucket. All of the sequencing files belonging to a particular sequencing run will be stored in a zipped folder with a unique name corresponding only to that sequencing run (for example 211220_APHA_run_n1041). The user will need to download the zipped folder via AWS Management Console interface onto the VM hosting the pipeline and unzip the folder using the unzip command. The unzipped folder then needs to be copied to the /home/$USER/WGS_Data directory.

2) Prior to running the pipeline it is absolutely imperative to edit the names of the newly sequenced files. This procedure will remove sections of the file name that were appended by the Illumina sequencing software which are redundant and will impede the pipeline from correctly recognizing that it has to process these sequencing files. Crucially, the file name renaming procedure will leave the unique isolate identifier unchanged in each of the file names. The renaming of the files is done using the changeAPHArunName.sh script: 
-	Copy the changeAPHArunName.sh script from /home/$USER to /home/$USER/WGS_Data/211220_APHA_run_n1041 (the folder name with the sequenced files in this example is 211220_APHA_run_n1041 but this will be different for every sequencing run and thus needs to be edited accordingly).
-	Open the terminal and access the /home/$USER/WGS_Data/211220_APHA_run_n1041 directory using the cd /home/$USER/WGS_Data/211220_APHA_run_n1041 command.
-	Execute the bash changeAPHArunName.sh command to change the names of all of the sequence files stored in the /home/$USER/WGS_Data/211220_APHA_run_n1041 directory.
-	For example, S02387-19_S193_R1_001.fastq.gz will be changed to S02387-19_R1.fastq.gz 
-	The changeAPHArunName.sh should be deleted from the /home/$USER/WGS_Data/211220_APHA_run_n1041 directory. 

3) In the terminal, use the command line to launch a pipeline run by typing:
cd /home/$USER/nextflow
./nextflow run SCE3_pipeline_update.nf

When prompted by a message in the terminal window the user should enter the name of the folder with the sequence files (from the above example this would be: “211220_APHA_run_n1041”), and press the Enter key. This will launch the pipeline.

In the SCE3_pipeline_update.nf script, values for some global variables are assigned. These include the paths to the directories of the tools used during the analysis namely: FastQC, KmerID, SeqSero2, MOST, Shovill, Quast, SISTR and SRST2. Additionally, there are paths to the WGS_Data and WGS_Results directories which house the pipeline input and output files, respectively. These computer paths can be updated if circumstances change; for example if a new version of KmerID is installed in a different directory. 

Specific tasks performed by the SCE3_pipeline_update.nf script include: 

Operating within the Nextflow framework, the pipeline will analyse all of the isolates with multiple software in parallel, utilising all available CPU cores that the VM possesses. Therefore, multiple isolates will be analysed with the same software concurrently and multiple software will be analysing the same isolate concurrently, with the exceptions specified below. The order as to which isolates will be analysed with which software and which software will be analysing which isolate is totally random except for the specified exceptions. Crucially, the unique isolate identifiers specified in the sequencing file names are consistently applied and utilized by the pipeline script throughout the running of the pipeline thus mitigating any risk of a mix-up between any of the analysed isolates.   

•	Isolates will be analysed concurrently by FastQC, shovill, KmerID, SeqSero2, and MOST.  
•	Analysis of an isolate with Quast will only commence upon successful completion of the genome assembly with shovill for that isolate.
•	Analysis of an isolate with SISTR will only commence upon successful completion of the genome assembly with shovill for that isolate.
•	Analysis of an isolate with SRST2 will only commence upon successful completion of serotyping with MOST for that isolate and only if the isolate’s serotype was identified as either Typhimurium or Enteritidis or Gallinarum or Pullorum or Idikan or Kedougou or Java or Paratyphi.

For each of the pipeline component software, the results are stored under the directory /home/$USER/WGS_Results/<runName>/<isolateName>/<toolName>. The results are written out “live” as the pipeline progresses.

When all the samples have been processed, the python script summaryTable_Python3.py is called automatically. This script extracts and combines and presents the information for all the samples and individual pipeline component software results into a CSV table which is named runName_SummaryTable.csv stored at location /home/$USER/WGS_Results/<runName>.

4) A message in the style of “Completed at: 25-Dec-2020 18:23:48” will appear in the terminal window upon successful completion of a pipeline run. 
In the directory /home/$USER/WGS_Results there will be a folder with the same name as was entered into the terminal window upon launching of the pipeline (such as “211220_APHA_run_n1041”) and hence the same name as the name of the folder that was copied into /home/$USER/WGS_Data. This folder in /home/$USER/WGS_Results will contain a series of sub-folders, each of which will be named after the unique isolate identifier for each of the files that were sequenced in the sequencing run. Each sub-folder will contain pipeline outputs specific only to the isolate that the sub-folder was named after. The folder in /home/$USER/WGS_Results will also contain a summary table CSV document named 211220_APHA_run_n1041_SummaryTable.csv (the “211220_APHA_run_n1041” prefix will be replaced by the name of the current sequencing run) that summarizes the outputs of a pipeline run.

The rows in the summary table relate to the strains in the sequencing run and the columns to the relevant information extracted from the results of the analysis tools; these are the following from left to right:
1.	StrainID. Name or identifier of the strain provided by the FASTQ files
2.	Consensus. A score showing the number of programs in agreement between, MOST, SeqSero and SISTR. Scores 1-3 stand for: 1 – only one program identified the sample, 2 – two programs are in agreement on the sample’s identity and one is not, and 3 – all three programs are in agreement on the identity of the sample 
3.	#ReadsR1. Number of DNA fragments that the R1 FASTQ file contains.
4.	GC%R1. Percentage of the GC content of the fragments of the R1 FASTQ.
5.	R1Kmerid. Percentage of similarity of the sample according to the KmerID analysis tool, using the R1 FASTQ, to the two most similar organisms.
6.	Contamination flag  
7.	MOST. Serotype according to MOST
8.	Most_light. The traffic light Quality Control system: The script validates the results based on coverage metrics and writes a cut-off value standard based on the "Traffic light system".
-	The "GREEN traffic light" indication is assigned if the max percentage non consensus depth < 15% and Complete pileup= TRUE and Minimum consensus depth > 2 and Percentage coverage =100 and ST not "Failed(incomplete locus coverage)
-	The "AMBER traffic light" indication is assigned if the there is no exact fit which matches either GREEN or RED
-	The RED traffic light indication is assigned if the Complete pileup= FAIL or Percentage coverage < 100 or ST is "Failed(incomplete locus coverage)
9.	st is the sequence type otherwise known as the allelic profile 
10.	MLST.
11.	MLST mean cov. Mean coverage for the 7 core MLST genes 
12.	Seqsero. Serotype according to Seqsero
13.	SS comment. Additional comment for Seqsero tool output.
14.	N50
15.	Serogroup as defined by SISTR
16.	Serovar as defined by SISTR
17.	serovar_antigen as defined by SISTR
18.	serovar_cgmlst (core genome MLST) as defined by SISTR
19.	vaccine. To differentiate vaccine from field S. Typhimurium, S. Enteritidis, S. Gallinarum and S. Pulorrum strains.
20.	mono. Differentiation of monophasic S. 13,23:i:- from biphasic S. Kedougou (1,13,23:i:l,w) and S. Idikan (1,13,23:i:1,5) isolates
21.	sseJ. Differentiation of S. Paratyphi B from S. Paratyphi B Java variants based on the detection of sseJ gene, present only in the genomes of S. Paratyphi B Java isolates.

5) Prior to any subsequent pipeline runs, the “211220_APHA_run_n1041” folder will need to be deleted from /home/$USER/WGS_Data and from /home/$USER/WGS_Results. Moreover, a folder called “work” that collates Nextflow temporary files during a pipeline run will need to be deleted from the /home/$USER/nextflow directory. These operations will not be performed automatically and need to performed by the user.

  




![image](https://user-images.githubusercontent.com/6979169/153393500-b2313500-9dc0-4883-bcb9-9d9ef65f734c.png)
![image](https://user-images.githubusercontent.com/6979169/153393680-a6f42c9d-ade7-4390-8c52-5b34837a0ebb.png)

  
