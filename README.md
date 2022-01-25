[![Build Status](https://apha.teamcity.com/app/rest/builds/buildType:(id:NextflowSerotypingPipeline_Pipeline)/statusIcon)](https://apha.teamcity.com/viewType.html?buildTypeId=NextflowSerotypingPipeline_Pipeline&guest=1)

Introduction

The Salmonella serotyping pipeline described has been developed to replace the classical serotyping process at the APHA by a method based on extracting the serotyping information from the whole genome sequencing (WGS) data of the Salmonella isolates. The pipeline benefits from several publicly available serotyping tools developed by other institutions which are combined in order to increase the reliability of the results. This document describes a newer version of the serotyping pipeline that was re-coded from the previous version of the pipeline to be compatible with Nextflow programming language and the Amazon Web Services (AWS) Scientific Computing Environment version 3 (SCEv3). Nextflow is a programming language that allows running different software packages in parallel, using all available Central Processing Units (CPUs). This highly efficient large scale parallelization results in improved (shorter) time to result per sample and for all samples combined. Software written in different programming languages (such as Bash, Python, or Perl), as is the case with the Salmonella serotyping pipeline, is adapted to become fully operational within the Nextflow framework; furthermore, Nextflow is highly compatible to efficiently operate within the SCEv3 architecture that the APHA is currently in the process of implementing.

Installation
To install the Nextflow Salmonella serotyping pipeline on blank state VM:
1)	Download the install-salmonella-pipeline-SCEv3.sh script (attached) into /home/$USER
2)	Open terminal and type cd /home/$USER
3)	In the terminal type sh install-salmonella-pipeline-SCEv3.sh
4)	When prompted press the Enter key or type “y” or “yes”
5)	Installation process takes about 20 minutes
6)	The VM will reboot after the installation has been completed
7)	No other action will be required by the user after completion of the installation wrapper script and the pipeline will be ready for utilization   

Utilization
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
