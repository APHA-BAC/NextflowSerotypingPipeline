[![Build Status](https://apha.teamcity.com/app/rest/builds/buildType:(id:NextflowSerotypingPipeline_Pipeline)/statusIcon)](https://apha.teamcity.com/viewType.html?buildTypeId=NextflowSerotypingPipeline_Pipeline&guest=1)

# Introduction

Salmonella whole genome sequencing (WGS) serotyping pipeline was developed by APHA and is written in [Nextflow](https://www.nextflow.io/). The pipeline compares outputs from several publicly available serotyping tools to increase performance and includes bespoke differentiation of strains important for veterinary surveillance.

# Installation
To install the Nextflow salmonella serotyping pipeline:

First clone the rpository into a local directory:

```
git clone https://github.com/APHA-BAC/NextflowSerotypingPipeline.git

```
And then run the following commands while in the repository:

```
$ cd install
$ bash install.bash

```

This script installs the following pipeline component software:
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
- `SISTR`
- `SRST2`
- `sratoolkit`

# Running the pipeline

To run the pipeline on a batch of samples, the raw `.fastq.gz` files must be stored in `~/wgs-reads/<runID>`.  Each read-pair sample is represented by a pair of files named `*_R1.fastq.gz` and `*_R2.fastq.gz`. For example, to batch two samples named `salmonella_a` and `salmonella_b`, a `~/wgs-reads/<runID>` directory containing four files is required: `salmonella_a_R1.fastq.gz`, `salmonella_a_R2.fastq.gz`,  `salmonella_b_R1.fastq.gz` and `salmonella_b_R2.fastq.gz`.

Then, to run the pipeline from the terminal call:
```
$ nextflow run SCE3_pipeline_update.nf --runID <runID>

```

Pipeline output is stored in  `~/wgs-results/<runID>/` and contains:
- TODO

## Run from docker

**Note:** While running from the terminal is the easiest method for developers and data analysts, the pipeline can also be run from docker. This method has the benefit of working across platforms while guaranteeing consistency with automated tests (see below).

A docker image containing all required dependencies is provided [here](https://hub.docker.com/r/jguzinski/salmonella-seq).

When running the pipeline on data stored locally, you can run the following command from the plate directory of the repository which will download the docker image if it is not downloaded
and will then run the pipeline in a docker container on your locally stored fastq files:
```
$ cd plate
$ python process_plate.py --local 1 --runID <runID>

```
# Run options

These are the run options for process_plate.py:

- --s3_uri allows you to run the plate on data stored on an AWS s3 bucket
- --reads-dir Sets the directory for where the directory containing the reads are stored. Default is: "~/wgs-reads"
- --results-dir Sets the directory for where the directory containing the results are stored. Default is: "~/wgs-results"
- --image Select a specific docker image to use. Default is: "jguzinski/salmonella-seq:prod"
- --local Set to 1 if your reads are in a local directory. Default is: 0
- --runID The name of the run which will also be the name of the directory for the results. Only needed if --local is set to 1

# Pipeline Algorithm

The pipelines processes data in three stages, as shown below. During the preprocessing stage; low quality bases and adapter sequences are removed from the fastq sample files and then subsampled to a maximum of 3M reads. Following this, the analysis stage runs multiple serotyping tools in parallel.
This strategy has been demonstrated to provide more accurate serovar detection than any tools running individually.
Bespoke typing of vaccine strains and particular servoars of interest to APHA is included as part of this analysis.
Outputs from each tool are compared in the consensus call step.
The final postprocessing stage assigns an `Outcome` to each sample by analysing data gathered during the analysis stage. The following `Outcome`s are used to signify subsequent lab processing steps:

- **Pass**: The sample contains Salmonella and has a serovar assigned
- **Check Required**: Further scrutiny of the data is required
- **Inconclusive**: The sample contains insufficient data volumes for analysis, is contaminated, or has low assembly quality.

![image](https://user-images.githubusercontent.com/6979169/154251677-43d55d28-24bb-4dc2-9def-61322ba58629.png)

# Validation

This pipeline has been internally validated and approved by the APHA validation team against a large dataset in excess of 3,500 samples as well as non-Salmonella isolates to determine sensitivity and specificity.


# Automated Tests

The automated tests provided here ensure the software runs as expected. If you make changes to the algorithm, it is **strongly** reccomended that you run these tests to verify the pipeline is behaving as intended. The tests are also automatically run by [TEAMCITY](https://apha.teamcity.com/viewType.html?buildTypeId=NextflowSerotypingPipeline_Pipeline&guest=1) on each pull-request.

The automated tests assert the correct serovar is assigned to known samples. These are called `inclusivity` tests. To run an inclusivity test, call:
```
$ bash -e test_isolates/jobs/inclusivity.bash 0

```

# Release Process

To release a new version of the software, the `master` branch needs only to be merged into the `prod` branch. To perform this merge, a pull-request from the `master` branch into the `prod` branch needs to be made. Approval of pull-requests to `prod` is made by the CODEOWNER (Liljana Petrovska). The CODEOWNER is responsible for ensuring the code conforms to the reliability tests. A positive test result is required for approval.

To release a new version of the software:
1. A developer makes a pull-request from the `master` to the `prod` branch. The CODEOWNER is automatically notified by e-mail.
1. The CODEOWNER ensures the reliability tests pass on the `master` branch and reviews the code changes.
1. The CODEOWNER approves the pull-request if they are satisfied, or requests changes.
1. The dev merges the `master` branch into `prod`
1. Following approval, the developer tags the current head of `master` as the next version. Versions are numbered incrementally with integers, for example `v1`, `v2`, etc. This can be performed by navigating to the github `master` branch and selecting `Create a release`

![image](https://user-images.githubusercontent.com/6979169/153393500-b2313500-9dc0-4883-bcb9-9d9ef65f734c.png)
![image](https://user-images.githubusercontent.com/6979169/153393680-a6f42c9d-ade7-4390-8c52-5b34837a0ebb.png)
