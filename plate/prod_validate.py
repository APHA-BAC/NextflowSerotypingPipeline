import process_plate
import pandas as pd
import os
import argparse
import shutil
from shutil import copy2
import subprocess
import logging

DEFAULT_IMAGE = "jguzinski/salmonella-seq:master"
DEFAULT_READS_DIR = os.path.expanduser("~/wgs-reads/")
READS_250_LOCATION = os.path.expanduser("~/mnt/Salmonella/BAC3_NGS_Archive/Validation/250_validation_Nov23/")
EXPECTED_RESULTS_PATH = "~/NextflowSerotypingPipeline/validation250/250_validation_EXPECTED_summary.csv"


def copy2_verbose(src, dst):
    print('Copying {0}'.format(src))
    copy2(src,dst)

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

def load_summary_table(csv_path):
    """ Returns a DataFrame containing Consensus and unique Isolate_ID columns """
    # Load
    df = pd.read_csv(csv_path)

    # Validate
    columns = df.columns.to_list()
    if not 'Isolate_ID' in columns:
        raise Exception(f"Isolate_ID column missing: {csv_path}")

    if not 'Consensus' in columns:
        raise Exception(f"Consensus column missing: {csv_path}")

    if not df.Isolate_ID.is_unique:
        raise Exception(f"Isolate_ID column is not unique: {csv_path}")

    return df

def analyse_results(expected_csv_path, actual_csv_path):
    """ Returns a merged DataFrame containing an Outcome column that indicates consistency """
    # Load
    expected_df = load_summary_table(expected_csv_path)[["Isolate_ID", "Consensus", "LIMS_Status", "LIMS_Reason"]]
    actual_df = load_summary_table(actual_csv_path)

    # Rename Columns
    expected_df = expected_df.rename(columns={"Consensus": "ExpectedConsensus", "LIMS_Status": "ExpectedStatus", "LIMS_Reason": "ExpectedReason"})
    actual_df = actual_df.rename(columns={"Consensus": "ActualConsensus", "LIMS_Status": "ActualStatus", "LIMS_Reason": "ActualReason"})

    # Join
    merged = expected_df.merge(actual_df, on='Isolate_ID', how='left')

    # Evaluate
    merged['Outcome'] = 'unset'
    merged.loc[merged['ActualConsensus']==merged['ExpectedConsensus'], 'Outcome'] = 'success'
    merged.loc[merged['ActualConsensus']!=merged['ExpectedConsensus'], 'Outcome'] = 'fail'
    merged.loc[merged['ActualConsensus'].isna(), 'success'] = 'missing'
    
    merged.loc[merged['ActualStatus']==merged['ExpectedStatus'], 'Outcome'] = 'success'
    merged.loc[merged['ActualStatus']!=merged['ExpectedStatus'], 'Outcome'] = 'fail'
    merged.loc[merged['ActualStatus'].isna(), 'success'] = 'missing'
    

    # Outcome
    return merged

def validate(image, expected_path, reads_path, fastq_path):
    
    # Copy the 250 isolates and create all the necessary directories
    full_reads_path = os.path.expanduser("~/wgs-reads/validation_test")
    
    shutil.copytree(fastq_path, full_reads_path, copy_function=copy2_verbose)
    
    full_results_path = "~/wgs-results/" + "validation_test/"
    full_results_path = os.path.expanduser(full_results_path)

    if os.path.exists(full_results_path):
        raise Exception("Results path already exists; Path must not exist: ", full_results_path)
    
    # Start running the 250 through the pipeline
    run(['python', 'process_plate.py', '-r', 'validation_test', '--image', image])
    actual_csv_path = full_results_path + "validation_test_SummaryTable_plusLIMS.csv"

    # Compare results of the run against the expected results
    merged = analyse_results(expected_path, actual_csv_path)
    outcome_csv_path = os.path.expanduser("~/NextflowSerotypingPipeline/validation250/validation_outcome.csv")
    merged.to_csv(outcome_csv_path)

    if (merged["Outcome"] == "success").all():
        print("Successful validation! :)")
    else:
        print("Validation failed :(")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Validate a docker image")
    parser.add_argument("--image", default=DEFAULT_IMAGE, help="Docker image")

    args = parser.parse_args()

    validate(args.image, EXPECTED_RESULTS_PATH, DEFAULT_READS_DIR, READS_250_LOCATION)
