import process_plate
import pandas as pd
import os
import argparse

DEFAULT_IMAGE = "jguzinski/salmonella-seq:master"
DEFAULT_READS_DIR = os.path.expanduser('~/mnt/Salmonella/BAC3_NGS_Archive/Salmonella/validation_panel_Feb2022/')
DEFAULT_RESULTS_DIR = os.path.expanduser('~/wgs-results/validation_test/')
DEFAULT_EXPECTED_CSV_PATH = '../validation250/validation250_Feb2022_EXPECTED_SummaryTable_plusLIMS.csv'
DEFAULT_OUTCOME_PATH = './outcome.csv'

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
    expected_df = load_summary_table(expected_csv_path)[["Isolate_ID", "Consensus", "LIMS_Status", "LIMS_Reason", "MOST_Light"]]
    actual_df = load_summary_table(actual_csv_path)

    # Rename Columns
    expected_df = expected_df.rename(columns={"Consensus": "ExpectedConsensus", "LIMS_Status": "ExpectedStatus", "LIMS_Reason": "ExpectedReason", "MOST_Light": "ExpectedLight"})
    actual_df = actual_df.rename(columns={"Consensus": "ActualConsensus", "LIMS_Status": "ActualStatus", "LIMS_Reason": "ActualReason", "MOST_Light": "ActualLight"})

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

def validate(
    reads_path,
    results_path,
    expected_csv_path,
    outcome_csv_path,
    image
):
    # if os.path.exists(results_path):
    #     raise Exception("Results path already exists; Path must not exist: ", results_path)

    # os.makedirs(results_path)

    # process_plate.run_pipeline(reads_path, results_path, "validation_test", image)
    actual_csv_path = results_path + "/validation_test_SummaryTable_plusLIMS.csv"

    merged = analyse_results(expected_csv_path, actual_csv_path)

    merged.to_csv(outcome_csv_path)

    if (merged["Outcome"] == "success").all():
        print("Successful validation! :)")
    else:
        print("Validation failed :(")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Validate a docker image")
    parser.add_argument("--reads", default=DEFAULT_READS_DIR, help="Directory containing fastq reads")
    parser.add_argument("--results", default=DEFAULT_RESULTS_DIR, help="Output directory for writing sequenced results")
    parser.add_argument("--expected", default=DEFAULT_EXPECTED_CSV_PATH, help="Path to csv that contains expected results")
    parser.add_argument("--outcome", default=DEFAULT_OUTCOME_PATH, help="Output outcome csv")
    parser.add_argument("--image", default=DEFAULT_IMAGE, help="Docker image")

    args = parser.parse_args()

    validate(os.path.expanduser(args.reads), 
        os.path.expanduser(args.results), 
        os.path.expanduser(args.expected), 
        os.path.expanduser(args.outcome), 
        args.image
        )
