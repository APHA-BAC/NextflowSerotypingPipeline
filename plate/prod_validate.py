import process_plate
import pandas as pd
import os
import argparse

DEFAULT_IMAGE = "jguzinski/salmonella-seq:master"
DEFAULT_READS_DIR = os.path.expanduser('~/mnt/Salmonella/BAC3_NGS_Archive/Salmonella/Validation_panel_1/')
DEFAULT_RESULTS_DIR = os.path.expanduser('~/wgs-results/validation_test/')
DEFAULT_EXPECTED_CSV_PATH = '../validation250/validation250_fastpTrimmed_SummaryTable.csv'
DEFAULT_OUTCOME_PATH = './outcome.csv'

def load_summary_table(csv_path):
    """ Returns a DataFrame containing Consensus and unique StrainID columns """
    # Load
    df = pd.read_csv(csv_path)

    # Validate
    columns = df.columns.to_list()
    if not 'StrainID' in columns:
        raise Exception(f"StrainID column missing: {csv_path}")

    if not 'Consensus' in columns:
        raise Exception(f"Consensus column missing: {csv_path}")

    if not df.StrainID.is_unique:
        raise Exception(f"StrainID column is not unique: {csv_path}")

    return df

def analyse_results(expected_csv_path, actual_csv_path):
    """ Returns a merged DataFrame containing an Outcome column that indicates consistency """
    # Load
    expected_df = load_summary_table(expected_csv_path)[["StrainID", "Consensus"]]
    actual_df = load_summary_table(actual_csv_path)

    # Rename Columns
    expected_df = expected_df.rename(columns={"Consensus": "ExpectedConsensus"})
    actual_df = actual_df.rename(columns={"Consensus": "ActualConsensus"})

    # Join
    merged = expected_df.merge(actual_df, on='StrainID', how='left')

    # Evaluate
    merged['Outcome'] = 'unset'
    merged.loc[merged['ActualConsensus']==merged['ExpectedConsensus'], 'Outcome'] = 'success'
    merged.loc[merged['ActualConsensus']!=merged['ExpectedConsensus'], 'Outcome'] = 'fail'
    merged.loc[merged['ActualConsensus'].isna(), 'success'] = 'missing'

    # Outcome
    return merged

def validate(
    reads_path,
    results_path,
    expected_csv_path,
    outcome_csv_path,
    image
):
    # process_plate.run_pipeline(reads_path, results_path, "validation_test", image)
    actual_csv_path = results_path + "/validation_test_SummaryTable.csv"

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
