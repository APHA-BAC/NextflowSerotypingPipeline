import process_plate
import pandas as pd
import os

pd.set_option("display.max_rows", None, "display.max_columns", None)

DEFAULT_IMAGE = "jguzinski/salmonella-seq:master"

def load_summary_table(csv_path):
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
    # Load
    expected_df = load_summary_table(expected_csv_path)[["StrainID", "Consensus"]]
    actual_df = load_summary_table(actual_csv_path)

    # Rename Columns
    expected_df = expected_df.rename(columns={"Consensus": "ExpectedConsensus"})
    actual_df = actual_df.rename(columns={"Consensus": "ActualConsensus"})

    # Join
    merged = expected_df.merge(actual_df, on='StrainID')

    # Evaluate
    merged['Outcome'] = 'unset'
    merged.loc[merged['ActualConsensus']==merged['ExpectedConsensus'], 'Outcome'] = 'success'
    merged.loc[merged['ActualConsensus']!=merged['ExpectedConsensus'], 'Outcome'] = 'fail'
    merged.loc[merged['ActualConsensus'].isna(), 'success'] = 'missing'

    # Outcome
    return merged

def validate(
    expected_csv_path
    image=DEFAULT_IMAGE
):
    reads_path = "/home/joshuapotter/wgs-reads/validation_test"
    results_path = "/home/joshuapotter/wgs-results/validation_test"

    #process_plate.run_pipeline(reads_path, results_path, "validation_test", image)

    results_csv_path = "../validation250/validation250_fastpTrimmed_SummaryTable.csv"
    expected_csv_path = "../validation250/validation_test_SummaryTable.csv"

    merged = analyse_results(expected_csv_path, results_csv_path)

    if (merged["Outcome"] == "success").all():
        print("Successful validation! :)")
    else:
        print("Validation failed :(")

validate()
quit()

expected_dict = {}

for i, row in expected_df.iterrows():
    strain_id = row["StrainID"]
    consensus = row["Consensus"]
    if strain_id in expected_dict:
        raise Exception(f"Duplicate strain ID in expected_dict {strain_id}")
    expected_dict[strain_id] = consensus



# print(results_df)

results_dict = {}

for i, row in results_df.iterrows():
    strain_id = row["StrainID"]
    consensus = row["Consensus"]
    if strain_id in results_dict:
        raise Exception(f"Duplicate strain ID in results_dict {strain_id}")
    results_dict[strain_id] = consensus

# print(results_dict)

outcomes = []

for strain_id, consensus in expected_dict.items():
    if strain_id not in results_dict:
        outcome = "MISSING"
        actual = "missing"
    elif results_dict[strain_id] == consensus:
        outcome = "SUCCESS"
        actual = results_dict[strain_id]
    else:
        outcome = "FAIL"
        actual = results_dict[strain_id]
    outcomes.append({"StrainID":strain_id,
        "Outcome":outcome,
        "Expected":consensus,
        "Actual":actual})

outcome_df = pd.DataFrame(outcomes).sort_values("StrainID")

if (outcome_df["Outcome"] == "SUCCESS").all():
    print("Successful validation! :)")
else:
    print("Validation failed :(")