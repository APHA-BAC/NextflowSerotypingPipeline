import process_plate
import pandas as pd
import os

pd.set_option("display.max_rows", None, "display.max_columns", None)

def load_summary_table(csv_path):
    # Load
    df = pd.read_csv(results_csv_path)

    # Validate
    columns = df.columns.to_list()
    if not 'StrainID' in columns:
        raise Exception(f"StrainID column missing: {csv_path}")

    if not 'Consensus' in columns:
        raise Exception(f"Consensus column missing: {csv_path}")

    if not df.StrainID.is_unique:
        raise Exception(f"StrainID column is not unique: {csv_path}")

    return df

# process_plate.run_pipeline("/home/joshuapotter/wgs-reads/validation_test", "/home/joshuapotter/wgs-results/validation_test", "validation_test", "jguzinski/salmonella-seq:master")

results_csv_path = "../validation250/validation250_fastpTrimmed_SummaryTable.csv"
expected_csv_path = "../validation250/validation_test_SummaryTable.csv"

# Load
expected_df = load_summary_table(results_csv_path)[["StrainID", "Consensus"]]
results_df = load_summary_table(expected_csv_path)

# Rename Columns
expected_df = expected_df.rename(columns={"Consensus": "ExpectedConsensus"})
results_df = results_df.rename(columns={"Consensus": "ActualConsensus"})
actual_df = results_df

# Join
merged = expected_df.merge(results_df, on='StrainID')

# Evaluate
merged['Outcome'] = 'unset'
merged.loc[merged['ActualConsensus']==merged['ExpectedConsensus'], 'Outcome'] = 'success'
merged.loc[merged['ActualConsensus']!=merged['ExpectedConsensus'], 'Outcome'] = 'fail'
merged.loc[merged['ActualConsensus'].isna(), 'success'] = 'missing'

# Outcome
if (merged["Outcome"] == "success").all():
    print("Validation successful")
else:
    print("Validation failed :(")


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