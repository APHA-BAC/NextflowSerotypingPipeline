import process_plate
import pandas as pd
import os

# process_plate.run_pipeline("/home/joshuapotter/wgs-reads/validation_test", "/home/joshuapotter/wgs-results/validation_test", "validation_test", "jguzinski/salmonella-seq:master")

results_dir = "/home/joshuapotter/wgs-results/validation_test"

expected_df = pd.read_csv(os.path.expanduser("~/NextflowSerotypingPipeline/validation250/validation250_fastpTrimmed_SummaryTable.csv"))

# print(expected_df)

expected_dict = {}

for i, row in expected_df.iterrows():
    strain_id = row["StrainID"]
    consensus = row["Consensus"]
    if strain_id in expected_dict:
        raise Exception(f"Duplicate strain ID in expected_dict {strain_id}")
    expected_dict[strain_id] = consensus

results_df = pd.read_csv(os.path.join(results_dir ,"validation_test_SummaryTable.csv"))

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