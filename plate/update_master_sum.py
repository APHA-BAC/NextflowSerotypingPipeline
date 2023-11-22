import pandas as pd


def update_summary(master_sum, new_sum):
	master_df = pd.read_csv(master_sum)
	new_df = pd.read_csv(new_sum)

	master_df = pd.concat([master_df, new_df])
	master_df.to_csv(master_sum, index=False)

