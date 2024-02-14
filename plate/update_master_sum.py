import pandas as pd
import re


def update_summary(master_sum,new_sum):
	master_df=pd.read_csv(master_sum)
	new_df=pd.read_csv(new_sum)
	sum_info=re.search("(\d{6})_APHA_(.*_\d{4})",str(new_sum))
	date=sum_info.group(1)
	plate_id=sum_info.group(2)

	day=date[0:2]
	month=date[2:4]
	year=date[4:6]
	date=day+"/"+month+"/"+year

	new_df.loc[0:,'Date']=date
	new_df.loc[0:,'Plate_ID']=plate_id

	master_df=pd.concat([master_df,new_df])
	master_df.to_csv(master_sum,index=False)