#!/usr/local/bin/python

""" assert_first_csv_row.py: 
		Assert the first entry of a specified column in a csv file against a provided value. 
"""

import pandas as pd
import argparse

def assert_first_csv_row(filename, field, value):
	df = pd.read_csv(filename, dtype=str, keep_default_na=False)

	if field not in df:
		raise Exception('CSV does not contain field: %s (%s)'%(field, filename))

	if not len(df):
		raise Exception('No rows found in %s'%filename)

	if not df[field][0] == value:
		raise Exception("First entry in %s field mismatch \nGot: %s \nExpected: %s" % (field, df[field][0], value))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="""
		Assert the first entry of a specified column in a csv file against a provided value.
		exitcode 1 on failure and 0 otherwise
	""")
	parser.add_argument('filename', help='path to csv')
	parser.add_argument('field', help='name of column to check')
	parser.add_argument('value', help='what the column should read')

	args = parser.parse_args()

	assert_first_csv_row(**vars(args))
