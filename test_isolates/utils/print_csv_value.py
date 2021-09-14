#!/usr/local/bin/python

""" print_csv_value.py: 
		Prints an item in a CSV 
"""

import pandas as pd
import argparse

def print_csv_value(filename, field, row):
	df = pd.read_csv(filename, dtype=str, keep_default_na=False)

	if field not in df:
		raise Exception('CSV does not contain field: %s (%s)'%(field, filename))

	if len(df) < row+1:
		raise Exception('Row number larger than number of rows')

	print(df[field][row])

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="""
		Assert the first entry of a specified column in a csv file against a provided value.
		exitcode 1 on failure and 0 otherwise
	""")
	parser.add_argument('filename', help='path to csv')
	parser.add_argument('row', type=int, help='which row to access')
	parser.add_argument('field', help='name of column to check')

	args = parser.parse_args()

	print_csv_value(**vars(args))
