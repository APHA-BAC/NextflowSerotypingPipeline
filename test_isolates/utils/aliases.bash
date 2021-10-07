#!/bin/bash
#
#================================================================
# aliases.bash
#================================================================
#
#% DESCRIPTION
#%    A number of aliases useful for testing
shopt -s expand_aliases

alias nextflowtest="bash test_isolates/utils/nextflowtest.bash"

alias print_path_to_summary_table="sh test_isolates/utils/print_path_to_summary_table.sh"

alias assert_first_csv_row="python test_isolates/utils/assert_first_csv_row.py"

alias combine_fastq="python test_isolates/utils/combine_fastq.py"

alias print_csv_value="python test_isolates/utils/print_csv_value.py"
