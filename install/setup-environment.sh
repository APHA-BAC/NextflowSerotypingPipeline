set -e

#################################################################################################
# Remaining setup
mkdir $HOME/wgs-reads
mkdir $HOME/wgs-results

#################################################################################################
# Summary table script

mkdir $HOME/summary # Changed to $HOME for SC3 - We don't want summary installed in /home/summary

cp ./summaryTable_reworked.py $HOME/summary
cp ./ebgs.csv $HOME/summary
cp ./lims_rules/lims_rules.py $HOME/summary
cp ./lims_rules/serogroup_lookup_table.tsv $HOME/summary
cp ./lims_rules/subgenus_lookup_table.tsv $HOME/summary
