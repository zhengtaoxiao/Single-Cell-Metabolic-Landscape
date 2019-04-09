#!/usr/bin/env bash

###### download melanoma dataset #####

url='ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/suppl/GSE72056_melanoma_single_cell_revised_v2.txt.gz'
wget $url || curl $url -o GSE72056_melanoma_single_cell_revised_v2.txt.gz
gunzip GSE72056_melanoma_single_cell_revised_v2.txt.gz

#rename corrupt gene symbols due to incorrect conversion by Excel
# 1-Mar -> MARCH1
# 1-Dec -> DEC1
# 1-Sep -> SEPT1
# SEPT15 should be SEP15


command="sed"
for (( i = 15; i > 0; i-- )); do
	command="$command -e \"s/$i-Mar/MARCH$i/g\" -e \"s/$i-Sep/SEPT$i/g\" -e \"s/$i-Dec/DEC$i/g\" -e \"s/SEPT15/SEP15/g\""
done

eval "cat GSE72056_melanoma_single_cell_revised_v2.txt | $command > GSE72056_melanoma_single_cell_corrected.txt"
rm GSE72056_melanoma_single_cell_revised_v2.txt
rm GSE72056_melanoma_single_cell_revised_v2.txt.gz

###### download head and neck dataset #####
url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz"
wget $url || curl $url -o GSE103322_HNSCC_all_data.txt.gz
gunzip GSE103322_HNSCC_all_data.txt.gz


##move the data to directory
if [ ! -d "dataset" ]
then
  mkdir dataset
fi

mv {GSE103322_HNSCC_all_data.txt,GSE72056_melanoma_single_cell_corrected.txt} dataset/
