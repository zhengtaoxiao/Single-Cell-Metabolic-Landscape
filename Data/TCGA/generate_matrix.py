#!/bin/env python
import sys
import pandas as pd
import os
import gzip

usage = """
#read the expression files and sample table , then generate expression matrix file
#usage:python3 generate_exp.py sample_sheet_file tumor_exp.txt tumor_type.txt
#inputFile: sample_sheet_file
#outputFile: tumor_exp.txt is for saving the expression matrix.
#            tumor_type.txt is annotation of the tumor type for each patient
"""
if(len(sys.argv)!=4): print(usage)

sample_sheet_file = sys.argv[1]
out_file = sys.argv[2]
tumor_type_file = sys.argv[3]
sample_sheet = pd.read_table(sample_sheet_file,header=0)
#filter the FPKM-UQ files
selection = sample_sheet["File Name"].str.contains('FPKM-UQ')
sample_sheet = sample_sheet[selection]
#select tumor and normal
selection = sample_sheet["Sample Type"].str.contains("Tumor|Normal|Metastatic")
sample_sheet = sample_sheet[selection]
#get the data and output as data.matrix
expression_dict = {}
for row in sample_sheet.iterrows():
	file_path = os.path.join(row[1]["File ID"],row[1]["File Name"])
	name = row[1]["Sample ID"]
	expression_dict[name] = []
	with gzip.open(file_path) as fin:
		for i in fin:
			expression_dict[name].append(i.decode('utf-8').strip().split())
#check if all file same gene name list
gene_name_list = [i[0] for i in expression_dict[name]]
for v in expression_dict.values():
	tmp = [i[0] for i in v]
	if tmp != gene_name_list:
		sys.stderr.write("gene name is not same!")
		sys.exit()
#write the data.matrix
data_frame = pd.DataFrame(index=gene_name_list,columns=list(expression_dict.keys()),dtype=str)
#fill the data
for i in expression_dict.keys():
	data_frame[i] = [j[1] for j in expression_dict[i]]
#convert to float, and remove the genes with zeros in all samples
data_frame = data_frame.astype(float)
mask = (data_frame<=0).all(axis=1)
data_frame[~mask].to_csv(out_file,sep="\t",index_label="GeneID")
# #out the tumortype
sample_sheet.to_csv(tumor_type_file,sep="\t",index=False,columns=["Sample ID","Sample Type"])

