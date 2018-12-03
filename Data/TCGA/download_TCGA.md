## step1 get the "Manifest File" and "Sample Sheet File" from:
    https://portal.gdc.cancer.gov/
    
    # These two files are provided in HNSC folder.
## step2 download the GDC tools from
    https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
## step3 filtering the FPKM-UQ file using the following command:
    head -n1 gdc_manifest_20180330_170153.txt > gdc_manifest_20180330_170153_FPKM-UQ.txt

    grep "FPKM-UQ" gdc_manifest_20180330_170153.txt >> gdc_manifest_20180330_170153_FPKM-UQ.txt
## step4 using the gdc tools download the expression for each patient:
    gdc-client download -m gdc_manifest_20180330_170153_FPKM-UQ.txt
## step5 extract the expression data for each patient and save as the matrix file
    python generate_matrix.py gdc_sample_sheet.2018-03-30.tsv exp.txt type.txt
    
