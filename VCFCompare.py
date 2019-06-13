#!/usr/local/bin python
import sys
import pandas as pd
import allel #pip install scikit-allel to install this module for analysis of large scale genetic variation data

golden_vcf = allel.vcf_to_dataframe(sys.argv[1], ['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT'], alt_number=1) # Storing the information in the VCF file into a dataframe
query_vcf = allel.vcf_to_dataframe(sys.argv[2], ['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT'], alt_number=1)
golden_variants_list = golden_vcf.values.tolist() #Each variant info will be saved onto a list, creating a list of list
query_variants_list = query_vcf.values.tolist()

#Totals
Truth_Total = len(golden_variants_list)
print (Truth_Total)
Query_Total = len(query_variants_list)
print (Query_Total)

#Calls
TP = [x for x in query_variants_list if x in golden_variants_list]
print (len(TP))
FP = [x for x in query_variants_list if x not in golden_variants_list]
print (len(FP))
FN = [x for x in golden_variants_list if x not in query_variants_list]
print (len(FN))
#print(sys.argv)
