#!/usr/local/bin python
import sys
#import argparse
import pandas as pd
import allel #pip install scikit-allel to install this module for analysis of large scale genetic variation data
import csv

# parser = argparse.ArgumentParser(description='VCFCompare V2.0')
# parser.add_argument('-G', type=str,   required=True,  metavar='<str>',                  help="* ref.fa")
# parser.add_argument('-Q', type=int,   required=True,  metavar='<int>',                  help="* read length")
# parser.add_argument('-O', type=str,   required=True,  metavar='<str>',                  help="* output prefix")

golden_vcf = allel.vcf_to_dataframe(sys.argv[1], ['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT', 'QUAL'], alt_number=1) # Storing the information in the VCF file into a dataframe
query_vcf = allel.vcf_to_dataframe(sys.argv[2], ['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT', 'QUAL'], alt_number=1)
golden_variants_list = golden_vcf.values.tolist() #Each variant info will be saved onto a list, creating a list of list
query_variants_list = query_vcf.values.tolist()

header = ['TRUTH.TOTAL', 'TP', 'FP', 'FN', 'QUERY.TOTAL', 'Recall', 'Precision']
list = []

#Totals
Truth_Total = len(golden_variants_list)
print ("Total Truth VCF records: " + str(Truth_Total))
Truth_Pass_Total = len([x for x in golden_variants_list if x[4] > 2]) #Filter variants by quality score. Ask how what is the minimum QUAL a variant must have to get a 'PASS'
print ("Non-reference VCF records: " + str(Truth_Pass_Total))
Query_Total = len(query_variants_list)
print ("Total Query VCF records: " + str(Query_Total))
Query_Pass_Total = len([x for x in query_variants_list if x[4] > 0.003])
print ("Non-reference VCF records: " + str(Query_Pass_Total))


#Calls
TP = [x for x in query_variants_list if x in golden_variants_list]
len_TP = len(TP)
FP = [x for x in query_variants_list if x not in golden_variants_list]
len_FP = len(FP)
FN = [x for x in golden_variants_list if x not in query_variants_list]
len_FN = len(FN)

#Calculations
Recall = len_TP/(len_TP+len_FN)
Precision = len_TP/(len_TP+len_FP)

list.extend((Truth_Total, len_TP,len_FP,len_FN, Query_Total, Recall, Precision))

with open("test.csv", "w", newline='') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(header)
    writer.writerow(list)
