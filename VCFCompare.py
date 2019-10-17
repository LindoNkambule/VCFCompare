#!/usr/local/bin python
import os
import argparse
import allel
import csv

def main():
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='VCFCompare V1.0')
    parser.add_argument("-t", "--truth", help="Truth VCF")
    parser.add_argument("-q", "--query", help="Query VCF")
    parser.add_argument("-o", "--output", help="Output CSV file")
    args = parser.parse_args()

    #Error handling
    if not args.truth:
        raise Exception("Please specify an output prefix using -t or --truth")
    if not args.query:
        raise Exception("Please specify an output prefix using -q or --query ")
    if not args.output:
        raise Exception("Please specify an output prefix using -o or --output")

    if not os.path.exists(args.truth):
        raise Exception("Input file %s does not exist." % args.truth)
    if not os.path.exists(args.query):
        raise Exception("Input file %s does not exist." % args.query)

    truth = allel.vcf_to_dataframe(args.truth, ['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT'], alt_number=1)
    query = allel.vcf_to_dataframe(args.query, ['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT'], alt_number=1)
    truth_list = truth.values.tolist()
    query_list = query.values.tolist()

    #Output file
    header = ['Type', 'TRUTH.TOTAL', 'TP', 'FP', 'FN', 'QUERY.TOTAL', 'Recall', 'Precision']
    snv = []
    indel = []

    #Truth SNVs and INDELs
    tsnv = []
    tindel = []
    for variant in truth_list:
        ref = variant[2]
        alt = variant[3]
        if (len(str(ref)) > 1 or len(str(alt)) > 1):
            tindel.append(variant)
        else:
            tsnv.append(variant)

    truthSNVs = len(tsnv)
    truthINDELs = len(tindel)

    #Query SNVs and INDELs
    qsnv = []
    qindel = []
    for variant in query_list:
        ref = variant[2]
        alt = variant[3]
        if (len(str(ref)) > 1 or len(str(alt)) > 1):
            qindel.append(variant)
        else:
            qsnv.append(variant)

    querySNVs = len(qsnv)
    queryINDELs = len(qindel)

    #Totals
    Truth_Total = len(truth_list)
    print ("Total Truth VCF records: " + str(Truth_Total))
    Query_Total = len(query_list)
    print ("Total Query VCF records: " + str(Query_Total))

    #Calls, Recall, and Precision
    #1.SNVs
    stp = [x for x in qsnv if x in tsnv]
    lenstp = len(stp)
    sfp = [x for x in qsnv if x not in tsnv]
    lensfp = len(sfp)
    sfn = [x for x in tsnv if x not in qsnv]
    lensfn = len(sfn)

    snvRecall = lenstp/(lenstp+lensfn)
    snvPrecision = lenstp/(lenstp+lensfp)

    #2.INDELs
    itp = [x for x in qindel if x in tindel]
    lenitp = len(itp)
    ifp = [x for x in qindel if x not in tindel]
    lenifp = len(ifp)
    ifn = [x for x in tindel if x not in qindel]
    lenifn = len(ifn)

    indelRecall = lenitp/(lenitp+lenifn)
    indelPrecision = lenitp/(lenitp+lenifp)

    snv.extend(('SNV', truthSNVs, lenstp, lensfp, lensfn, querySNVs, snvRecall, snvPrecision))
    indel.extend(('INDEL', truthINDELs, lenitp, lenifn, lenifn, queryINDELs, indelRecall, indelPrecision))

    #list.extend((Truth_Total, len_TP,len_FP,len_FN, Query_Total, Recall, Precision))

    csvfile = args.output + ".csv"
    with open(csvfile, "w", newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(header)
        writer.writerow(snv)
        writer.writerow(indel)
    f.close()

if __name__ == '__main__':
   main()
