#!/usr/bin/env python
import os
import argparse
import csv
import fun

def main():
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='VCFCompare V1.0')
    parser.add_argument("-t", "--truth", help="Truth VCF")
    parser.add_argument("-q", "--query", help="Query VCF")
    parser.add_argument("-o", "--output", help="Output CSV file prefix")
    args = parser.parse_args()

    #Error handling
    if not args.truth:
        raise Exception("Please specify the truth VCF file using -t or --truth")
    if not args.query:
        raise Exception("Please specify the query VCF file using -q or --query ")
    if not args.output:
        raise Exception("Please specify an output prefix using -o or --output")

    if not os.path.exists(args.truth):
        raise Exception("Input file {} does not exist.".format(args.truth))
    if not os.path.exists(args.query):
        raise Exception("Input file {} does not exist.".format(args.query))

    truth = fun.vcfExtract(args.truth)
    query = fun.vcfExtract(args.query)
    truth_list = fun.vcfDFtoList(truth)
    query_list = fun.vcfDFtoList(query)

    #Output file
    header = ['Type', 'TRUTH.TOTAL', 'TP', 'FP', 'FN', 'QUERY.TOTAL', 'Recall', 'Precision']
    snv = []
    indel = []

    #Truth SNVs and INDELs
    tsnv = []
    tindel = []
    fun.snvINDELlists(tsnv, tindel, truth_list)
    truthSNVs = len(tsnv)
    truthINDELs = len(tindel)

    #Query SNVs and INDELs
    qsnv = []
    qindel = []
    fun.snvINDELlists(qsnv, qindel, query_list)
    querySNVs = len(qsnv)
    queryINDELs = len(qindel)

    #Totals
    Truth_Total = len(truth_list)
    print ("Total Truth VCF records: {}".format(Truth_Total))
    Query_Total = len(query_list)
    print ("Total Query VCF records: {}".format(Query_Total))

    #Calls, Recall, and Precision
    #1.SNVs
    stp, sfp, sfn = fun.variantCalls(tsnv, qsnv)
    lenstp = len(stp)
    lensfp = len(sfp)
    lensfn = len(sfn)

    snvRecall = lenstp/(lenstp+lensfn)
    snvPrecision = lenstp/(lenstp+lensfp)

    #2.INDELs
    itp, ifp, ifn = fun.variantCalls(tindel, qindel)
    lenitp = len(itp)
    lenifp = len(ifp)
    lenifn = len(ifn)

    indelRecall = lenitp/(lenitp+lenifn)
    indelPrecision = lenitp/(lenitp+lenifp)

    snv.extend(('SNV', truthSNVs, lenstp, lensfp, lensfn, querySNVs, snvRecall, snvPrecision))
    indel.extend(('INDEL', truthINDELs, lenitp, lenifn, lenifn, queryINDELs, indelRecall, indelPrecision))

    csvfile = args.output + ".csv"
    with open(csvfile, "w", newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(header)
        writer.writerow(snv)
        writer.writerow(indel)
    f.close()

if __name__ == '__main__':
    main()
