#!/usr/bin/env python
import os
import argparse
import csv
from fun import infoExtract, createLists, concordance, variantSeparation


def main():
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='VCFCompare V1.0')
    parser.add_argument("truth", help="Truth VCF")  # required positional arg
    parser.add_argument("query", help="Query VCF")  # required positonal arg
    parser.add_argument("-o", "--out", help="Output CSV file prefix")  # optional arg
    parser.add_argument("-t", "--type", help="The type of variant concordance: SNV or INDEL")  # optional arg
    args = parser.parse_args()

    # Error handling
    if not args.out:
        raise Exception("Please specify an output prefix using -o or --output")

    if not os.path.exists(args.truth):
        raise Exception("Input truth file {} does not exist.".format(args.truth))
    if not os.path.exists(args.query):
        raise Exception("Input query file {} does not exist.".format(args.query))

    truth_list = infoExtract(args.truth).alleles()
    query_list = infoExtract(args.query).alleles()


    # Totals
    Truth_Total = len(truth_list)
    print ("Total Truth VCF records: {}".format(Truth_Total))
    Query_Total = len(query_list)
    print ("Total Query VCF records: {}".format(Query_Total))

    # Output file
    header = ['Type', 'TRUTH.TOTAL', 'TP', 'FP', 'FN', 'QUERY.TOTAL', 'Recall', 'Precision']

    variants = variantSeparation(truth_list, query_list)

    if args.type == "SNV":
        snv = variants.SNVs()
        csvfile = args.out + ".SNV.csv"
    elif args.type == "INDEL":
        indel = variants.INDELs()
        csvfile = args.out + ".INDEL.csv"
    else:
        # NB!!! PARALLELIZE THE FOLOWING TWO TO SAVE SAVE
        snv = variants.SNVs()
        indel = variants.INDELs()
        csvfile = args.out + ".csv"

    with open(csvfile, "w", newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(header)
        if args.type == "SNV":
            writer.writerow(snv)
        elif args.type == "INDEL":
            writer.writerow(indel)
        else:
            writer.writerow(snv)
            writer.writerow(indel)
    f.close()


if __name__ == '__main__':
    main()
