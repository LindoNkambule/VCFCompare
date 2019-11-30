#!/usr/bin/env python
import os
import argparse
import csv
from fun import infoExtract, createLists, concordance


class variantSeparation:
    def __init__(self, truth, query):
        self.truth = truth
        self.query = query


    def SNVs(self):
        snv = []
        # Truth and Query
        tsnv = createLists(self.truth).snvList()
        truthSNVs = len(tsnv)
        qsnv = createLists(self.query).snvList()
        querySNVs = len(qsnv)

        # Calls, Recall, and Precision
        stp, sfp, sfn = concordance(tsnv, qsnv).variantCalls()
        lenstp = len(stp)
        lensfp = len(sfp)
        lensfn = len(sfn)

        try:
            snvRecall = lenstp/(lenstp+lensfn)
            snvPrecision = lenstp/(lenstp+lensfp)
            snv.extend(('SNV', truthSNVs, lenstp, lensfp, lensfn, querySNVs, snvRecall, snvPrecision))
        except Exception:
            pass
        return snv


    def INDELs(self):
        indel = []
        # Truth and Query
        tindel = createLists(self.truth).indelList()
        truthINDELs = len(tindel)
        qindel = createLists(self.query).indelList()
        queryINDELs = len(qindel)

        # Calls, Recall, and Precision
        itp, ifp, ifn = concordance(tindel, qindel).variantCalls()
        lenitp = len(itp)
        lenifp = len(ifp)
        lenifn = len(ifn)

        try:
            indelRecall = lenitp/(lenitp+lenifn)
            indelPrecision = lenitp/(lenitp+lenifp)
            indel.extend(('INDEL', truthINDELs, lenitp, lenifn, lenifn, queryINDELs, indelRecall, indelPrecision))
        except Exception:
            pass
        return indel


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
