#!/usr/bin/env python


class infoExtract:
    def __init__(self, vcf):
        self.vcf = vcf


    def alleles(self):
        import allel
        vcfInfo = allel.vcf_to_dataframe(self.vcf, ['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT'], alt_number=1)
        vcfList = vcfInfo.values.tolist()
        return vcfList


class createLists():
    def __init__(self, variantsList):
        self.variantsList = variantsList


    def snvList(self):
        snvList = []
        for variant in self.variantsList:
            ref = len(str(variant[2]))
            alt = len(str(variant[3]))
            if(ref == 1 and alt == 1):
                snvList.append(variant)
        return snvList


    def indelList(self):
        indelList = []
        for variant in self.variantsList:
            ref = len(str(variant[2]))
            alt = len(str(variant[3]))
            if (ref > 1 or alt > 1):
                indelList.append(variant)
        return indelList


    def snvINDELlists(self):
        snvList = []
        indelList = []
        for variant in self.variantsList:
            ref = len(str(variant[2]))
            alt = len(str(variant[3]))
            if (ref > 1 or alt > 1):
                indelList.append(variant)
            else:
                snvList.append(variant)
        return snvList, indelList

class concordance():
    def __init__(self, truth, query):
        self.truth = truth
        self.query = query


    def variantCalls(self):
        Truth_set = set(map(tuple, self.truth))
        Query_Set = set(map(tuple, self.query))
        TPs = Truth_set.intersection(Query_Set)
        FPs = Query_Set.difference(Truth_set)
        FNs = Truth_set.difference(Query_Set)
        return TPs, FPs, FNs
