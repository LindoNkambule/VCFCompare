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
