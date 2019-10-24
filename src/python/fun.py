def vcfExtract(vcf):
    import allel
    vcfInfo = allel.vcf_to_dataframe(vcf, ['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT'], alt_number=1)
    return vcfInfo
    #extract information from vcf to a df

def vcfDFtoList(vcfdf):
    vcfList = vcfdf.values.tolist()
    return vcfList
    #convert df to a list

def snvINDELlists(snvList, indelList, variantsList):
    for variant in variantsList:
        ref = len(str(variant[2]))
        alt = len(str(variant[3]))
        if (ref > 1 or alt > 1):
            indelList.append(variant)
        else:
            snvList.append(variant)
    #separate SNVs and INDELs into separate lists
    #this function takes (1) two empty SNV and INDEL lists and (2) a list with variants, and separates the variants according to size (SNVs and INDELs)

def variantCalls(truth, query):
    TPs = [x for x in query if x in truth]
    FPs = [x for x in query if x not in truth]
    FNs = [x for x in truth if x not in query]
    return TPs, FPs, FNs
