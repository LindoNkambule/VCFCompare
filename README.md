# VCFCompare
A Python program for evaluating site-level concordance of a query VCF against a truth VCF.

This tool compares two variant callsets against each other and produces a CSV file with summary metrics. The summary metric CSV file contains:
* Variant type: SNV or INDEL
* Total number of variants in the truth and query VCF files
* Total true-positive, false-positive, and false-negative calls
* Recall and Precision

## Usage
```
python3 VCFCompare.py --truth truth.vcf --query query.vcf --output output
```
