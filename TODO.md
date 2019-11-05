Todo list for VCFCompare
====

## Tasks

- [ ] 1. Parallelize the comparison for SNVs and INDELs to save time
- [ ] 2. Add png output for upset plots. The current output format is PDF
- [x] 3. Currently, the program assumes both the truth and query VCF files contain both SNVs and INDELs. If one or both (SNVs and INDELs) are not present in the two files, you will get a ZeroDivisionError
- [ ] 4. Write a function for handling no-ALT alleles
- [ ] 5. Write a function for genotype concordance
