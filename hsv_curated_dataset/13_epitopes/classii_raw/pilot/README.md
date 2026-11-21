# Pilot submission — NOT part of the analysis

`905068_NetMHCIIpan.xls` is the format-validation run: HSV-1 gD (US6) only, 4 alleles
(DRB1_0101, DRB1_1501, DRB3_0301, HLA-DQA10501-DQB10201).

It is kept here, in a subdirectory, for two reasons:

1. It is the record of how the output format was verified before the full runs.
2. It must NOT be globbed into the analysis. It covers one protein instead of 40, its
   allele set is not the IEDB reference panel (DRB3_0301 is not in the panel), and two of
   its alleles also appear in run 1 — combining it with the full runs would duplicate
   those peptide/allele pairs.

`step29_classii_predictions.py` reads `13_epitopes/classii_raw/*.xls`, which does not
recurse, so this file is excluded by location.
