## Regression tests

*E. coli* K-12 proteome downloaded from Uniprot using the following query

> organism:"escherichia coli strain k12" AND reviewed:yes AND organism:"Escherichia coli (strain K12) [83333]"

No particular reason to use this dataset, but it presumably covers many
cases to regression test the module and any future changes.

From DisEMBL v1.4, output files were generated via the following commands

```bash
python DisEMBL.py 8 8 4 1.2 1.4 1.2 test/ecolik12.faa scores > test/ecolik12.faa.disembl_v1.4.scores
python DisEMBL.py 8 8 4 1.2 1.4 1.2 test/ecolik12.faa > test/ecolik12.faa.disembl_v1.4

python DisEMBL.py 8 8 4 1.2 1.4 1.2 test/secy.faa scores > test/secy.faa.disembl_v1.4.scores
python DisEMBL.py 8 8 4 1.2 1.4 1.2 test/secy.faa > test/secy.faa.disembl_v1.4
```
