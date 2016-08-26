#!/usr/bin/env python

from itertools import islice
import re

import numpy as np
import pandas as pd
import Bio.SeqIO

import disembl

def test_raw_scores(old_output_fn, records_fn):
    records = Bio.SeqIO.parse(records_fn, 'fasta')
    old_output = open(old_output_fn, 'r+')

    # initial header/comment lines
    for _ in range(8):
        next(old_output, False)

    # columns to compare
    data_columns = ['COILS', 'REM465', 'HOTLOOPS']

    # first sequence record
    currecord = next(records, False)
    while(currecord):
        # +1 to account for the header/comment line
        print(len(currecord.seq))
        old = pd.DataFrame(list(islice(old_output, len(currecord.seq)+1))[1:])

        old = old[0].str.split(expand=True)
        try:
            old.columns = ['RESIDUE', 'COILS', 'REM465', 'HOTLOOPS']
        except ValueError:
            print(old)
            exit()

        # print(currecord.id)
        # print(old.head())
        print(old.tail(1))

        # calculation on new sequence
        new = disembl.calc_disembl(currecord, mode='scores')
        new.rename(columns=lambda x: x.upper(), inplace=True)

        # check that the sequences match in identity
        if not np.array_equal(old['RESIDUE'], new['RESIDUE']):
            print(currecord.id)
            def get_seq(x): return "".join(x['RESIDUE'].tolist())
            print("OLD SEQ:", get_seq(old))
            print("NEW SEQ:", get_seq(new))
            raise ValueError("reference and sequence do not match")

        # Check all outputs
        for col in data_columns:
            # `atol` becuase of the precision to which `old` was output
            if not np.allclose(pd.to_numeric(old[col]), new[col], atol=1e-04, rtol=1):
                print("OLD\n", old.tail())
                print("NEW\n", new.tail(), flush=True)
                raise AssertionError("`old` and `new` don't match for %s " % col)

        currecord = next(records, False)


if __name__ == '__main__':
    test_raw_scores("ecoli_k12.faa.disembl_v1.4.scores", "ecoli_k12.faa")
