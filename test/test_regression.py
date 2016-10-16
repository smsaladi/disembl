#!/usr/bin/env python

from itertools import islice
import re
import subprocess

import numpy as np
import pandas as pd
import Bio.SeqIO

import disembl

def raw_score_compare(old_output_fn, records_fn):
    """
    """
    records = Bio.SeqIO.parse(records_fn, 'fasta')
    old_output = open(old_output_fn, 'r+')

    re_remove = re.compile('[^FIVWMLCHYAGNRTPDEQSK]')

    # initial header/comment lines
    for _ in range(8):
        next(old_output, False)

    # columns to compare
    data_columns = ['COILS', 'REM465', 'HOTLOOPS']

    # first sequence record
    currecord = next(records, False)
    while(currecord):
        # +1 to account for the header/comment line
        seq_len = len(currecord.seq)
        currecord.seq = re_remove.sub('', str(currecord.seq))

        old = pd.DataFrame(list(islice(old_output, len(currecord.seq)+1))[1:])

        old = old[0].str.split(expand=True)
        try:
            old.columns = ['RESIDUE', 'COILS', 'REM465', 'HOTLOOPS']
        except ValueError:
            print(old)
            exit()

        # calculation on new sequence
        try:
            new = disembl.calc_disembl(currecord, mode='scores')
        except TypeError as err:
            print(err, currecord.seq[10:])
            exit()
        new.rename(columns=lambda x: x.upper(), inplace=True)

        # check that the sequences match in identity
        # if there were undetermined characters, the sequences will not match
        if not np.array_equal(old['RESIDUE'], new['RESIDUE']) and \
                len(currecord.seq) == seq_len:
            print(currecord.id)
            def get_seq(x): return "".join(x['RESIDUE'].tolist())
            print("OLD SEQ:", get_seq(old))
            print("NEW SEQ:", get_seq(new))
            raise ValueError("reference and sequence do not match")

        # Check all outputs
        for col in data_columns:
            # `atol` becuase of the precision to which `old` was printed
            if not np.allclose(pd.to_numeric(old[col]), new[col], atol=1e-04, rtol=1):
                print("OLD\n", old.tail())
                print("NEW\n", new.tail(), flush=True)
                raise AssertionError("`old` and `new` don't match for %s " % col)

        currecord = next(records, False)

    return

def disembl_fmt_compare(old_output_fn, records_fn):
    """
    """
    old_output = open(old_output_fn, 'r+')

    oldline = next(old_output)
    while oldline[0] != '>':
        oldline = next(old_output)

    with subprocess.Popen(['DisEMBL.py', records_fn, '--old_filter',
                           '--old_invalid'],
                         stdin = subprocess.PIPE,
                         stdout = subprocess.PIPE,
                         universal_newlines=True) as p:
        data_start = False
        for newline in p.communicate(input='')[0].split('\n'):
            if not newline.strip() or newline[0] == '>' or data_start:
                data_start = True
                if newline.rstrip() != oldline.rstrip():
                    print("ERR")
                    print("NEW", newline)
                    print("OLD", oldline)
                    raise AssertionError("Printed lines do not match")
                try:
                    oldline = next(old_output)
                except StopIteration:
                    return
    return

def test_ecolik12_raw():
    raw_score_compare("test/ecoli_k12.faa.disembl_v1.4.scores", "test/ecoli_k12.faa")
    return

def test_ecolik12_formatted():
    disembl_fmt_compare("test/ecoli_k12.faa.disembl_v1.4", "test/ecoli_k12.faa")
    return

def test_secy_raw():
    raw_score_compare("test/secy.faa.disembl_v1.4.scores", "test/secy.faa")
    return

def test_secy_formatted():
    disembl_fmt_compare("test/secy.faa.disembl_v1.4", "test/secy.faa")
    return
