#!/usr/bin/env python
"""
This version of DisEMBL represents a complete rewriting of the code. The
regions identified as disordered do not change.

Much of the calculation is done with NumPy/Scipy instead of previously written,
custom C code since it ends up being faster and significantly simplifies
maintainance.

The original calculation will be available via '--old' but will be removed in
further versions.

Authors:

Originally -
Copyright (C) 2004 Rune Linding & Lars Juhl Jensen - EMBL

Rewrite (2.0+) -
Copyright (C) 2016 Shyam Saladi - Caltech


The DisEMBL is licensed under the GNU General Public License
(http://www.opensource.org/licenses/gpl-license.php)
"""

import sys

import subprocess
import argparse
import io

import numpy as np
import pandas as pd
import scipy.signal
import Bio.SeqIO

default_params = {
    'smooth_frame': 8,
    'peak_frame': 8,
    'join_frame': 4,
    'fold_coils': 1.2,
    'fold_hotloops': 1.4,
    'fold_rem465': 1.2,
    'expect_coils': 0.43,
    'expect_hotloops': 0.086,
    'expect_rem465': 0.50
}


def JensenNet(sequence, NN_bin='/Users/saladi/disembl/disembl/disembl'):
    with subprocess.Popen([NN_bin],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          universal_newlines=True) as p:
        stdout_data = p.communicate(input=sequence + '\n')[0]

    df = pd.read_csv(io.StringIO(stdout_data), header=None, comment='%',
                     names=['COILS', 'HOTLOOPS', 'REM465'], delimiter='\t')

    return df.COILS.tolist(), df.HOTLOOPS.tolist(), df.REM465.tolist()


def getSlices(NNdata, fold, join_frame, peak_frame, expect_val):
    slices = []
    inSlice = 0
    for i in range(len(NNdata)):
        if inSlice:
            if NNdata[i] < expect_val:
                if maxSlice >= fold*expect_val:
                    slices.append([beginSlice, endSlice])
                inSlice = 0
            else:
                endSlice += 1
                if NNdata[i] > maxSlice:
                    maxSlice = NNdata[i]
        elif NNdata[i] >= expect_val:
            beginSlice = i
            endSlice = i
            inSlice = 1
            maxSlice = NNdata[i]
    if inSlice and maxSlice >= fold*expect_val:
        slices.append([beginSlice, endSlice])

    i = 0
    while i < len(slices):
        if i+1 < len(slices) and slices[i+1][0]-slices[i][1] <= join_frame:
            slices[i] = [slices[i][0], slices[i+1][1]]
            del slices[i+1]
        elif slices[i][1]-slices[i][0]+1 < peak_frame:
            del slices[i]
        else:
            i += 1
    return slices


def get_all_slices(coils, rem465, hotloops,
                   fold_coils=default_params['fold_coils'],
                   expect_coils=default_params['expect_coils'],
                   fold_rem465=default_params['fold_rem465'],
                   expect_rem465=default_params['expect_rem465'],
                   fold_hotloops=default_params['fold_hotloops'],
                   expect_hotloops=default_params['expect_hotloops'],
                   join_frame=default_params['join_frame'],
                   peak_frame=default_params['peak_frame']):
        return {
            'coils': getSlices(coils, fold_coils, join_frame,
                               peak_frame, expect_coils),
            'rem465': getSlices(rem465, fold_rem465, join_frame,
                                peak_frame, expect_rem465),
            'hotloops': getSlices(hotloops, fold_hotloops, join_frame,
                                  peak_frame, expect_hotloops)
            }


def reportSlicesTXT(slices, sequence):
    if slices == []:
        s = sequence.lower()
    else:
        if slices[0][0] > 0:
            s = sequence[0:slices[0][0]].lower()
        else:
            s = ''
        for i in range(len(slices)):
            if i > 0:
                sys.stdout.write(', ')
            sys.stdout.write(str(slices[i][0]+1) + '-' + str(slices[i][1]+1))
            s = s + sequence[slices[i][0]:(slices[i][1]+1)].upper()
            if i < len(slices)-1:
                s = s + sequence[(slices[i][1]+1):(slices[i+1][0])].lower()
            elif slices[i][1] < len(sequence)-1:
                s = s + sequence[(slices[i][1]+1):(len(sequence))].lower()
    print('')
    print(s)

    return


def calc_disembl(protseq, smooth_frame=default_params['smooth_frame'],
                 coils=True, hotloops=True, rem465=True, old_filter=False):
    """Calculate the output of the DisEMBL neural network for a sequence

    Parameters
    ----------
    protseq : str
        Protein sequence for which to calculate DisEMBL scores

    smooth_frame : int,
        Number of residues forward and backward to smooth other (window size
        would be smooth_frame * 2 + 1)

    coils : Optional[bool]
        Whether to provide coils output

    hotloops : Optional[bool]
        Whether to provide hotloops output

    rem465 : Optional[bool]
        Whether to provide rem465 output

    old_filter : Optional[bool]
        The old filter (from TISEAN) simply replaces the first [smooth_frame]
            and last [smooth_frame] values with the raw calcuations

    Returns
    -------
    pd.DataFrame
        Columns corresponding to the sequence as well as the properties
        specified (e.g. hotloops, coils)

    Raises
    ------
    None
    """

    protseq = protseq.upper()
    COILS_raw, HOTLOOPS_raw, REM465_raw = JensenNet(protseq)

    pred = {'sequence': list(protseq)}

    if coils:
        pred['coils'] = scipy.signal.savgol_filter(COILS_raw,
                            window_length=smooth_frame*2+1,
                            polyorder=2, deriv=0, mode='interp')
        if old_filter:
            pred['coils'][:smooth_frame] = REM465_raw[:smooth_frame]
            pred['coils'][-smooth_frame:] = REM465_raw[-smooth_frame:]

    if hotloops:
        pred['hotloops'] = scipy.signal.savgol_filter(HOTLOOPS_raw,
                            window_length=smooth_frame*2+1,
                            polyorder=2, deriv=0, mode='interp')
        if old_filter:
            pred['hotloops'][:smooth_frame] = REM465_raw[:smooth_frame]
            pred['hotloops'][-smooth_frame:] = REM465_raw[-smooth_frame:]

    if rem465:
        pred['rem465'] = scipy.signal.savgol_filter(REM465_raw,
                            window_length=smooth_frame*2+1,
                            polyorder=2, deriv=0, mode='interp')
        if old_filter:
            pred['rem465'][:smooth_frame] = REM465_raw[:smooth_frame]
            pred['rem465'][-smooth_frame:] = REM465_raw[-smooth_frame:]

    return pd.DataFrame(pred)


def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser(
        description='Calculate DisEMBL disorder prediction.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('inputfile',
                        type=str,
                        default=sys.stdin,
                        help='File of protein coding sequences to predict on.')

    parser.add_argument('--inputformat',
                        type=str,
                        default='fasta',
                        help='File format of [protein_seq].'
                             'Any Bio.SeqIO-readable supported')

    parser.add_argument('--smooth_frame',
                        type=int,
                        default=default_params['smooth_frame'],
                        help='smooth_frame (window = smooth_frame * 2 + 1)')

    parser.add_argument('--peak_frame',
                        type=int,
                        default=default_params['peak_frame'],
                        help='peak_frame')

    parser.add_argument('--join_frame',
                        type=int,
                        default=default_params['join_frame'],
                        help='join_frame')

    parser.add_argument('--fold_coils',
                        type=float,
                        default=default_params['fold_coils'],
                        help='fold_coils')

    parser.add_argument('--fold_hotloops',
                        type=float,
                        default=default_params['fold_hotloops'],
                        help='fold_hotloops')

    parser.add_argument('--fold_rem465',
                        type=float,
                        default=default_params['fold_rem465'],
                        help='fold_rem465')

    parser.add_argument('--expect_coils',
                        type=float,
                        default=default_params['expect_coils'],
                        help='expect_coils')

    parser.add_argument('--expect_hotloops',
                        type=float,
                        default=default_params['expect_hotloops'],
                        help='expect_hotloops')

    parser.add_argument('--expect_rem465',
                        type=float,
                        default=default_params['expect_rem465'],
                        help='expect_rem465')

    parser.add_argument('--old_filter',
                        action='store_true',
                        help="Simply replace the first [smooth_frame] and last"
                             " [smooth_frame] values, instead of using "
                             "Use 'interp' from scipy.signal.savgol_filter")

    parser.add_argument('--mode',
                        type=str,
                        choices=['default', 'scores'],
                        default='default',
                        help='mode: default or scores which will give scores'
                             'per residue in TAB seperated format')

    parser.add_argument('--quiet',
        action='store_true',
        help='Supress printing the DisEMBL banner')

    args = parser.parse_args()

    if not args.quiet:
        print(' ____  _     _____ __  __ ____  _     ',
              '|  _ \(_)___| ____|  \/  | __ )| |    ',
              '| | | | / __|  _| | |\/| |  _ \| |    ',
              '| |_| | \__ \ |___| |  | | |_) | |___ ',
              '|____/|_|___/_____|_|  |_|____/|_____|',
              '# ',
              '# Original:',
              '# Copyright (C) 2004 - Rune Linding & Lars Juhl Jensen',
              '# EMBL Biocomputing Unit - Heidelberg - Germany',
              '# ',
              '# Rewrite (v2.0+):',
              '# Copyright (C) 2016 - Shyam Saladi',
              '# California Institute of Technology - Pasadena - CA - USA',
              '# ', sep="\n", file=sys.stderr)

    for record in Bio.SeqIO.parse(args.inputfile, args.inputformat):
        sequence = str(record.seq).upper()

        preds = calc_disembl(sequence,
                             smooth_frame=args.smooth_frame,
                             old_filter=args.old_filter)

        if args.mode == 'default':
            slices = get_all_slices(coils=preds.coils.tolist(),
                                    rem465=preds.rem465.tolist(),
                                    hotloops=preds.hotloops.tolist(),
                                    fold_coils=args.fold_coils,
                                    expect_coils=args.expect_coils,
                                    fold_rem465=args.fold_rem465,
                                    expect_rem465=args.expect_rem465,
                                    fold_hotloops=args.fold_hotloops,
                                    expect_hotloops=args.expect_hotloops,
                                    join_frame=args.join_frame,
                                    peak_frame=args.peak_frame)
            sys.stdout.write('> ' + record.id + '_COILS ')
            reportSlicesTXT(slices['coils'], sequence)
            sys.stdout.write('> ' + record.id + '_REM465 ')
            reportSlicesTXT(slices['rem465'], sequence)
            sys.stdout.write('> ' + record.id + '_HOTLOOPS ')
            reportSlicesTXT(slices['hotloops'], sequence)
            print('\n')
        elif args.mode == 'scores':
            # Format header row as close to original DisEMBL as reasonable
            preds.rename(columns={"sequence": "residue"}, inplace=True)
            preds.columns = [x.upper() for x in preds.columns]
            print(preds.to_csv(sep=" ", index=False))
    return


if __name__ == '__main__':
    main()
