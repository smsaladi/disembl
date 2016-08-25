#!/usr/bin/env python
"""
 ____  _     _____ __  __ ____  _        ____      ______
|  _ \(_)___| ____|  \/  | __ )| |      /___  )   (  __  )
| | | | / __|  _| | |\/| |  _ \| |         / /    | |  | |
| |_| | \__ \ |___| |  | | |_) | |___     / /_    | |__| |
|____/|_|___/_____|_|  |_|____/|_____|  /_____|(_)(______)

Originally -
Copyright (C) 2004 Rune Linding & Lars Juhl Jensen - EMBL

Rewrite (2.0+) -
Copyright (C) 2016 Nadine Bradbury & Shyam Saladi - Caltech

The DisEMBL is licensed under the GNU General Public License
(http://www.opensource.org/licenses/gpl-license.php)
"""

import sys
import argparse
import ctypes
from pkg_resources import resource_filename

from numpy.ctypeslib import load_library, ndpointer
import numpy as np
import pandas as pd
import scipy.signal
import Bio.SeqIO
import Bio.SeqRecord

libdisembl = load_library("libdisembl", resource_filename(__name__, '.'))
"""Load C library that does the heavy-lifting of the neural network
"""
libdisembl.predict_seq.argtypes = [ctypes.c_char_p,
                                   ndpointer(dtype=np.float32),
                                   ndpointer(dtype=np.float32),
                                   ndpointer(dtype=np.float32)]
libdisembl.predict_seq.restype = ctypes.c_void_p

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
"""Default parameters for identifying regions of disorder

Specified here to facilitate using as default parameters both through argparse
as well as through direct calls to calc_disembl.
"""

def JensenNet(seq):
    """Calculate hotloop, coil, and REM465 propensities

    Wraps predict_seq from libdisembl.c

    Parameters
    ----------
    seq : str
        Protein sequence for which propensities will be calculated

    Returns
    -------
    pd.DataFrame
        coils, hotloops, & rem465 propensities

    Raises
    ------
    None
    """
    # create arrays for disembl output
    coils = np.zeros(len(seq), dtype=np.float32)
    rem465 = np.zeros_like(coils)
    hotloops = np.zeros_like(coils)

    # http://stackoverflow.com/a/37888716/2320823
    libdisembl.predict_seq(str.encode(seq), rem465, hotloops, coils)

    return pd.DataFrame({'residue': list(seq),
                         'coils': coils,
                         'rem465': rem465,
                         'hotloops': hotloops})

def getSlices(NNdata, fold, join_frame, peak_frame, expect_val):
    slices = []

    inSlice = False
    beginSlice = 0
    endSlice = 0
    maxSlice = NNdata[0]

    ## here we find windows with all values > expect_val and at least one
    ## value >= fold*expect_val, storing the first and last residue positions
    ## in slices [(,),(,),(,)]

    for i, nn_value in enumerate(NNdata):
        if inSlice:
            if nn_value < expect_val:
                if maxSlice >= fold*expect_val:
                    slices.append([beginSlice, endSlice])
                inSlice = 0
            else:
                endSlice += 1
                if nn_value > maxSlice:
                    maxSlice = nn_value
        elif nn_value >= expect_val:
            beginSlice = i
            endSlice = i
            inSlice = True
            maxSlice = nn_value
    if inSlice and maxSlice >= fold*expect_val:
        slices.append([beginSlice, endSlice])

    ## if the distance between successive windows is <= join_frame,
    ## consolidate the two frames into a single one
    ## after consolidating, if the length of a window is < peak_frame
    ## then delete the window

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


def reportSlicesTXT(slices, sequence):
    if slices == []:
        s = sequence.lower()
    else:
        if slices[0][0] > 0:
            s = sequence[0:slices[0][0]].lower()
        else:
            s = ''
        for i, _ in enumerate(slices):
            if i > 0:
                print(', ', end='')
            print(str(slices[i][0]+1) + '-' + str(slices[i][1]+1), end='')
            s = s + sequence[slices[i][0]:(slices[i][1]+1)].upper()
            if i < len(slices)-1:
                s = s + sequence[(slices[i][1]+1):(slices[i+1][0])].lower()
            elif slices[i][1] < len(sequence)-1:
                s = s + sequence[(slices[i][1]+1):(len(sequence))].lower()
    print('\n', s)

    return


def calc_disembl_raw(protseq, smooth_frame, calc_coils, calc_hotloops,
                     calc_rem465, old_filter):
    """Calculate the output of the DisEMBL neural network for a sequence

    Parameters
    ----------
    protseq : str
        Protein sequence for which to calculate DisEMBL scores

    smooth_frame : int,
        Number of residues forward and backward to smooth other (window size
        would be smooth_frame * 2 + 1)

    calc_coils : Optional[bool]
        Whether to provide coils output

    calc_hotloops : Optional[bool]
        Whether to provide hotloops output

    calc_rem465 : Optional[bool]
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

    pred = JensenNet(protseq.upper())

    if calc_coils:
        if old_filter:
            temp_beg = np.copy(pred['coils'].values[:smooth_frame])
            temp_end = np.copy(pred['coils'].values[-smooth_frame:])

        pred['coils'] = scipy.signal.savgol_filter(pred['coils'],
                            window_length=smooth_frame*2+1,
                            polyorder=2, deriv=0, mode='interp')
        if old_filter:
            pred['coils'][:smooth_frame] = temp_beg
            pred['coils'][-smooth_frame:] = temp_end

    if calc_hotloops:
        if old_filter:
            temp_beg = np.copy(pred['hotloops'].values[:smooth_frame])
            temp_end = np.copy(pred['hotloops'].values[-smooth_frame:])

        pred['hotloops'] = scipy.signal.savgol_filter(pred['hotloops'],
                            window_length=smooth_frame*2+1,
                            polyorder=2, deriv=0, mode='interp')
        if old_filter:
            pred['hotloops'][:smooth_frame] = temp_beg
            pred['hotloops'][-smooth_frame:] = temp_end

    if calc_rem465:
        if old_filter:
            temp_beg = np.copy(pred['rem465'].values[:smooth_frame])
            temp_end = np.copy(pred['rem465'].values[-smooth_frame:])

        pred['rem465'] = scipy.signal.savgol_filter(pred['rem465'],
                            window_length=smooth_frame*2+1,
                            polyorder=2, deriv=0, mode='interp')
        if old_filter:
            pred['rem465'][:smooth_frame] = temp_beg
            pred['rem465'][-smooth_frame:] = temp_end

    return pred


def calc_disembl(sequence, mode='summary', print_output=False,
                 calc_coils=True, calc_hotloops=True, calc_rem465=True,
                 smooth_frame=default_params['smooth_frame'],
                 fold_coils=default_params['fold_coils'],
                 expect_coils=default_params['expect_coils'],
                 fold_rem465=default_params['fold_rem465'],
                 expect_rem465=default_params['expect_rem465'],
                 fold_hotloops=default_params['fold_hotloops'],
                 expect_hotloops=default_params['expect_hotloops'],
                 join_frame=default_params['join_frame'],
                 peak_frame=default_params['peak_frame'],
                 old_filter=False):
    """
    summary or scores
    """
    if isinstance(sequence, Bio.SeqRecord.SeqRecord):
        preds = calc_disembl_raw(sequence.seq, smooth_frame=smooth_frame,
                    calc_coils=calc_coils, calc_hotloops=calc_coils,
                    calc_rem465=calc_coils, old_filter=old_filter)
        seq_id = sequence.id
    else:
        preds = calc_disembl_raw(sequence, smooth_frame=smooth_frame,
                    calc_coils=calc_coils, calc_hotloops=calc_coils,
                    calc_rem465=calc_coils, old_filter=old_filter)
        seq_id = 'disembl_seq'

    if mode == 'scores':
        if print_output:
            print(preds.to_string())
        return preds
    elif mode == 'summary':
        slices = {
            'coils': getSlices(preds.coils, fold_coils, join_frame,
                               peak_frame, expect_coils),
            'rem465': getSlices(preds.rem465, fold_rem465, join_frame,
                                peak_frame, expect_rem465),
            'hotloops': getSlices(preds.hotloops, fold_hotloops, join_frame,
                                  peak_frame, expect_hotloops)
            }

        if print_output:
            print('> ', seq_id, '_COILS ', sep='', end='')
            reportSlicesTXT(slices['coils'], sequence)
            print('> ', seq_id, '_REM465 ', sep='', end='')
            reportSlicesTXT(slices['rem465'], sequence)
            print('> ', seq_id, '_HOTLOOPS ', sep='', end='')
            reportSlicesTXT(slices['hotloops'], sequence)
            print('\n')

        return slices

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

    if args.mode == 'default':
        mode = 'summary'
    elif args.mode == 'scores':
        mode = 'scores'
    else:
        raise ValueError("Mode unrecognized. See help.")

    for record in Bio.SeqIO.parse(args.inputfile, args.inputformat):
        record.seq = record.seq.upper()

        calc_disembl(record, mode=mode, print=True,
                     smooth_frame=args.smooth_frame,
                     fold_coils=args.fold_coils,
                     expect_coils=args.expect_coils,
                     fold_rem465=args.fold_rem465,
                     expect_rem465=args.expect_rem465,
                     fold_hotloops=args.fold_hotloops,
                     expect_hotloops=args.expect_hotloops,
                     join_frame=args.join_frame,
                     peak_frame=args.peak_frame,
                     old_filter=args.old_filter)

    return

if __name__ == '__main__':
    main()
