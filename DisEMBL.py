#!/usr/bin/env python
"""
This version of DisEMBL represents a complete rewriting of the code. The regions
identified as disordered do not change.

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

import fpformat
import sys
import tempfile
from os import system, popen3
from string import *

import os
import argparse

import pandas as pd
import scipy.signal
import Bio.SeqIO


def JensenNet(sequence, NN_bin='~/DisEMBL-1.4/disembl'):
    outFile = tempfile.mktemp()
    inFile = tempfile.mktemp()
    open(inFile,'w').write(sequence+'\n')
    system(NN_bin + '< ' + inFile +' > ' + outFile)
    REM465 = []
    COILS = []
    HOTLOOPS = []
    resultsFile = open(outFile,'r')
    results = resultsFile.readlines()
    resultsFile.close()
    for result in results:
        coil = float(fpformat.fix(split(result)[0],6))
        COILS.append(coil)
        hotloop = float(fpformat.fix(split(result)[1],6))
        HOTLOOPS.append(hotloop)
        rem465 = float(fpformat.fix(split(result)[2],6))
        REM465.append(rem465)
    os.remove(inFile)
    os.remove(outFile)
    return COILS, HOTLOOPS, REM465


def SavitzkyGolay(window, derivative, datalist, SG_bin='~/DisEMBL-1.4/sav_gol'):
    if len(datalist) < 2*window:
        window = len(datalist)/2
    elif window == 0:
        window = 1
    stdin, stdout, stderr = popen3(SG_bin + ' -V0 -D' + str(derivative) + ' -n' + str(window)+','+str(window))
    for data in datalist:
        stdin.write(`data`+'\n')
    try:
        stdin.close()
    except:
        print stderr.readlines()
    results = stdout.readlines()
    stdout.close()
    SG_results = []
    for result in results:
        f = float(fpformat.fix(result,6))
        if f < 0:
            SG_results.append(0)
        else:
            SG_results.append(f)
    return SG_results


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
            slices[i] = [ slices[i][0], slices[i+1][1] ]
            del slices[i+1]
        elif slices[i][1]-slices[i][0]+1 < peak_frame:
            del slices[i]
        else:
            i += 1
    return slices


def reportSlicesTXT(slices, sequence):
    if slices == []:
        s = lower(sequence)
    else:
        if slices[0][0] > 0:
            s = lower(sequence[0:slices[0][0]])
        else:
            s = ''
        for i in range(len(slices)):
            if i > 0:
                sys.stdout.write(', ')
            sys.stdout.write( str(slices[i][0]+1) + '-' + str(slices[i][1]+1) )
            s = s + upper(sequence[slices[i][0]:(slices[i][1]+1)])
            if i < len(slices)-1:
                s = s + lower(sequence[(slices[i][1]+1):(slices[i+1][0])])
            elif slices[i][1] < len(sequence)-1:
                s = s + lower(sequence[(slices[i][1]+1):(len(sequence))])
    print ''
    print s

    return


def runDisEMBLpipeline():
    try:
        smooth_frame = int(sys.argv[1])
        peak_frame = int(sys.argv[2])
        join_frame = int(sys.argv[3])
        fold_coils = float(sys.argv[4])
        fold_hotloops = float(sys.argv[5])
        fold_rem465 = float(sys.argv[6])
        file = str(sys.argv[7])
        try:
            mode = sys.argv[8]
        except:
            mode = 'default'
    except:
        print '\nDisEMBL.py smooth_frame peak_frame join_frame fold_coils fold_hotloops fold_rem465 sequence_file [mode]\n'
        print 'A default run would be: ./DisEMBL.py 8 8 4 1.2 1.4 1.2  fasta_file'
        print 'Mode: "default"(nothing) or "scores" which will give scores per residue in TAB seperated format'
        raise SystemExit
    db = open(file,'r')
    print ' ____  _     _____ __  __ ____  _       _  _  _'
    print '|  _ \(_)___| ____|  \/  | __ )| |     / || || |'
    print '| | | | / __|  _| | |\/| |  _ \| |     | || || |_'
    print '| |_| | \__ \ |___| |  | | |_) | |___  | ||__   _|'
    print '|____/|_|___/_____|_|  |_|____/|_____| |_(_) |_|'
    print '# Copyright (C) 2004 - Rune Linding & Lars Juhl Jensen '
    print '# EMBL Biocomputing Unit - Heidelberg - Germany        '
    print '#'
    for cur_record in Bio.SeqIO.parse(db, "fasta"):
	    sequence = upper(str(cur_record.seq))
            # Run NN
            COILS_raw, HOTLOOPS_raw, REM465_raw = JensenNet(sequence)
            # Run Savitzky-Golay
            REM465_smooth = SavitzkyGolay(smooth_frame,0,REM465_raw)
            COILS_smooth = SavitzkyGolay(smooth_frame,0,COILS_raw)
            HOTLOOPS_smooth = SavitzkyGolay(smooth_frame,0,HOTLOOPS_raw)
            if mode == 'default':
                sys.stdout.write('> '+cur_record.id+'_COILS ')
                reportSlicesTXT( getSlices(COILS_smooth, fold_coils, join_frame, peak_frame, 0.43), sequence )
                sys.stdout.write('> '+cur_record.id+'_REM465 ')
                reportSlicesTXT( getSlices(REM465_smooth, fold_rem465, join_frame, peak_frame, 0.50), sequence )
                sys.stdout.write('> '+cur_record.id+'_HOTLOOPS ')
                reportSlicesTXT( getSlices(HOTLOOPS_smooth, fold_hotloops, join_frame, peak_frame, 0.086), sequence )
                sys.stdout.write('\n')
            elif mode == 'scores':
                sys.stdout.write('# RESIDUE COILS REM465 HOTLOOPS\n')
                for i in range(len(REM465_smooth)):
                    sys.stdout.write(sequence[i]+'\t'+fpformat.fix(COILS_smooth[i],5)+'\t'+fpformat.fix(REM465_smooth[i],5)+'\t'+fpformat.fix(HOTLOOPS_smooth[i],5)+'\n')
            else:
                sys.stderr.write('Wrong mode given: '+mode+'\n')
                raise SystemExit
    db.close()

    return


def calc_disembl(protseq, smooth_frame,
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

    protseq = upper(protseq)
    COILS_raw, HOTLOOPS_raw, REM465_raw = JensenNet(protseq)

    pred = {'sequence': protseq}

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
        help='File format of [protein_seq]. Any Bio.SeqIO-readable format is supported')

    parser.add_argument('--smooth_frame',
        type=int,
        default=8,
        help='smooth_frame (window length = smooth_frame * 2 + 1)')

    parser.add_argument('--peak_frame',
        type=int,
        default=8,
        help='peak_frame')

    parser.add_argument('--join_frame',
        type=int,
        default=4,
        help='join_frame')

    parser.add_argument('--fold_coils',
        type=float,
        default=1.2,
        help='fold_coils')

    parser.add_argument('--fold_hotloops',
        type=float,
        default=1.4,
        help='fold_hotloops')

    parser.add_argument('--fold_rem465',
        type=float,
        default=1.2,
        help='fold_rem465')

    parser.add_argument('--expect_coils',
        type=float,
        default=0.43,
        help='expect_coils')

    parser.add_argument('--expect_hotloops',
        type=float,
        default=0.086,
        help='expect_hotloops')

    parser.add_argument('--expect_rem465',
        type=float,
        default=0.50,
        help='expect_rem465')

    parser.add_argument('--new_filter',
        action='store_true',
        help="Don't simply replace the first [smooth_frame] and last "
             "[smooth_frame] values. "
             "Use 'interp' from scipy.signal.savgol_filter")

    parser.add_argument('--mode',
        type=str,
        choices=['default', 'scores'],
        default='default',
        help='mode: default or scores which will give scores per residue'
             'in TAB seperated format')

    parser.add_argument('--quiet',
        action='store_true',
        help='Supress printing the DisEMBL banner')

    parser.add_argument('--old',
        action='store_true',
        help='Completely run the old code')

    args = parser.parse_args()

    if args.old:
        sys.argv = [sys.argv[0],
                    args.smooth_frame, args.peak_frame,
                    args.join_frame, args.fold_coils,
                    args.fold_hotloops, args.fold_rem465,
                    args.inputfile, args.mode]

        runDisEMBLpipeline()
        return


    if not args.quiet:
        print "\n".join(
             (' ____  _     _____ __  __ ____  _     ',
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
              '# ')) #, sep="\n", file=sys.stderr)

    for record in Bio.SeqIO.parse(args.inputfile, args.inputformat):
        sequence = upper(str(record.seq))

        preds = calc_disembl(sequence, args.smooth_frame)

        if args.mode == 'default':
            sys.stdout.write('> ' + record.id + '_COILS ')
            reportSlicesTXT(getSlices(preds.coils.tolist(),
                                      args.fold_coils,
                                      args.join_frame,
                                      args.peak_frame,
                                      args.expect_coils), sequence)
            sys.stdout.write('> ' + record.id + '_REM465 ')
            reportSlicesTXT(getSlices(preds.coils.tolist(),
                                      args.fold_rem465,
                                      args.join_frame,
                                      args.peak_frame,
                                      args.expect_rem465), sequence)
            sys.stdout.write('> ' + record.id + '_HOTLOOPS ')
            reportSlicesTXT(getSlices(preds.coils.tolist(),
                                      args.fold_hotloops,
                                      args.join_frame,
                                      args.peak_frame,
                                      args.expect_hotloops), sequence)
            print('\n')
        elif args.mode == 'scores':
            # Format header row as close to original DisEMBL as reasonable
            preds.rename(columns={"sequence": "residue"}, inplace=True)
            preds.columns = [upper(x) for x in preds.columns]
            print(preds.to_csv(sep="\t", index=False))
    return


if __name__ == '__main__':
    main()
