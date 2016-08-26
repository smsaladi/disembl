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

import re
import ctypes
import warnings
from pkg_resources import resource_filename

from numpy.ctypeslib import load_library, ndpointer
import numpy as np
import pandas as pd
import scipy.signal
import Bio.SeqRecord

libdisembl = None
"""C library that does the heavy-lifting (loaded upon 1st use)
"""

re_remove = None
"""Initialize regex used to remove invalid characters
"""

default_params = {
    'handle_invalid': False,
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

def JensenNet(seq, handle_invalid):
    """Calculate hotloop, coil, and REM465 propensities

    Wraps predict_seq from libdisembl.c. Loads library upon first call.

    Parameters
    ----------
    seq : str
        Protein sequence for which propensities will be calculated

    handle_invalid : bool
        If `True`, attempt to handle invalid residues.

    Returns
    -------
    pd.DataFrame
        coils, hotloops, & rem465 propensities

    Raises
    ------
    None
    """

    if libdisembl is None:
        global libdisembl
        libdisembl = load_library("libdisembl",
                                  resource_filename(__name__, '.'))

        libdisembl.predict_seq.argtypes = [ctypes.c_char_p,
                                           ndpointer(dtype=np.float32),
                                           ndpointer(dtype=np.float32),
                                           ndpointer(dtype=np.float32)]
        libdisembl.predict_seq.restype = ctypes.c_void_p

        global re_remove
        re_remove = re.compile('[^FIVWMLCHYAGNRTPDEQSK]')

    # seq_len = len(seq)
    if len(seq) <= 20:
        # seq = "A" * 20 + seq + "A" * 20
        warnings.warn(("%s...%s: Sequence shorter than Neutral network window."
                      " Effectively padding with K") % (seq[:7], seq[-7:]))

    # handle invalid characters
    missing_pos = [match.span()[1]-1 for match in re_remove.finditer(seq)]
    seq = re_remove.sub('', seq)

    # create arrays for disembl output
    coils = np.empty(len(seq), dtype=np.float32)
    rem465 = np.empty_like(coils)
    hotloops = np.empty_like(coils)

    # http://stackoverflow.com/a/37888716/2320823
    libdisembl.predict_seq(str.encode(seq), rem465, hotloops, coils)

    if missing_pos and handle_invalid:
        warnings.warn("%s...%s: Unknown residues. Interpolating." %
                      (seq[:7], seq[-7:]))
        # If there are multiple missing, doing this in order will keep
        # indicies correct
        for pos in missing_pos:
            seq = seq[:pos] + 'X' + seq[pos:]
            try:
                # Average the values that are immediately surrounding
                coils = np.insert(coils, obj=pos,
                            values=np.mean((coils[pos-1], coils[pos])))
                rem465 = np.insert(rem465, obj=pos,
                            values=np.mean((rem465[pos-1], rem465[pos])))
                hotloops = np.insert(hotloops, obj=pos,
                            values=np.mean((hotloops[pos-1], hotloops[pos])))
            except IndexError as err:
                if pos == 0:
                    coils = np.insert(coils, obj=pos, values=coils[pos])
                    rem465 = np.insert(rem465, obj=pos, values=rem465[pos])
                    hotloops = np.insert(hotloops, obj=pos,
                                         values=hotloops[pos])
                elif pos == coils.size:
                    coils = np.append(coils, values=coils[-1])
                    rem465 = np.append(rem465, values=rem465[-1])
                    hotloops = np.append(hotloops, values=hotloops[-1])
                else:
                    print(pos, coils.size, flush=True)
                    raise IndexError(err)

    # Remove Alanine padding is performed
    # if seq_len <= 20:
    #     seq = seq[20:-20]
    #     coils = coils[20:-20]
    #     rem465 = rem465[20:-20]
    #     hotloops = hotloops[20:-20]

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
    print('\n', s, sep='')

    return


def calc_disembl_raw(protseq, handle_invalid, smooth_frame, calc_coils,
                     calc_hotloops, calc_rem465, old_filter):
    """Calculate the output of the DisEMBL neural network for a sequence

    Parameters
    ----------
    protseq : str
        Protein sequence for which to calculate DisEMBL scores

    handle_invalid : bool
        If `True`, attempt to handle invalid residues.

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

    pred = JensenNet(protseq.upper(), handle_invalid=handle_invalid)

    # For short sequences, shorten the window
    # For even length sequences, there is no smoothing peformed becuase
    # the window length encompassed the entire sequence
    # This seems arbitrary but is based on the original code/logic
    if len(pred.index) < smooth_frame*2+1:
        if len(pred.index) % 2 == 0:
            return pred
        else:
            smooth_frame = len(pred.index) // 2
    elif smooth_frame == 0:
        smooth_frame = 1

    to_calculate = []

    if calc_coils:
        to_calculate.append('coils')
    if calc_hotloops:
        to_calculate.append('hotloops')
    if calc_rem465:
        to_calculate.append('rem465')

    for name in to_calculate:
        if old_filter:
            temp_beg = np.copy(pred.loc[:, name].values[:smooth_frame])
            temp_end = np.copy(pred.loc[:, name].values[-smooth_frame:])

        pred[name] = scipy.signal.savgol_filter(pred[name],
                            window_length=smooth_frame*2+1,
                            polyorder=2, deriv=0, mode='interp')
        if old_filter:
            pred.loc[:, name].values[:smooth_frame] = temp_beg
            pred.loc[:, name].values[-smooth_frame:] = temp_end

    return pred


def calc_disembl(sequence, mode='summary', print_output=False,
                 handle_invalid=default_params['handle_invalid'],
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
                 old_filter=True):
    """Wrapper function to calculate the DisEMBL disorder predictions for a
    sequence of interest

    This is necessary to provide backwards compatibility such that predictions
    can be made from scripts/DisEMBL.py (installed into the path by `pip`) as
    well as by importing the module [new functionality].

    Below, XXX refers to the options that correspond to `coils`, `hotloops`,
    and `rem465`.

    Parameters
    ----------
    protseq : str or Bio.SeqRecord.SeqRecord
        Protein sequence for which to calculate DisEMBL scores

    mode : Optional[str] ('summary' or 'scores')
        Whether to return a summary of the disordered regions or the raw
        predictions.

    print_output : Optional[bool]
        If `True`, output will be printed to stdout.  If protseq is provided
        as a `str` and mode is `summary`, the sequence name will be the first
        and last 7 characters separated by dots.

    handle_invalid : Optional[bool]
        If `True`, attempt to handle invalid residues by interpolating
        neighbors. Passed eventually to `JensenNet`

    calc_XXX : Optional[bool]
        Whether to provide XXX output. Passed to `calc_disembl_raw`.

    smooth_frame : int
        Number of residues forward and backward to smooth other (window size
        would be smooth_frame * 2 + 1). Passed to `calc_disembl_raw`.

    fold_XXX : Optional[float]
        How much greater a peak value must be over its surrounding residues for
        a window to be identified as disordered. Passed to `getSlices`.

    expect_XXX : Optional[float]
        Minimum value used to identify regions of disorder. Passed to
        'getSlices`.

    join_frame : int
        If the distance between successive, putative regions is <= `join_frame`,
        the two are consolidated. Passed to `getSlices`.

    peak_frame : int
        After consolidating, if the length of a putative region is <
        `peak_frame`, the region is discarded. Passed to `getSlices`.

    old_filter : Optional[bool]
        The old filter (from TISEAN) simply replaces the first [smooth_frame]
        and last [smooth_frame] values with the raw calcuations. Passed to
        `calc_disembl_raw`.

    Returns
    -------
    pd.DataFrame
        Columns corresponding to the sequence as well as the properties
        specified (e.g. hotloops, coils)

    Raises
    ------
    None
    """

    if isinstance(sequence, Bio.SeqRecord.SeqRecord):
        seq_id = sequence.id
        sequence = str(sequence.seq)
    elif isinstance(sequence, str):
        seq_id = sequence[:7] + "..." + sequence[-7:]
    elif isinstance(sequence, Bio.Seq.Seq):
        ValueError("Pass protein sequence as Bio.SeqRecord.SeqRecord or string")
    else:
        ValueError("Unknown sequence type provided")

    preds = calc_disembl_raw(sequence, handle_invalid=handle_invalid,
                calc_coils=calc_coils, calc_hotloops=calc_coils,
                calc_rem465=calc_coils, smooth_frame=smooth_frame,
                old_filter=old_filter)

    if mode == 'scores':
        if print_output:
            preds.rename(columns=lambda x: x.upper(), inplace=True)
            preds = preds[['RESIDUE', 'COILS', 'REM465', 'HOTLOOPS']]
            print(preds.to_string(index=False))
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
            print()

        return slices
    else:
        raise ValueError("`mode` unrecognized: must be `scores` or `summary`")
