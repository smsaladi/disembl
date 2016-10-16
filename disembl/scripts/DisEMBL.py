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

import Bio.SeqIO
import disembl

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
                        default=disembl.default_params['smooth_frame'],
                        help='smooth_frame (window = smooth_frame * 2 + 1)')

    parser.add_argument('--peak_frame',
                        type=int,
                        default=disembl.default_params['peak_frame'],
                        help='peak_frame')

    parser.add_argument('--join_frame',
                        type=int,
                        default=disembl.default_params['join_frame'],
                        help='join_frame')

    parser.add_argument('--fold_coils',
                        type=float,
                        default=disembl.default_params['fold_coils'],
                        help='fold_coils')

    parser.add_argument('--fold_hotloops',
                        type=float,
                        default=disembl.default_params['fold_hotloops'],
                        help='fold_hotloops')

    parser.add_argument('--fold_rem465',
                        type=float,
                        default=disembl.default_params['fold_rem465'],
                        help='fold_rem465')

    parser.add_argument('--expect_coils',
                        type=float,
                        default=disembl.default_params['expect_coils'],
                        help='expect_coils')

    parser.add_argument('--expect_hotloops',
                        type=float,
                        default=disembl.default_params['expect_hotloops'],
                        help='expect_hotloops')

    parser.add_argument('--expect_rem465',
                        type=float,
                        default=disembl.default_params['expect_rem465'],
                        help='expect_rem465')

    parser.add_argument('--old_invalid',
                        action='store_true',
                        default=not disembl.default_params['handle_invalid'],
                        help="Don't attempt to handle invalid residues by"
                             "interpolating neighbors.")

    parser.add_argument('--old_filter',
                        action='store_true',
                        help="Simply replace the first [smooth_frame] and last"
                             " [smooth_frame] values, instead of using "
                             "Use 'interp' from scipy.signal.savgol_filter")

    parser.add_argument('--mode',
                        type=str,
                        choices=['summary', 'scores'],
                        default='summary',
                        help='Default (regions identified with upper/lowercase'
                             'sequences or raw scores in tab delimited format')

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
              '# Copyright (C) 2016 - Nadine Bradbury & Shyam Saladi',
              '# California Institute of Technology - Pasadena - CA - USA',
              '# ', sep="\n", file=sys.stderr)

    if(args.inputfile == '-'):
        args.inputfile = sys.stdin

    for record in Bio.SeqIO.parse(args.inputfile, args.inputformat):
        record.seq = record.seq.upper()

        disembl.calc_disembl(record, mode=args.mode, print_output=True,
                             handle_invalid=not args.old_invalid,
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
