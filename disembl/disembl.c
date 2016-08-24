/***************************************************************/
/*  ____  _     _____ __  __ ____  _        ____      _____   */
/* |  _ \(_)___| ____|  \/  | __ )| |      /___  )   (  __  ) */
/* | | | | / __|  _| | |\/| |  _ \| |         / /    | |  | | */
/* | |_| | \__ \ |___| |  | | |_) | |___     / /_    | |__| | */
/* |____/|_|___/_____|_|  |_|____/|_____|  /_____|(_)(______) */
/* DisEMBL is Copyright (C) 2003                               */
/* Lars Juhl Jensen & Rune Linding - EMBL                      */
/* Licensed under GPL                                          */
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Uncomment to have scores printed to standard output */
#define OLD

/* Define size of the alphabet */
#define NA 21

/* Define max number of neurons */
#define MW 41
#define MH 30

#include "russel.h"
#include "bfactor.h"
#include "missing.h"

static float sigmoid[256] = {
  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
  0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000001, 0.000001, 0.000001,
  0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000002, 0.000002, 0.000002,
  0.000002, 0.000003, 0.000003, 0.000003, 0.000004, 0.000004, 0.000005, 0.000005,
  0.000006, 0.000007, 0.000008, 0.000009, 0.000010, 0.000011, 0.000013, 0.000015,
  0.000017, 0.000019, 0.000021, 0.000024, 0.000028, 0.000031, 0.000035, 0.000040,
  0.000045, 0.000051, 0.000058, 0.000066, 0.000075, 0.000085, 0.000096, 0.000109,
  0.000123, 0.000140, 0.000158, 0.000180, 0.000203, 0.000231, 0.000261, 0.000296,
  0.000335, 0.000380, 0.000431, 0.000488, 0.000553, 0.000626, 0.000710, 0.000804,
  0.000911, 0.001032, 0.001170, 0.001325, 0.001501, 0.001701, 0.001927, 0.002183,
  0.002473, 0.002801, 0.003173, 0.003594, 0.004070, 0.004610, 0.005220, 0.005911,
  0.006693, 0.007577, 0.008577, 0.009708, 0.010987, 0.012432, 0.014064, 0.015906,
  0.017986, 0.020332, 0.022977, 0.025957, 0.029312, 0.033086, 0.037327, 0.042088,
  0.047426, 0.053403, 0.060087, 0.067547, 0.075858, 0.085099, 0.095349, 0.106691,
  0.119203, 0.132964, 0.148047, 0.164516, 0.182426, 0.201813, 0.222700, 0.245085,
  0.268941, 0.294215, 0.320821, 0.348645, 0.377541, 0.407333, 0.437823, 0.468791,
  0.500000, 0.531209, 0.562177, 0.592667, 0.622459, 0.651355, 0.679179, 0.705785,
  0.731059, 0.754915, 0.777300, 0.798187, 0.817574, 0.835484, 0.851953, 0.867036,
  0.880797, 0.893309, 0.904651, 0.914901, 0.924142, 0.932453, 0.939913, 0.946597,
  0.952574, 0.957912, 0.962673, 0.966914, 0.970688, 0.974043, 0.977023, 0.979668,
  0.982014, 0.984094, 0.985936, 0.987568, 0.989013, 0.990292, 0.991423, 0.992423,
  0.993307, 0.994089, 0.994780, 0.995390, 0.995930, 0.996406, 0.996827, 0.997199,
  0.997527, 0.997817, 0.998073, 0.998299, 0.998499, 0.998675, 0.998830, 0.998968,
  0.999089, 0.999196, 0.999290, 0.999374, 0.999447, 0.999512, 0.999569, 0.999620,
  0.999665, 0.999704, 0.999739, 0.999769, 0.999797, 0.999820, 0.999842, 0.999860,
  0.999877, 0.999891, 0.999904, 0.999915, 0.999925, 0.999934, 0.999942, 0.999949,
  0.999955, 0.999960, 0.999965, 0.999969, 0.999972, 0.999976, 0.999979, 0.999981,
  0.999983, 0.999985, 0.999987, 0.999989, 0.999990, 0.999991, 0.999992, 0.999993,
  0.999994, 0.999995, 0.999995, 0.999996, 0.999996, 0.999997, 0.999997, 0.999997,
  0.999998, 0.999998, 0.999998, 0.999998, 0.999999, 0.999999, 0.999999, 0.999999,
  0.999999, 0.999999, 0.999999, 0.999999, 0.999999, 1.000000, 1.000000, 1.000000,
  1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000
};


/*
 *  Perform feed forward evaluation of a neural network on a sequence window.
 */
static float
feed_forward(int const *s, float const w[], int nw, int nh) {

    float h[MH], o[2], x;
    int i, j;

    /* Shift input window to match network window size */
    s += (MW-nw)/2;

    /* Feed input values to hidden layer making use of sparse encoding */
    for (i = 0; i < nh; ++i) {
        x = w[(NA*nw+1)*(i+1)-1];

        for (j = 0; j < nw; ++j)
            x += w[(NA*nw+1)*i+NA*j+s[j]];

        if (x <= -16)
            h[i] = 0;
        else if (x >= 16)
            h[i] = 1;
        else
            h[i] = sigmoid[(int)(8*x+128)];
    }

    /* Feed hidden layer values to output layer */
    for (i = 0; i <= 1; ++i) {
        x = w[(NA*nw+1)*nh+(nh+1)*(i+1)-1];

        for (j = 0; j < nh; ++j)
            x += w[(NA*nw+1)*nh+(nh+1)*i+j]*h[j];

        if (x <= -16)
            o[i] = 0;
        else if (x >= 16)
            o[i] = 1;
        else
            o[i] = sigmoid[(int)(8*x+128)];
    }

    /* Combine the scores from the two output neurons */
    return((o[0]+1-o[1])/2);
}


/*
 *  Calculate scores for a sequence window.
 */
static void
predict(int const *s, float *sm, float *sb, float *sr) {

    *sr = feed_forward(s, r19_1, 19, 30) +
       feed_forward(s, r19_2, 19, 30) +
       feed_forward(s, r19_3, 19, 30) +
       feed_forward(s, r19_4, 19, 30) +
       feed_forward(s, r19_5, 19, 30);
    *sr = 0.07387214+0.8020778*(*sr)/5;

    *sb = feed_forward(s, b41_1, 41, 5) +
       feed_forward(s, b41_2, 41, 5) +
       feed_forward(s, b41_3, 41, 5) +
       feed_forward(s, b41_4, 41, 5) +
       feed_forward(s, b41_5, 41, 5);
    *sb = 0.08016882+0.6282424*(*sb)/5;
    *sb *= *sr;

    *sm = feed_forward(s, m9_1, 9, 30) +
       feed_forward(s, m9_2, 9, 30) +
       feed_forward(s, m9_3, 9, 30) +
       feed_forward(s, m9_4, 9, 30) +
       feed_forward(s, m9_5, 9, 30) +
       feed_forward(s, m21_1, 21, 30) +
       feed_forward(s, m21_2, 21, 30) +
       feed_forward(s, m21_3, 21, 30) +
       feed_forward(s, m21_4, 21, 30) +
       feed_forward(s, m21_5, 21, 30);
    *sm /= 10;

    return;
}


/*
 *  Calculate scores for an entire sequence
 */
void
predict_seq(char const *seq, float *sm_arr, float *sb_arr, float *sr_arr) {
    // This ordering specifies how characters (residues)
    // are turned into numbers to be interpreted by the
    // neural network (Shyam)
    char *alphabet = "FIVWMLCHYAGNRTPDEQSK";

    char *p;
    int i, j, s[MW];

    // Fill with default values
    for (i = 0; i < MW; ++i)
        s[i] = NA-1;

    if ((int)strlen(seq) < (MW-1)/2) {
        printf("Sequence must be at least %d residues long", (MW-1)/2);
        exit(1);
    }

    for (i = 0; i < (int)strlen(seq); ++i) {
        p = strchr(alphabet, seq[i]);
        // If character (residue) is not found in alphabet (e.g. X),
        // it's skipped over in the calculation
        // Could be an issue with mantaining sequence-to-score register (Shyam)
        if (p != NULL) {
            // Move characters over in array
            for (j = 1; j < MW; ++j)
                s[MW-j] = s[MW-j-1];
            // Assign new character
            s[0] = p - alphabet;

            // start score calculation when window has been filled up
            if (i >= (MW-1)/2) {
                predict(s,
                        &sm_arr[i-(MW-1)/2],
                        &sb_arr[i-(MW-1)/2],
                        &sr_arr[i-(MW-1)/2]);

#ifdef OLD
                printf("%d\t%f\t%f\t%f\n",
                       i-(MW-1)/2,
                       sr_arr[i-(MW-1)/2],
                       sb_arr[i-(MW-1)/2],
                       sm_arr[i-(MW-1)/2]);
#endif

            }
        }
    }

    // Last bit of sequence
    for (i = (int)strlen(seq)-(MW-1)/2; i < (int)strlen(seq); ++i) {
        // Move characters over in array
        for (j = 1; j < MW; ++j)
            s[MW-j] = s[MW-j-1];
        // Assign new character (placeholder)
        s[0] = NA-1;

        predict(s, &sm_arr[i], &sb_arr[i], &sr_arr[i]);

#ifdef OLD
        printf("%d\t%f\t%f\t%f\n", i, sr_arr[i], sb_arr[i], sm_arr[i]);
#endif

    }

    return;
}


int main(int ARGC, char *ARGV[]) {
    // Read in single sequence line from standard input
    char *buffer = NULL;
    size_t len;
    int seq_len = getline(&buffer, &len, stdin);
    if (seq_len == -1)
        printf("No line read...\n");

    // Calculate scores
    float sm_arr[seq_len], sb_arr[seq_len], sr_arr[seq_len];
    predict_seq(buffer, sm_arr, sb_arr, sr_arr);

  return(0);
}
