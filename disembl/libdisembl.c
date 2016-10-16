/***************************************************************/
/*  ____  _     _____ __  __ ____  _        ____      ______   */
/* |  _ \(_)___| ____|  \/  | __ )| |      /___  )   (  __  )  */
/* | | | | / __|  _| | |\/| |  _ \| |         / /    | |  | |  */
/* | |_| | \__ \ |___| |  | | |_) | |___     / /_    | |__| |  */
/* |____/|_|___/_____|_|  |_|____/|_____|  /_____|(_)(______)  */
/* DisEMBL is Copyright (C) 2003                               */
/* Lars Juhl Jensen & Rune Linding - EMBL                      */
/*                                                             */
/* Modifications for v2.0+                                     */
/* Copyright (C) 2016 Shyam Saladi - Caltech                   */
/*                                                             */
/* Licensed under GPL                                          */
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Uncomment to have scores printed to standard output */
// #define OLD

/* Define size of the alphabet */
#define NA 21

/* Define max number of neurons */
#define MW 41
#define MH 30

#include "include/russel.h"
#include "include/bfactor.h"
#include "include/missing.h"
#include "include/libdisembl.h"


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

    if (s[(MW-1)/2] == NA-1) {
      return;
    }

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
predict_seq(char const *seq, float sm_arr[], float sb_arr[], float sr_arr[]) {
    // This ordering specifies how characters (residues)
    // are turned into numbers to be interpreted by the
    // neural network (Shyam)
    char *alphabet = "FIVWMLCHYAGNRTPDEQSK";

    char *p;
    int i, j, s[MW];

    // Fill with default values
    for (i = 0; i < MW; ++i)
        s[i] = NA-1;
    // Sequences less than (MW-1)/2 are effectively padded with `K`
    // This does not seem like the ideal treatment (Shyam)
    // if ((int)strlen(seq) < (MW-1)/2) {
    //     printf("Sequence must be at least %d residues long", (MW-1)/2);
    //     exit(1);
    // }

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
