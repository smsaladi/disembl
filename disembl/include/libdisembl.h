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

#indef LIBDISEMBL_H
#define LIBDISEMBL_H

void predict_seq(char const *seq, float *sm_arr, float *sb_arr, float *sr_arr);

#endif
