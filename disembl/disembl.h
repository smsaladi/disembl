#indef DISEMBL_HEADER
#define DISEMBL_HEADER

void predict_seq(char const *seq, float *sm_arr, float *sb_arr, float *sr_arr);

void predict(int const *s, float *sm, float *sb, float *sr);


#endif DISEMBL_HEADER
