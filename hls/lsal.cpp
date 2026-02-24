#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lsal.h"

#define N 32
#define M 65536
#define M_p (M + 2 * (N - 1))

#define match 2
#define mismatch -1
#define gap_row -1
#define gap_col -1

#define DIR_D 1
#define DIR_U 2
#define DIR_L 3
#define DIR_NONE 0

void lsal_compute_matrices_aug(char *q,
							   char *d,
							   int *max_idx,
							   char *direction)
{
#pragma HLS TOP name=lsal_compute_matrices_aug
#pragma HLS INTERFACE s_axilite port=max_idx
#pragma HLS INTERFACE m_axi port=d bundle=hp0 offset=slave
#pragma HLS INTERFACE m_axi port=q bundle=hp1 offset=slave
#pragma HLS BIND_STORAGE variable=q type=ram_t2p impl=autosrl
#pragma HLS BIND_STORAGE variable=d type=ram_t2p impl=autosrl
#pragma HLS BIND_STORAGE variable=direction type=ram_t2p impl=autosrl

    int8_t max_similarity = 0;
    int max_idx_tmp = 0;

    char q_buf[N], d_buf[N];
#pragma HLS ARRAY_PARTITION variable=q_buf dim=1 complete
#pragma HLS ARRAY_PARTITION variable=d_buf dim=1 complete

    int8_t buf_curr[N], buf_prev_1[N + 1], buf_prev_2[N + 1];
#pragma HLS ARRAY_PARTITION variable=buf_curr dim=1 complete
#pragma HLS ARRAY_PARTITION variable=buf_prev_1 dim=1 complete
#pragma HLS ARRAY_PARTITION variable=buf_prev_2 dim=1 complete

    char dir_buf[N];
#pragma HLS ARRAY_PARTITION variable=dir_buf dim=1 complete

    int max_row_buf[N];
    int8_t max_value_buf[N];
#pragma HLS ARRAY_PARTITION variable=max_row_buf dim=1 complete
#pragma HLS ARRAY_PARTITION variable=max_value_buf dim=1 complete

    q: memcpy(q_buf, q, N * sizeof(char));
    d: memcpy(d_buf, d, N * sizeof(char));

    sim_prev_1: memset(buf_prev_1, 0, (N + 1) * sizeof(int8_t));
    sim_prev_2: memset(buf_prev_2, 0, (N + 1) * sizeof(int8_t));

    max_idx: memset(max_row_buf, 0, N * sizeof(int));
    max_val: memset(max_value_buf, 0, N * sizeof(int8_t));

    Round: for (int row = 0; row < (M + N - 1); row++) {
#pragma HLS PIPELINE
        Off: for (int col = 0; col < N; col++) {
        	char d_char = d_buf[row + N - 1 - col];
        	char q_char = q_buf[col];
        	int8_t buf_D = buf_prev_2[col];
        	int8_t buf_U = buf_prev_1[col + 1];
        	int8_t buf_L = buf_prev_1[col];
        	int8_t max_value = max_value_buf[col];

        	int8_t score = (d_char == q_char) ? match : mismatch;

            int8_t D = buf_D + score;
            int8_t U = buf_U + gap_row;
            int8_t L = buf_L + gap_col;

            int8_t best1, best2, best;
            char dir1, dir2, dir;

            if (D > 0) {
            	best1 = D;
            	dir1 = DIR_D;
            } else {
            	best1 = 0;
            	dir1 = DIR_NONE;
            }

            if (U > L) {
            	best2 = U;
            	dir2 = DIR_U;
            } else {
            	best2 = L;
            	dir2 = DIR_L;
            }

            if (best1 > best2) {
            	best = best1;
            	dir = dir1;
            } else {
            	best = best2;
            	dir = dir2;
            }

            buf_curr[col] = best;
            dir_buf[col] = dir;

            if (best > max_value) {
               	max_value_buf[col] = best;
               	max_row_buf[col] = row;
            }
        }

        sim_2_prev: memcpy(buf_prev_2 + 1, buf_prev_1 + 1, N * sizeof(int8_t));
        sim_1_prev: memcpy(buf_prev_1 + 1, buf_curr, N * sizeof(int8_t));

        dir: memcpy(direction + (N * row), dir_buf, N * sizeof(char));

        memcpy(d_buf, d_buf + 1, (N - 1) * sizeof(char));
        d_buf[N - 1] = d[row + N];
    }

    max_final: for (int col = 0; col < N; col++) {
#pragma HLS UNROLL
    	if (max_similarity <= max_value_buf[col]) {
    		max_similarity = max_value_buf[col];
    		max_idx_tmp = max_row_buf[col] * N + col;
    	}
    }

    *max_idx = max_idx_tmp;
}
