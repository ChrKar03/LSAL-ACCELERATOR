#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TEST
#define TEST 0
#endif

const int match = 2;
const int mismatch = -1;
const int gap_row = -1;
const int gap_col = -1;

int max(int a, int b) { return a > b ? a : b; }
int min(int a, int b) { return a < b ? a : b; }

void lsal_compute_matrices_o(char *q, char *d, size_t *max_idx, int *similarity, char *direction, size_t N, size_t M) {
    int max_similarity = 0;
    *max_idx = 0;

    for (size_t row = 0; row < M; row++) {
        for (size_t col = 0; col < N; col++) {
            size_t idx = row * N + col;

            int score = (d[row] == q[col]) ? match : mismatch;

            int D = (row > 0 && col > 0) ? similarity[idx - N - 1] + score : score;
            int U = (row > 0) ? similarity[idx - N] + gap_row : gap_row;
            int L = (col > 0) ? similarity[idx - 1] + gap_col : gap_col;

            int best = 0;
            int dir = '-';

            if (D > best) { best = D; dir = 'D'; }
            if (U > best) { best = U; dir = 'U'; }
            if (L > best) { best = L; dir = 'L'; }

            similarity[idx] = best;
            direction[idx] = dir;
            
            if (best > max_similarity) {
                max_similarity = best;
                *max_idx = idx;
            }
        }
    }
}

void lsal_print_similarity(const char *q, const char *d, int *similarity) {
    size_t N = strlen(q);
    size_t M = strlen(d);

    printf("\n   |");
    for (size_t j = 0; j < N; j++) {
        printf("%3c ", q[j]);
    }
    printf("\n---+");
    for (size_t j = 0; j < N; j++) {
        printf("----");
    }
    printf("\n");

    for (size_t i = 0; i < M; i++) {
        printf("%2c |", d[i]);
        for (size_t j = 0; j < N; j++) {
            printf("%3d ", similarity[i * N + j]);
        }
        printf("\n");
    }
}

void lsal_print_direction(const char *q, const char *d, char *direction) {
    size_t N = strlen(q);
    size_t M = strlen(d);

    printf("\n   |");
    for (size_t j = 0; j < N; j++) {
        printf("%3c ", q[j]);
    }
    printf("\n---+");
    for (size_t j = 0; j < N; j++) {
        printf("----");
    }
    printf("\n");

    for (size_t i = 0; i < M; i++) {
        printf("%2c |", d[i]);
        for (size_t j = 0; j < N; j++) {
            printf("%3c ", direction[i * N + j]);
        }
        printf("\n");
    }
}

void lsal_traceback(const char *q, const char *d, int *similarity, char *direction, size_t max_idx) {
    size_t N = strlen(q);

    char aligned_d[512], aligned_q[512];
    int strpos = 511;
    aligned_d[strpos] = '\0';
    aligned_q[strpos] = '\0';
    strpos--;
    
    int idx = max_idx;
    int row = idx / N;
    int col = idx % N;

    while (row >= 0 && col >= 0 && similarity[idx] > 0) {
        if (direction[idx] == 'D') {
            aligned_q[strpos] = q[col];
            aligned_d[strpos] = d[row];
            row--; col--;
        } else if (direction[idx] == 'U') {
            aligned_q[strpos] = '-';
            aligned_d[strpos] = d[row];
            row--;
        } else if (direction[idx] == 'L') {
            aligned_q[strpos] = q[col];
            aligned_d[strpos] = '-';
            col--;
        }

        idx = row * N + col;
        strpos--;
    }
    
    printf("\nAligned Sequences:\n");
    printf("Q: %s\n", &aligned_q[strpos + 1]);
    printf("D: %s\n", &aligned_d[strpos + 1]);
}

void init_random_buf(char *buf, size_t n) {
    static const char *choices = "ATGC";

    for (size_t i = 0; i < n; i++) {
        buf[i] = choices[rand() % 4];
    }
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <query_length> <database_length>\n", argv[0]);
        return 1;
    }

    int qlen = atoi(argv[1]);
    int dlen = atoi(argv[2]);

    char *q = (char *) calloc(qlen + 1, sizeof(char));
    char *d = (char *) calloc(dlen + 1, sizeof(char));

    init_random_buf(q, qlen);
    init_random_buf(d, dlen);

    int *similarity = (int *) calloc(qlen * dlen, sizeof(int));
    char *direction = (char *) calloc(qlen * dlen, sizeof(char));
    
    size_t max_idx;
    
    #if TEST
        printf("Q: %s\nD: %s\n\n", q, d);
    #endif

    #if TEST == 0
    size_t num_iter = 10;
    clock_t total_time = clock();
    
    for (size_t i = 0; i < num_iter; i++)
    #endif

        lsal_compute_matrices_o(q, d, &max_idx, similarity, direction, qlen, dlen);

    #if TEST == 0
    total_time = clock() - total_time;
    total_time /= num_iter;
    
    double total_time_secs = (double) total_time / (double) CLOCKS_PER_SEC; 
    #endif

    #if TEST
    printf("Max score at: (%lu, %lu)\n", max_idx / qlen, max_idx % qlen);
    
    lsal_print_similarity(q, d, similarity);
    lsal_print_direction(q, d, direction);

    lsal_traceback(q, d, similarity, direction, max_idx);
    #endif

    #if TEST == 0
    printf("Exection Time: %lfs\n", total_time_secs);
    #endif

    free(similarity);
    free(direction);
    free(q);
    free(d);

    return 0;
}
