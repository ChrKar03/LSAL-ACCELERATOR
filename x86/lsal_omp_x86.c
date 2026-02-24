#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TEST
#define TEST 0
#endif

#define TILE_SIZE 4096

const int match = 2;
const int mismatch = -1;
const int gap_row = -1;
const int gap_col = -1;

int max(int a, int b) { return a > b ? a : b; }
int min(int a, int b) { return a < b ? a : b; }

void lsal_compute_matrices_omp(const char *q, const char *d, size_t *max_idx, int *similarity, char *direction, size_t N, size_t M)
{
    int global_max = 0;
    size_t global_max_idx = 0;

    int num_threads;
    #pragma omp parallel
    {
        #pragma omp single
        num_threads = omp_get_num_threads();
    }

    // Allocate thread-local maxima
    int thread_max[num_threads];
    size_t thread_max_idx[num_threads];
    memset(thread_max, 0, sizeof(thread_max));
    memset(thread_max_idx, 0, sizeof(thread_max_idx));

    // Loop over tiles in wavefront (anti-diagonal) fashion
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        for (int round = 0; round < (N + M - 1 + TILE_SIZE - 1) / TILE_SIZE; round++) {
            #pragma omp for schedule(static)
            for (int tile_row = 0; tile_row < (int)((N + TILE_SIZE - 1) / TILE_SIZE); tile_row++) {
                int tile_col = round - tile_row;
                if (tile_col < 0 || tile_col >= (int)((M + TILE_SIZE - 1) / TILE_SIZE)) continue;

                int row_start = tile_row * TILE_SIZE;
                int col_start = tile_col * TILE_SIZE;

                int row_end = min(row_start + TILE_SIZE, M);
                int col_end = min(col_start + TILE_SIZE, N);

                for (int row = row_start; row < row_end; row++) {
                    for (int col = col_start; col < col_end; col++) {
                        int idx = row * N + col;

                        int score = (d[row] == q[col]) ? match : mismatch;
                        int D = (row > 0 && col > 0) ? similarity[(row - 1) * N + (col - 1)] + score : score;
                        int U = (row > 0) ? similarity[(row - 1) * N + col] + gap_row : gap_row;
                        int L = (col > 0) ? similarity[row * N + (col - 1)] + gap_col : gap_col;

                        int best = 0;
                        char dir = '-';
                        if (D > best) { best = D; dir = 'D'; }
                        if (U > best) { best = U; dir = 'U'; }
                        if (L > best) { best = L; dir = 'L'; }

                        similarity[idx] = best;
                        direction[idx] = dir;

                        if (best > thread_max[tid]) {
                            thread_max[tid] = best;
                            thread_max_idx[tid] = idx;
                        }
                    }
                }
            }
        }
    }

    // Final reduction outside parallel region
    for (int i = 0; i < num_threads; i++) {
        if (thread_max[i] > global_max) {
            global_max = thread_max[i];
            global_max_idx = thread_max_idx[i];
        }
    }

    *max_idx = global_max_idx;
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

    char *q = calloc(qlen + 1, sizeof(char));
    char *d = calloc(dlen + 1, sizeof(char));
    
    init_random_buf(q, qlen);
    init_random_buf(d, dlen);

    
    int *similarity = calloc(qlen * dlen, sizeof(int));
    char *direction = calloc(qlen * dlen, sizeof(char));
    
    size_t max_idx;
    
    #if TEST
    printf("Q: %s\nD: %s\n\n", q, d);
    #endif
    
    #if TEST == 0
    size_t num_iter = 10;

    double omp_time = omp_get_wtime();
    
    for (size_t i = 0; i < num_iter; i++)
    #endif
        lsal_compute_matrices_omp(q, d, &max_idx, similarity, direction, qlen, dlen);

    #if TEST == 0
    omp_time = omp_get_wtime() - omp_time;
    omp_time /= num_iter;
    #endif

    #if TEST
    printf("Q: %s\nD: %s\n\n", q, d);
    #endif
    
    #if TEST
    printf("Max score at: (%lu, %lu)\n", max_idx / qlen, max_idx % qlen);
    
    lsal_print_similarity(q, d, similarity);
    lsal_print_direction(q, d, direction);

    lsal_traceback(q, d, similarity, direction, max_idx);
    #endif

    #if TEST == 0
    printf("Execution Time: %lfs\n", omp_time);
    #endif

    free(similarity);
    free(direction);
    free(q);
    free(d);

    return 0;
}