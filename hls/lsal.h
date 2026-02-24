#ifndef LSAL_H
#define LSAL_H

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ap_int.h>

extern "C" {
void lsal_compute_matrices_aug(char *q, char *d, int *max_idx, char *direction);
}
#endif
