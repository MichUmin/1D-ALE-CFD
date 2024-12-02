#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "global.h"
#include "polynomial_reconstruction.h"

#define IN
#define OUT
#define INOUT

#define MATRIX_SIZE (SPACE_ORDER - 1)

double Matrix[MATRIX_SIZE][MATRIX_SIZE];

void factorizeLU() {
    // factorize A=LU using Crout's Algorithm
    for (int p=0; p<MATRIX_SIZE; p++) {
        // divide row of U by the diagonal entry of L
        double diagonal_entry = Matrix[p][p];
        if (diagonal_entry == 0.0) {
            printf("Matrix[%d][%d] = 0\n", p, p);
            abort();
        }
        for (int i = (p+1); i<MATRIX_SIZE; i++) {
            Matrix[p][i] /= diagonal_entry;
        }
        for (int i = (p+1); i < MATRIX_SIZE; i++) {
            for (int j = (p+1) ; j < MATRIX_SIZE; j++) {
                Matrix[i][j] -= Matrix[i][p]*Matrix[p][j];
                // Matrix[j][i] -= Matrix[j][p]*Matrix[p][i];
            }
        }
    }
}

void solveLU(OUT double x[MATRIX_SIZE], IN double b[MATRIX_SIZE]) {
    // solve Lx' = b
    for (int i = 0; i < MATRIX_SIZE; i++) {
        x[i] = b[i];
        for (int j = 0; j < i; j++) {
            x[i] -= Matrix[i][j] * x[j];
        }

        if (i == MATRIX_SIZE-1) {
            if (Matrix[i][i] == 0.0) {
                printf("LU[%d][%d] = 0 during solving LU\n", i, i);
                x[i] = 0.0;
            } else {
                x[i] /= Matrix[i][i];
            }
        } else {
            #ifdef DEBUG
                assert(Matrix[i][i] != 0.0);
            #endif
            x[i] /= Matrix[i][i];
        }
        // printf("%lf\n", x[i]);
    }

    // solve Ux = x'
    for (int i = (MATRIX_SIZE - 1); i >= 0; i--) {
        for (int j = i+1; j < MATRIX_SIZE; j++) {
            x[i] -= Matrix[i][j] * x[j];
        }
        x[i] /= 1.0; // U has 1s on the diagonal
        // printf("%lf\n", x[i]);
    }
}


void FindPolynomial(double values[SPACE_ORDER], double points[SPACE_ORDER], double result[SPACE_ORDER]) {
    bool found_centre = false;
    double RHS[SPACE_ORDER-1];
    for (int i = 0; i < (SPACE_ORDER - 1); i++) {
      RHS[i] = 0.0;
    }
    for (int i = 0; i < SPACE_ORDER; i++) {
        if (points[i] == 0.0) {
            result[0] = values[i];
            for (int j = 0; j < (SPACE_ORDER - 1); j++) {
                RHS[j] -= values[i];
            }
            found_centre = true;
        } else {
            int index = i;
            if (found_centre) {
                index--;
            }
            double power = points[i];
            for (int ord = 1; ord < SPACE_ORDER; ord++) {
                Matrix[index][ord-1] = power;
                power *= points[i];
            }
            RHS[index] += values[i];
        }
    }
    factorizeLU();
    solveLU(result+1, RHS);
}

/*
int main(void) {
  double values1[SPACE_ORDER], values2[SPACE_ORDER], points1[SPACE_ORDER], points2[SPACE_ORDER];
  double result[SPACE_ORDER];
  values1[0] = 1.0;
  values1[1] = 2.0;
  values1[2] = 3.0;
  values2[0] = 1.0;
  values2[1] = 2.0;
  values2[2] = 1.0;
  points1[0] = 0.0;
  points1[1] = 2.0;
  points1[2] = 3.0;
  points2[0] = -0.5;
  points2[1] = 0.0;
  points2[2] = 0.75;
  FindPolynomial(values1, points1, result);
  printf("The polynomial is: %f + %f*x + %f*x^2\n", result[0], result[1], result[2]);
    FindPolynomial(values2, points2, result);
    printf("The polynomial is: %f + %f*x + %f*x^2\n", result[0], result[1], result[2]);
  return 0;
}
*/


