#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "global.h"
#include "polynomial_reconstruction.h"

#define IN
#define OUT
#define INOUT

#define EPSILON 0.000001
#define MATRIX_SIZE (SPACE_ORDER-1)

// extern double reconstruction_polynomials[NUM_CELLS + 2*NUM_GHOST_CELLS][NUM_VARIABLES][SPACE_ORDER];
double reconstruction_polynomials[NUM_CELLS + 2*NUM_GHOST_CELLS][NUM_VARIABLES][SPACE_ORDER];

double square(double x) {
    return x*x;
}

double SmoothnessIndicator(double polynomial[SPACE_ORDER], double width) {
    if (SPACE_ORDER == 2) {
        return ((width * width) * (polynomial[1] * polynomial[1]));
    } else if (SPACE_ORDER == 3) {
        double result = (width * width) * (polynomial[1] * polynomial[1]); // first derivative of the 2nd order part
        result += (4.0 + (1.0/3.0)) * (width*width*width*width) * (polynomial[2] * polynomial[2]); // first and second derivatives of the 3rd order part
        return result;
    } else {
        printf("TBD higher order WENO");
        abort();
    }
}

void factorizeLU(INOUT double Matrix[MATRIX_SIZE][MATRIX_SIZE]) {
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

void solveLU(IN double LU[MATRIX_SIZE][MATRIX_SIZE], OUT double x[MATRIX_SIZE], IN double b[MATRIX_SIZE]) {
    // solve Lx' = b
    for (int i = 0; i < MATRIX_SIZE; i++) {
        x[i] = b[i];
        for (int j = 0; j < i; j++) {
            x[i] -= LU[i][j] * x[j];
        }

        if (i == MATRIX_SIZE-1) {
            if (LU[i][i] == 0.0) {
                printf("LU[%d][%d] = 0 during solving LU\n", i, i);
                x[i] = 0.0;
            } else {
                x[i] /= LU[i][i];
            }
        } else {
            #ifdef DEBUG
                assert(LU[i][i] != 0.0);
            #endif
            x[i] /= LU[i][i];
        }
        // printf("%lf\n", x[i]);
    }

    // solve Ux = x'
    for (int i = (MATRIX_SIZE - 1); i >= 0; i--) {
        for (int j = i+1; j < MATRIX_SIZE; j++) {
            x[i] -= LU[i][j] * x[j];
        }
        x[i] /= 1.0; // U has 1s on the diagonal
        // printf("%lf\n", x[i]);
    }
}

void find_reconstruction_polynomial(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS]) {
    double Matrix[SPACE_ORDER-1][SPACE_ORDER-1];
    double RHS[NUM_VARIABLES][SPACE_ORDER-1];
    double IntermediatePolynomials[SPACE_ORDER][NUM_VARIABLES][SPACE_ORDER];
    double Weights[SPACE_ORDER];
    double temp[SPACE_ORDER];
    double SmoothnessIndicators[SPACE_ORDER];
    for (int cell = (SPACE_ORDER-1); cell <= (NUM_CELLS + 2*NUM_GHOST_CELLS - SPACE_ORDER); cell++) {
        // printf("WENO cell %d\n", cell);
        for (int var = 0; var < NUM_VARIABLES; var++) {
            reconstruction_polynomials[cell][var][0] = field_values[cell][var];
            for (int order = 1; order < SPACE_ORDER; order++) {
                reconstruction_polynomials[cell][var][order] = 0.0;
                IntermediatePolynomials[order][var][0] = field_values[cell][var];
            }
        }
        double main_centre = 0.5*(node_positions[cell] + node_positions[cell+1]);
        double main_width = node_positions[cell+1] - node_positions[cell];
        #if (SPACE_ORDER == 2)
            Weights[0] = 0.5;
            Weights[1] = 0.5;
        #elif (SPACE_ORDER == 3)
            Weights[0] = 0.2;
            Weights[1] = 0.6;
            Weights[2] = 0.2;
        #else
            printf("TBD higher order WENO");
            abort();
        #endif

        for (int template_num = 0; template_num < SPACE_ORDER; template_num++) {
            // printf("WENO template %d\n", template_num);
            int template_start = cell - (SPACE_ORDER - 1) + template_num;
            int this_cell = template_start;
            for (int i = 0; i < (SPACE_ORDER-1); i++) {
                if (this_cell == cell) {
                    this_cell++;
                }
                double this_centre = 0.5*(node_positions[this_cell] + node_positions[this_cell+1]);
                double distance = this_centre - main_centre;
                #ifdef DEBUG
                    assert(distance != 0.0);
                    assert(distance != -0.0);
                #endif
                Matrix[i][0] = distance;
                for (int ord = 1; ord < (SPACE_ORDER-1); ord++) {
                    Matrix[i][ord] = Matrix[i][ord-1] * distance;
                }
                for (int var=0; var < NUM_VARIABLES; var++) {
                    RHS[var][i] = field_values[this_cell][var] - field_values[cell][var];
                }
                this_cell++;
            }
            factorizeLU(Matrix);
            // printf("Matrix factorized\n");
            for (int var = 0; var < NUM_VARIABLES; var++) {
                // printf("WENO variable %d\n", var);
                solveLU(Matrix, temp, RHS[var]);
                for (int ord = 1; ord < SPACE_ORDER; ord++) {
                    IntermediatePolynomials[template_num][var][ord] = temp[ord-1];
                }
                // solveLU(Matrix, IntermediatePolynomials[template_num][var] + 1, RHS[var]);
            }
            SmoothnessIndicators[template_num] = SmoothnessIndicator(IntermediatePolynomials[template_num][NUM_VARIABLES-1], main_width);
            Weights[template_num] = Weights[template_num] / square((SmoothnessIndicators[template_num] + EPSILON));
        }
        double weight_sum = 0.0;
        for (int template_num = 0; template_num < SPACE_ORDER; template_num++) {
            weight_sum += Weights[template_num];
        }
        for (int template_num = 0; template_num < SPACE_ORDER; template_num++) {
            Weights[template_num] /= weight_sum;
            for (int var = 0; var < NUM_VARIABLES; var++) {
                for (int ord = 1; ord < SPACE_ORDER; ord++) {
                    reconstruction_polynomials[cell][var][ord] += Weights[template_num]*IntermediatePolynomials[template_num][var][ord];
                }
            }
        }
        // printf("\n");
    }
    for (int cell = 0; cell < (SPACE_ORDER -1); cell++) {
        for (int var = 0; var < NUM_VARIABLES; var++) {
            reconstruction_polynomials[cell][var][0] = field_values[cell][var];
            for (int order = 1; order < SPACE_ORDER; order++) {
                reconstruction_polynomials[cell][var][order] = 0.0;
            }
        }
    }
    for (int cell = (NUM_CELLS + 2*NUM_GHOST_CELLS - SPACE_ORDER + 1); cell < (NUM_CELLS + 2*NUM_GHOST_CELLS); cell++) {
        for (int var = 0; var < NUM_VARIABLES; var++) {
            reconstruction_polynomials[cell][var][0] = field_values[cell][var];
            for (int order = 1; order < SPACE_ORDER; order++) {
                reconstruction_polynomials[cell][var][order] = 0.0;
            }
        }
    }
}

/*
int main() {
    double A[MATRIX_SIZE][MATRIX_SIZE];
    double b[MATRIX_SIZE];
    double x[MATRIX_SIZE];

    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] = i*MATRIX_SIZE + j + 1.0;
        }
        b[i] = i;
        // x[i] = 0.0;
    }
    A[2][2] -= 1.0;
    // A[0][0] = 2.0;
    // A[0][1] = 4.0;
    // A[1][0] = -1.0;
    // A[1][1] = 6.0;
    //
    // b[0] = 1.0;
    // b[1] = 1.0;
    //
    printf("A\n");
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
    factorizeLU(A);
    solveLU(A, x, b);

    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
    printf("x\n");
    for (int i = 0; i < MATRIX_SIZE; i++) {
        printf("%lf\n", x[i]);
    }
}
*/

