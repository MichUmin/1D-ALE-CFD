#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>

#include "global.h"
#include "polynomial_reconstruction.h"

#define IN
#define OUT
#define INOUT

extern double reconstruction_polynomials[NUM_CELLS + 2*NUM_GHOST_CELLS][NUM_VARIABLES][SPACE_ORDER];

#define MATRIX_SIZE (SPACE_ORDER - 1)

double Matrix[MATRIX_SIZE][MATRIX_SIZE];

#define EPSILON 0.000001

double square(double x) {
    return x*x;
}

double power(double x, int n) {
    double result = 1;
    for (int i = 0; i < n; i++) {
        result *= x;
    }
    return result;
}

double max(double a, double b) {
    if (a > b) {
        return a;
    } else {
        return b;
    }
}

double min(double a, double b) {
    if (a < b) {
        return a;
    } else {
        return b;
    }
}

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


double SmoothnessIndicator(double polynomial[SPACE_ORDER], double width) {
    if (SPACE_ORDER == 2) {
        return ((width * width) * (polynomial[1] * polynomial[1]));
    } else if (SPACE_ORDER == 3) {
        double result = (width * width) * (polynomial[1] * polynomial[1]); // first derivative of the 2nd order part
        result += (4.0 + (1.0/3.0)) * (width*width*width*width) * (polynomial[2] * polynomial[2]); // first and second derivatives of the 3rd order part
        return result;
    } else if (SPACE_ORDER == 4) {
        double result = 36.0 * power(width, 6) * square(polynomial[3]);
        result += 3.0 * power(width, 6) * square(polynomial[3]);
        result += 4.0 * power(width, 4) * square(polynomial[2]);
        result += 9.0 / 80.0 * power(width, 6) * square(polynomial[3]);
        result += 1.0 / 3.0 * power(width, 4) * square(polynomial[2]);
        result += square(width) * square(polynomial[1]);
        result += 0.5 * power(width, 4) * polynomial[1] * polynomial[3];
        return result;
    } else {
        printf("TBD higher order WENO");
        abort();
    }
}

/*
double SmoothnessIndicator(double polynomial[SPACE_ORDER], double width) {
    if (SPACE_ORDER == 2) {
        return (width*(polynomial[1] * polynomial[1]));
    } else if (SPACE_ORDER == 3) {
        double result = (width * width) * (polynomial[1] * polynomial[1]); // first derivative of the 2nd order part
        result += (4.0*width + ((1.0/3.0)*width*width*width)) * (polynomial[2] * polynomial[2]); // first and second derivatives of the 3rd order part
        return result;
    } else if (SPACE_ORDER == 4) {
        double result = 36.0 * width * square(polynomial[3]);
        result += 3.0 * power(width, 3) * square(polynomial[3]);
        result += 4.0* width * square(polynomial[2]);
        result += 9.0 / 80.0 * power(width, 5) * square(polynomial[3]);
        result += 1.0 / 3.0 * power(width, 3) * square(polynomial[2]);
        result += width * square(polynomial[1]);
        result += 0.5 * power(width, 3) * polynomial[1] * polynomial[3];
        return result;
    } else {
        printf("TBD higher order WENO");
        abort();
    }
}
*/


void find_reconstruction_polynomial(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS]) {
    double IntermediatePolynomials[SPACE_ORDER][NUM_VARIABLES][SPACE_ORDER];
    double Weights[SPACE_ORDER];
    double distances[SPACE_ORDER];
    double values[NUM_VARIABLES][SPACE_ORDER];
    double SmoothnessIndicators[SPACE_ORDER];
    for (int cell = (SPACE_ORDER-1); cell <= (NUM_CELLS + 2*NUM_GHOST_CELLS - SPACE_ORDER); cell++) {
        // printf("WENO cell %d\n", cell);
        for (int var = 0; var < NUM_VARIABLES; var++) {
            for (int order = 0; order < SPACE_ORDER; order++) {
                reconstruction_polynomials[cell][var][order] = 0.0;
            }
        }
        double main_centre = 0.5*(node_positions[cell] + node_positions[cell+1]);
        double main_width = node_positions[cell+1] - node_positions[cell];
        #if (SPACE_ORDER == 2)
            Weights[0] = 0.5;
            Weights[1] = 0.5;
        #elif (SPACE_ORDER == 3)
            Weights[0] = 0.05;
            Weights[1] = 0.9;
            Weights[2] = 0.05;
        #elif (SPACE_ORDER == 4)
            if (field_values[cell][2] > 0.0) {
                Weights[1] = 0.3;
                Weights[2] = 0.69;
            } else {
                Weights[1] = 0.69;
                Weights[2] = 0.3;
            }
            Weights[0] = 0.005;
            Weights[3] = 0.005;
        #else
            printf("TBD higher order WENO");
            abort();
        #endif

        for (int template_num = 0; template_num < SPACE_ORDER; template_num++) {
            // printf("WENO template %d\n", template_num);
            int template_start = cell - (SPACE_ORDER - 1) + template_num;

            for (int i = 0; i < SPACE_ORDER; i++) {
                int this_cell = template_start + i;
                double this_centre = 0.5*(node_positions[this_cell] + node_positions[this_cell+1]);
                double distance = this_centre - main_centre;
                distances[i] = distance;
                for (int var = 0; var < NUM_VARIABLES; var++) {
                    values[var][i] = field_values[this_cell][var];
                }
            }
            for (int var = 0; var < NUM_VARIABLES; var++) {
                FindPolynomial(values[var], distances, IntermediatePolynomials[template_num][var]);
            }
            double smoothness_here = SmoothnessIndicator(IntermediatePolynomials[template_num][2], main_width);
            // smoothness_here = max(smoothness_here, SmoothnessIndicator(IntermediatePolynomials[template_num][NUM_VARIABLES - 1], main_width));
            Weights[template_num] = Weights[template_num] / power((smoothness_here + EPSILON), 4);
        }
        double weight_sum = 0.0;
        for (int template_num = 0; template_num < SPACE_ORDER; template_num++) {
            weight_sum += Weights[template_num];
        }
        for (int template_num = 0; template_num < SPACE_ORDER; template_num++) {
            Weights[template_num] /= weight_sum;
            for (int var = 0; var < NUM_VARIABLES; var++) {
                for (int ord = 0; ord < SPACE_ORDER; ord++) {
                    reconstruction_polynomials[cell][var][ord] += Weights[template_num]*IntermediatePolynomials[template_num][var][ord];
                }
            }
        }
        for (int var = 0; var < NUM_VARIABLES; var++) {
            double ZeroOrderError = reconstruction_polynomials[cell][var][0] - field_values[cell][var];
            double RelativeZeroOrderError = ZeroOrderError / field_values[cell][var];
            if ((RelativeZeroOrderError > 0.000001) || (RelativeZeroOrderError < -0.000001)) {
                printf("Error in cell %d, wrong 0th order term in WENO polynomial\n", cell);
            }
        }
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


