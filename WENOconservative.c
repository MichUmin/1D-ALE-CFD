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

double Matrix[SPACE_ORDER][SPACE_ORDER];

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
    for (int p=0; p<SPACE_ORDER; p++) {
        // divide row of U by the diagonal entry of L
        double diagonal_entry = Matrix[p][p];
        if (diagonal_entry == 0.0) {
            printf("Matrix[%d][%d] = 0\n", p, p);
            abort();
        }
        for (int i = (p+1); i<SPACE_ORDER; i++) {
            Matrix[p][i] /= diagonal_entry;
        }
        for (int i = (p+1); i < SPACE_ORDER; i++) {
            for (int j = (p+1) ; j < SPACE_ORDER; j++) {
                Matrix[i][j] -= Matrix[i][p]*Matrix[p][j];
                // Matrix[j][i] -= Matrix[j][p]*Matrix[p][i];
            }
        }
    }
}

void solveLU(OUT double x[SPACE_ORDER], IN double b[SPACE_ORDER]) {
    // solve Lx' = b
    for (int i = 0; i < SPACE_ORDER; i++) {
        x[i] = b[i];
        for (int j = 0; j < i; j++) {
            x[i] -= Matrix[i][j] * x[j];
        }

        if (i == SPACE_ORDER-1) {
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
    for (int i = (SPACE_ORDER - 1); i >= 0; i--) {
        for (int j = i+1; j < SPACE_ORDER; j++) {
            x[i] -= Matrix[i][j] * x[j];
        }
        x[i] /= 1.0; // U has 1s on the diagonal
        // printf("%lf\n", x[i]);
    }
}

void FindPolynomials(state values[SPACE_ORDER], double points[SPACE_ORDER+1], double centre, double result[NUM_VARIABLES][SPACE_ORDER]) {
    double RHS[NUM_VARIABLES][SPACE_ORDER];
    for (int i = 0; i < SPACE_ORDER; i++) {
        double width = points[i+1] - points[i];
        double dLeft = points[i] - centre;
        double dRight = points[i+1] - centre;
        double dLeftPower = dLeft;
        double dRightPower = dRight;
        for (int ord = 0; ord < SPACE_ORDER; ord++) {
            Matrix[i][ord] = (dRightPower - dLeftPower) / ((double)(ord+1));
            dLeftPower *= dLeft;
            dRightPower *= dRight;
        }
        for (int var = 0; var < NUM_VARIABLES; var++) {
            RHS[var][i] = values[i][var] * width;
        }
    }
    factorizeLU();
    for (int var = 0; var < NUM_VARIABLES; var++) {
        solveLU(result[var], RHS[var]);
    }
}


double SmoothnessIndicator(double polynomial[SPACE_ORDER], double width) {
    #ifdef DEBUG
        assert(width > 0.0);
    #endif
    if (SPACE_ORDER == 2) {
        return ((width * width) * (polynomial[1] * polynomial[1]));
    } else if (SPACE_ORDER == 3) {
        double result = (width * width) * (polynomial[1] * polynomial[1]); // first derivative of the 2nd order part
        result += (4.0 + (1.0/3.0)) * (width*width*width*width) * (polynomial[2] * polynomial[2]); // first and second derivatives of the 3rd order part
        #ifdef DEBUG
            assert(result >= 0.0);
        #endif
        return result;
    } else if (SPACE_ORDER == 4) {
        double result = 36.0 * power(width, 6) * square(polynomial[3]);
        result += 3.0 * power(width, 6) * square(polynomial[3]);
        result += 4.0 * power(width, 4) * square(polynomial[2]);
        result += 9.0 / 80.0 * power(width, 6) * square(polynomial[3]);
        result += 1.0 / 3.0 * power(width, 4) * square(polynomial[2]);
        result += square(width) * square(polynomial[1]);
        result += 0.5 * power(width, 4) * polynomial[1] * polynomial[3];
        #ifdef DEBUG
            assert(result >= 0.0);
        #endif
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
            Weights[0] = 0.025;
            Weights[1] = 0.95;
            Weights[2] = 0.025;
        #elif (SPACE_ORDER == 4)
            Weights[0] = 0.025;
            Weights[1] = 0.475;
            Weights[2] = 0.475;
            Weights[3] = 0.025;
        #else
            printf("TBD higher order WENO");
            abort();
        #endif

        for (int template_num = 0; template_num < SPACE_ORDER; template_num++) {
            // printf("WENO template %d\n", template_num);
            int template_start = cell - (SPACE_ORDER - 1) + template_num;

            FindPolynomials(field_values + template_start, node_positions + template_start, main_centre, IntermediatePolynomials[template_num]);

            double smoothness_here = 0.0;
            for (int var = 0; var < NUM_VARIABLES; var++) {
                smoothness_here = max(smoothness_here, SmoothnessIndicator(IntermediatePolynomials[template_num][var], main_width));
            }
            // double smoothness_here = SmoothnessIndicator(IntermediatePolynomials[template_num][NUM_VARIABLES-1], main_width);
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




