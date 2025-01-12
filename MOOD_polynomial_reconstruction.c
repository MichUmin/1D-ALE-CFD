#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>

#include "global.h"
#include "MOOD_polynomial_reconstruction.h"

#define IN
#define OUT
#define INOUT

extern double reconstruction_polynomials[NUM_CELLS + 2*NUM_GHOST_CELLS][NUM_VARIABLES][SPACE_ORDER];
extern bool do_reconstruction[NUM_CELLS + 2*NUM_GHOST_CELLS];

double Matrix[SPACE_ORDER][SPACE_ORDER];

#define EPSILON 0.00001
#define MY_INF 1.0e20

void find_trivial_reconstruction_polynomial(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS]);

void find_limiter_reconstruction_polynomial(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS]);

void find_high_order_reconstruction_polynomial(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS], int order);

void find_reconstruction_polynomial(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS], int order) {
    #ifdef DEBUG
        assert (order >= 1);
    #endif
    if (order == 1) {
        find_trivial_reconstruction_polynomial(field_values, node_positions);
    } else if (order == 2) {
        find_limiter_reconstruction_polynomial(field_values, node_positions);
    } else {
        find_high_order_reconstruction_polynomial(field_values, node_positions, order);
    }
}

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

void factorizeLU(int order) {
    // factorize A=LU using Crout's Algorithm
    for (int p=0; p < order; p++) {
        // divide row of U by the diagonal entry of L
        double diagonal_entry = Matrix[p][p];
        if (diagonal_entry == 0.0) {
            printf("Matrix[%d][%d] = 0\n", p, p);
            abort();
        }
        for (int i = (p+1); i < order; i++) {
            Matrix[p][i] /= diagonal_entry;
        }
        for (int i = (p+1); i < order; i++) {
            for (int j = (p+1) ; j < order; j++) {
                Matrix[i][j] -= Matrix[i][p]*Matrix[p][j];
                // Matrix[j][i] -= Matrix[j][p]*Matrix[p][i];
            }
        }
    }
}

void solveLU(OUT double x[SPACE_ORDER], IN double b[SPACE_ORDER], int order) {
    // solve Lx' = b
    for (int i = 0; i < order; i++) {
        x[i] = b[i];
        for (int j = 0; j < i; j++) {
            x[i] -= Matrix[i][j] * x[j];
        }

        if (i == order-1) {
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
    for (int i = (order - 1); i >= 0; i--) {
        for (int j = i+1; j < order; j++) {
            x[i] -= Matrix[i][j] * x[j];
        }
        x[i] /= 1.0; // U has 1s on the diagonal
        // printf("%lf\n", x[i]);
    }
}

void FindPolynomials(state values[SPACE_ORDER], double points[SPACE_ORDER+1], double centre, int order, double result[NUM_VARIABLES][SPACE_ORDER]) {
    double RHS[NUM_VARIABLES][order];
    for (int i = 0; i < order; i++) {
        double width = points[i+1] - points[i];
        double dLeft = points[i] - centre;
        double dRight = points[i+1] - centre;
        double dLeftPower = dLeft;
        double dRightPower = dRight;
        for (int j = 0; j < order; j++) {
            Matrix[i][j] = (dRightPower - dLeftPower) / ((double)(j+1));
            dLeftPower *= dLeft;
            dRightPower *= dRight;
        }
        for (int var = 0; var < NUM_VARIABLES; var++) {
            RHS[var][i] = values[i][var] * width;
        }
    }
    factorizeLU(order);
    for (int var = 0; var < NUM_VARIABLES; var++) {
        solveLU(result[var], RHS[var], order);
    }
}

void find_high_order_reconstruction_polynomial(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS], int order) {

    for (int cell = (SPACE_ORDER-1); cell <= (NUM_CELLS + 2*NUM_GHOST_CELLS - SPACE_ORDER); cell++) {
        if (do_reconstruction[cell]) {
            for (int var = 0; var < NUM_VARIABLES; var++) {
                for (int ord = 0; ord < SPACE_ORDER; ord++) {
                    reconstruction_polynomials[cell][var][ord] = 0.0;
                }
            }
            double main_centre = 0.5*(node_positions[cell] + node_positions[cell+1]);
            double main_width = node_positions[cell+1] - node_positions[cell];

            int template_start = cell - (order / 2);
            if ((order % 2) == 0) {
                if (field_values[cell][1] > 0.0) {// rho*v positive - fluid moving to the right
                    template_start += 1;
                }    
            }

            FindPolynomials(field_values + template_start, node_positions + template_start, main_centre, order, reconstruction_polynomials[cell]);
        }    
    }
    for (int cell = 0; cell < (SPACE_ORDER-1); cell++) {
        for (int var = 0; var < NUM_VARIABLES; var++) {
            reconstruction_polynomials[cell][var][0] = field_values[cell][var];
            for (int ord = 1; ord < SPACE_ORDER; ord++) {
                reconstruction_polynomials[cell][var][ord] = 0.0;
            }
        }
    }
    for (int cell = (NUM_CELLS + 2*NUM_GHOST_CELLS - SPACE_ORDER + 1); cell < (NUM_CELLS + 2*NUM_GHOST_CELLS); cell++) {
        for (int var = 0; var < NUM_VARIABLES; var++) {
            reconstruction_polynomials[cell][var][0] = field_values[cell][var];
            for (int ord = 1; ord < SPACE_ORDER; ord++) {
                reconstruction_polynomials[cell][var][ord] = 0.0;
            }
        }
    }
}

double minbee(double ratio) {
    if (ratio <= 0.0) {
        return 0.0;
    } else {
        if (ratio <= 1.0) {
            return ratio;
        } else {
            return (2.0 / (1.0 + ratio));
        }
    }
}

void find_limiter_reconstruction_polynomial(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS]) {
    for (int cell = 1; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS - 1); cell++) {
        if (do_reconstruction[cell]) {
            double centre = 0.5 * (node_positions[cell+1] + node_positions[cell]);
            double centreRight = 0.5 * (node_positions[cell+2] + node_positions[cell+1]);
            double centreLeft = 0.5 * (node_positions[cell] + node_positions[cell-1]);
            double dxRight = 0.5 * (node_positions[cell+2] - node_positions[cell]); // centreRight - centre;
            double dxLeft = 0.5 * (node_positions[cell+1] - node_positions[cell-1]); // centre - centreLeft;
            state derivativeLeft, derivativeRight;
            for (int i = 0; i < NUM_VARIABLES; i++) {
                derivativeLeft[i] = (field_values[cell][i] - field_values[cell-1][i])/dxLeft;
                derivativeRight[i] = (field_values[cell+1][i] - field_values[cell][i])/dxRight;
            }
            double limiter;
            if ((derivativeLeft[1] != 0.0) && (derivativeRight[1] != 0.0)) {
                double derivative_ratio = derivativeLeft[1] / derivativeRight[1];
                limiter = minbee(derivative_ratio);
            } else {
                limiter = 0.0;
            }
            #ifdef DEBUG
                assert(((limiter >= 0.0) && (limiter <= 1.0)) || (limiter != limiter));
            #endif
            // limiter = 0.0;
            for (int var = 0; var < NUM_VARIABLES; var++) {
                #ifdef DEBUG
                            assert((dxLeft + dxRight) > 0.00001);
                #endif
                double derivativeCentre = (dxRight*derivativeLeft[var] + dxLeft*derivativeRight[var]) / (dxLeft + dxRight);
                reconstruction_polynomials[cell][var][1] = limiter * derivativeCentre;
                reconstruction_polynomials[cell][var][0] = field_values[cell][var];
                for (int ord = 2; ord < SPACE_ORDER; ord++) {
                    reconstruction_polynomials[cell][var][ord] = 0.0;
                }
            }
        }    
    }
    for (int var = 0; var < NUM_VARIABLES; var++) {
        int last_cell = NUM_CELLS + 2*NUM_GHOST_CELLS - 1;
        reconstruction_polynomials[0][var][0] = field_values[0][var];
        reconstruction_polynomials[last_cell][var][0] = field_values[last_cell][var];
        for (int ord=1; ord < SPACE_ORDER; ord++) {
            reconstruction_polynomials[last_cell][var][ord] = 0.0;
            reconstruction_polynomials[0][var][ord] = 0.0;
        }
    }
}

void find_trivial_reconstruction_polynomial(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS]) {
    for (int cell = 0; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS); cell++) {
        if (do_reconstruction[cell]) {
            for (int var = 0; var < NUM_VARIABLES; var++) {
                reconstruction_polynomials[cell][var][0] = field_values[cell][var];
                for (int ord = 1; ord < SPACE_ORDER; ord++) {
                    reconstruction_polynomials[cell][var][ord] = 0.0;
                }
            }
        }    
    }
}
