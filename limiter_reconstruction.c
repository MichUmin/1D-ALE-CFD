#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "global.h"
#include "polynomial_reconstruction.h"

#define IN
#define OUT
#define INOUT

extern double reconstruction_polynomials[NUM_CELLS + 2*NUM_GHOST_CELLS][NUM_VARIABLES][SPACE_ORDER];

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

void find_reconstruction_polynomial(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS]) {
    for (int cell = 1; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS - 1); cell++) {
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
        if ((derivativeLeft[NUM_VARIABLES-1] != 0.0) && (derivativeRight[NUM_VARIABLES-1] != 0.0)) {
            double derivative_ratio = derivativeLeft[NUM_VARIABLES-1] / derivativeRight[NUM_VARIABLES-1];
            limiter = minbee(derivative_ratio);
        } else {
            limiter = 0.0;
        }
        #ifdef DEBUG
            assert((limiter >= 0.0) && (limiter <= 1.0));
        #endif
        // limiter = 0.0;
        for (int var = 0; var < NUM_VARIABLES; var++) {
            #ifdef DEBUG
                        assert((dxLeft + dxRight) > 0.00001);
            #endif
            double derivativeCentre = (dxRight*derivativeLeft[var] + dxLeft*derivativeRight[var]) / (dxLeft + dxRight);
            reconstruction_polynomials[cell][var][1] = limiter * derivativeCentre;
            reconstruction_polynomials[cell][var][0] = field_values[cell][var];
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
