#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdbool.h>

#include "global.h"
#include "polynomial_reconstruction.h"

#define IN
#define OUT
#define INOUT

double Tstart, Tend;

#define GAMMA 1.4
#define CFL 0.25

double reconstruction_polynomials[NUM_CELLS + 2*NUM_GHOST_CELLS][NUM_VARIABLES][SPACE_ORDER];

double node_position[TIME_ORDER+1][NUM_NODES + 2*NUM_GHOST_CELLS];
double node_velocity[NUM_NODES + 2*NUM_GHOST_CELLS];
double widths[NUM_CELLS + 2*NUM_GHOST_CELLS];

state field_in_cell[TIME_ORDER+1][NUM_CELLS + 2*NUM_GHOST_CELLS];

state left_reconstructed[NUM_CELLS + 2*NUM_GHOST_CELLS];
state right_reconstructed[NUM_CELLS + 2*NUM_GHOST_CELLS];



state flux_at_node[NUM_NODES + 2*NUM_GHOST_CELLS];

double cell_width(double nodes[NUM_NODES + 2*NUM_GHOST_CELLS], int index) {
    #ifdef DEBUG
        assert(index < (NUM_CELLS + 2*NUM_GHOST_CELLS));
    #endif
    return (nodes[index+1] - nodes[index]);
}

double pressure(state S) {
    double rho = S[0];
    #ifdef DEBUG
         assert(rho > 0.0);
    #endif
    double v = S[1]/rho;
    double rhoEpsilon = S[NUM_VARIABLES-1] - 0.5*v*S[1];
    #ifdef DEBUG
        assert(rhoEpsilon > 0.0);
    #endif
    double p = (GAMMA - 1)*rhoEpsilon;
    return p;
}

double speed_of_sound(state S) {
    double rho = S[0];
    #ifdef DEBUG
        assert(rho > 0.0);
    #endif
    double p = pressure(S);
    #ifdef DEBUG
        assert(p > 0.0);
    #endif
    return sqrt(GAMMA*(p/rho));
}

void wavespeed_estimate(IN state StateLeft, IN state StateRight, OUT double * sLeft, OUT double * sRight) {
    // very conservative wavespeed estimate
    double cLeft = speed_of_sound(StateLeft);
    double cRight = speed_of_sound(StateRight);
    double rhoLeft = StateLeft[0];
    double rhoRight = StateRight[0];
    double vLeft = StateLeft[1]/rhoLeft;
    double vRight = StateRight[1]/rhoRight;
    double v_max = fmax(fabs(vLeft), fabs(vRight));
    double sPlus = v_max + fmax(cLeft, cRight);
    *sRight = sPlus;
    *sLeft = (-1.0)*sPlus;
}

// void wavespeed_estimate(IN state StateLeft, IN state StateRight, OUT double * sLeft, OUT double * sRight) {
//     double cLeft = speed_of_sound(StateLeft);
//     double cRight = speed_of_sound(StateRight);
//     double rhoLeft = StateLeft[0];
//     double rhoRight = StateRight[0];
//     double vLeft = StateLeft[1]/rhoLeft;
//     double vRight = StateRight[1]/rhoRight;
//     *sRight = fmax(vLeft + cLeft, vRight + cRight);;
//     *sLeft = fmin(vLeft - cLeft, vRight - cRight);
// }

void physical_lagrnge_flux(IN state input, IN double node_v, OUT state result) {
    double rho = input[0];
    double v = input[1]/rho;
    double p = pressure(input);
    result[0] = rho*(v - node_v);
    result[1] = input[1]*(v - node_v) + p;
    result[NUM_VARIABLES-1] = input[NUM_VARIABLES-1]*(v - node_v) + p*v;
}

// void FORCE_flux(IN state Sleft, IN state Sright, IN double node_v, IN double dt, IN double dx, OUT state result) {
//     state Fleft, Fright, intermediate, Richtmyer;
//     physical_lagrnge_flux(Sleft, node_v, Fleft);
//     physical_lagrnge_flux(Sright, node_v, Fright);
//     for (int i = 0; i < NUM_VARIABLES; i++) {
//         intermediate[i] = 0.5 *(Sleft[i] + Sright[i]);
//         intermediate[i] -= 0.5 * (dt/dx) * (Fright[i] - Fleft[i]);
//     }
//     physical_lagrnge_flux(intermediate, node_v, Richtmyer);
//     for (int i = 0; i < NUM_VARIABLES; i++) {
//         result[i] = 0.5*Richtmyer[i];
//         result[i] += 0.25*((dx/dt)*(Sleft[i] - Sright[i]) + Fleft[i] + Fright[i]);
//     }
// }

void HLL_flux(IN state StateLeft, IN state StateRight, IN double node_v, OUT state result) {
    double sLeft, sRight;
    wavespeed_estimate(StateLeft, StateRight, &sLeft, &sRight);
    #ifdef DEBUG
        assert(sRight > sLeft);
    #endif
    if (node_v >= sRight) {
        physical_lagrnge_flux(StateRight, node_v, result);
    } else {
        if (node_v <= sLeft) {
            physical_lagrnge_flux(StateLeft, node_v, result);
        } else {
            state FLeft, FRight;
            physical_lagrnge_flux(StateLeft, node_v, FLeft);
            physical_lagrnge_flux(StateRight, node_v, FRight);
            for (int i = 0; i < NUM_VARIABLES; i++) {
                result[i] = (sRight*FLeft[i] - sLeft*FRight[i] + sLeft*sRight*(StateRight[i] - StateLeft[i])) / (sRight - sLeft);
            }
        }
    }
}


double find_dt(state current_values[NUM_CELLS + 2*NUM_GHOST_CELLS], double current_positions[NUM_NODES + 2*NUM_GHOST_CELLS]) {
    double result = DBL_MAX;
    int guilty_cell = -1;
    for (int i = NUM_GHOST_CELLS; i < (NUM_CELLS + NUM_GHOST_CELLS); i++) {
        double dx = cell_width(current_positions, i);
        double c = speed_of_sound(current_values[i]);
        double local_dt = (CFL) * (dx/c);
        if (local_dt < result) {
            result = local_dt;
            guilty_cell = i;
        }
    }
    // printf("dt due to cell %d\n", guilty_cell);
    return result;
}



// void find_reconstruction_polynomial_limiter(IN state field_values[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_positions[NUM_NODES + 2*NUM_GHOST_CELLS]) {
//     for (int cell = 0; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS); cell++) {
//         for (int var = 0; var < NUM_VARIABLES; var++) {
//             reconstruction_polynomials[cell][var][1] = 0.0;
//             reconstruction_polynomials[cell][var][0] = field_values[cell][var];
//         }
//     }
// }

void reconstruct(state in_the_middle[NUM_CELLS + 2*NUM_GHOST_CELLS], double nodes[NUM_NODES + 2*NUM_GHOST_CELLS]) {
    for (int i = NUM_GHOST_CELLS; i < (NUM_CELLS+NUM_GHOST_CELLS); i++) {
        double length = cell_width(nodes, i);
        double dx = 0.5*length;
        for (int j = 0; j < NUM_VARIABLES; j++) {
            left_reconstructed[i][j] = in_the_middle[i][j];
            right_reconstructed[i][j] = in_the_middle[i][j];
            double dx_power = dx;
            for (int k = 1; k < SPACE_ORDER; k++) {
                right_reconstructed[i][j] += dx_power*reconstruction_polynomials[i][j][k];
                if ((k % 2) == 0) {
                    left_reconstructed[i][j] += dx_power*reconstruction_polynomials[i][j][k];
                } else {
                    left_reconstructed[i][j] -= dx_power*reconstruction_polynomials[i][j][k];
                }
                dx_power *= dx;
            }
        }
    }
    for (int i_left = 0; i_left < NUM_GHOST_CELLS; i_left++) {
        for (int j = 0; j < NUM_VARIABLES; j++) {
            left_reconstructed[i_left][j] = in_the_middle[i_left][j];
            right_reconstructed[i_left][j] = in_the_middle[i_left][j];
            int i_right = i_left + NUM_GHOST_CELLS + NUM_CELLS;
            left_reconstructed[i_right][j] = in_the_middle[i_right][j];
            right_reconstructed[i_right][j] = in_the_middle[i_right][j];
        }
    }
}

// void find_node_volocities(state in_the_middle[]) {
//     for (int node = 1; node < NUM_NODES + 2*NUM_GHOST_CELLS - 1; node++) {
//         #ifdef DEBUG
//             assert(in_the_middle[node-1][0] > 0.0);
//             assert(in_the_middle[node][0] > 0.0);
//         #endif
//         double vLeft  = in_the_middle[node-1][1] / in_the_middle[node-1][0];
//         double vRight = in_the_middle[node][1]   / in_the_middle[node][0];
//         node_velocity[node] = 0.5 * (vLeft + vRight);
//     }
//     int last_node = NUM_NODES + 2*NUM_GHOST_CELLS - 1;
//     int last_cell = last_node - 1;
//     #ifdef DEBUG
//         assert(in_the_middle[last_cell][0] > 0.0);
//         assert(in_the_middle[0][0] > 0.0);
//     #endif
//     node_velocity[last_node] = in_the_middle[last_cell][1]  / in_the_middle[last_cell][0];
//     node_velocity[0] = in_the_middle[0][1]   / in_the_middle[0][0];
// }

#ifdef LAGRANGE
void find_node_volocities() {
     for (int node = 1; node < NUM_NODES + 2*NUM_GHOST_CELLS - 1; node++) {
         double rhoLeft = right_reconstructed[node-1][0];
         double rhoRight = left_reconstructed[node][0];
         #ifdef DEBUG
                 assert(rhoLeft > 0.0);
                 assert(rhoRight > 0.0);
         #endif

         double vLeft  = right_reconstructed[node-1][1] / rhoLeft;
         double vRight = left_reconstructed[node][1]   / rhoRight;
         double sqrtLeft = sqrt(rhoLeft);
         double sqrtRight = sqrt(rhoRight);
         node_velocity[node] = (sqrtLeft*vLeft + sqrtRight*vRight) / (sqrtLeft + sqrtRight);
     }
     int last_node = NUM_NODES + 2*NUM_GHOST_CELLS - 1;
     int last_cell = last_node - 1;
     #ifdef DEBUG
         assert(right_reconstructed[last_cell][0] > 0.0);
         assert(left_reconstructed[0][0] > 0.0);
     #endif
     node_velocity[last_node] = right_reconstructed[last_cell][1]  / right_reconstructed[last_cell][0];
     node_velocity[0] = left_reconstructed[0][1]   / left_reconstructed[0][0];
}
#endif

#ifdef EULER
    void find_node_volocities() {
        for (int node = 0; node < NUM_NODES + 2*NUM_GHOST_CELLS; node++) {
            node_velocity[node] = 0.0;
        }
    }
#endif

void compute_fluxes() {
    for (int node = 1; node < (NUM_NODES + 2*NUM_GHOST_CELLS - 1); node++) {
        HLL_flux(right_reconstructed[node-1], left_reconstructed[node], node_velocity[node], flux_at_node[node]);
        // printf("Flux at node %d done\n", node);
    }
    int last_node = NUM_NODES + 2*NUM_GHOST_CELLS - 1;
    for (int var=0; var < NUM_VARIABLES; var++) {
        flux_at_node[0][var] = flux_at_node[1][var];
        flux_at_node[last_node][var] = flux_at_node[last_node-1][var];
    }
}

void single_update(IN state field_old[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_old[NUM_NODES + 2*NUM_GHOST_CELLS], IN double widths[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double dt, OUT state field_new[NUM_CELLS + 2*NUM_GHOST_CELLS], OUT double node_new[NUM_NODES + 2*NUM_GHOST_CELLS]) {
    for (int cell = 0; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS); cell++) {
        double dx = widths[cell];
        for (int var = 0; var < NUM_VARIABLES; var++) {
            field_new[cell][var] = field_old[cell][var] + (dt/dx)*(flux_at_node[cell][var] - flux_at_node[cell + 1][var]);
        }
    }
    for (int node = 0; node < (NUM_NODES + 2*NUM_GHOST_CELLS); node++) {
        node_new[node] = node_old[node] + (dt * node_velocity[node]);
    }
    #ifdef DEBUG
    for (int node = 1; node < (NUM_NODES + 2*NUM_GHOST_CELLS); node++) {
        assert(node_new[node] > node_new[node-1]);
    }
    #endif
}

// void project_back() {
//     for (int cell = 0; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS); cell++) {
//         if ((node_position[TIME_ORDER][cell] == node_position[0][cell]) && (node_position[TIME_ORDER][cell+1] == node_position[0][cell+1])) {
//             // the cell did not move
//             for (int var = 0; var < NUM_VARIABLES; var++) {
//                 field_in_cell[0][cell][var] = field_in_cell[TIME_ORDER][cell][var];
//             }
//         } else {
//             // cell moved
//             printf("TBD implement remap\n");
//             abort();
//         }
//     }
// }

void integrate(IN double polynomial[NUM_VARIABLES][SPACE_ORDER], IN double left, IN double centre, IN double right, OUT state result) {
    for (int var = 0; var < NUM_VARIABLES; var++) {
        result[var] = 0.0;
    }
    #ifdef DEBUG
        assert(right >= left);
    #endif
    const double dLeft = left-centre;
    const double dRight = right-centre;
    double dLeftPower = 1.0;
    double dRightPower = 1.0;
    // int factorial = 1;
    // for (int ord = 0; ord < SPACE_ORDER; ord++) {
    for (int ord = 0; ord < 1; ord++) {
        // factorial *= (ord+1);
        dLeftPower *= dLeft;
        dRightPower *= dRight;
        for (int var = 0; var < NUM_VARIABLES; var++) {
            result[var] += (dRightPower * polynomial[var][ord] - dLeftPower * polynomial[var][ord]) / ((double)(ord+1));
        }
    }
}

void project_back() {
    bool do_remap = false;
    for (int node = 0; node < (NUM_NODES + 2*NUM_GHOST_CELLS); node++) {
        if (node_position[0][node] != node_position[TIME_ORDER][node]) {
            do_remap = true;
        }
    }
    if (!do_remap) {
        // printf("No remap\n");
        for (int cell = 0; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS); cell++) {
            for (int var = 0; var < NUM_VARIABLES; var++) {
                field_in_cell[0][cell][var] = field_in_cell[TIME_ORDER][cell][var];
            }
        }
    } else {
        // printf("Doing remap\n");
        find_reconstruction_polynomial(field_in_cell[TIME_ORDER], node_position[TIME_ORDER]);
        for (int new_cell = 1; new_cell < (NUM_CELLS + 2*NUM_GHOST_CELLS - 1); new_cell++) {
            double new_left = node_position[0][new_cell];
            double new_right = node_position[0][new_cell + 1];
            double new_centre = 0.5*(new_left + new_right);
            int last_before = 0;
            int first_after = (NUM_NODES + 2*NUM_GHOST_CELLS - 1);
            while (node_position[TIME_ORDER][last_before] <= new_left) {
                last_before++;
            }
            last_before -= 1;
            while (node_position[TIME_ORDER][first_after] >= new_right) {
                first_after--;
            }
            first_after += 1;
            #ifdef DEBUG
                assert(last_before >= 0);
                assert(first_after <= (NUM_NODES + 2*NUM_GHOST_CELLS - 1));
                assert(last_before < first_after);
                assert(node_position[TIME_ORDER][last_before] <= node_position[0][new_cell]);
                assert(node_position[TIME_ORDER][first_after] >= node_position[0][new_cell+1]);
            #endif

            if (first_after == last_before+1) {
                double old_centre = 0.5*(node_position[TIME_ORDER][last_before] + node_position[TIME_ORDER][first_after]);
                state contribution;
                integrate(reconstruction_polynomials[last_before], new_left, old_centre, new_right, contribution);
                for (int var = 0; var < NUM_VARIABLES; var++) {
                    field_in_cell[0][new_cell][var] = contribution[var];
                }
            } else {
                for (int var = 0; var < NUM_VARIABLES; var++) {
                    field_in_cell[0][new_cell][var] = 0.0;
                }
                double old_left, old_centre, old_right;
                state contribution;
                old_centre = 0.5*(node_position[TIME_ORDER][last_before] + node_position[TIME_ORDER][last_before+1]);
                old_right = node_position[TIME_ORDER][last_before+1];
                integrate(reconstruction_polynomials[last_before], new_left, old_centre, old_right, contribution);
                for (int var = 0; var < NUM_VARIABLES; var++) {
                    field_in_cell[0][new_cell][var] += contribution[var];
                }
                for (int old_cell = (last_before+1); old_cell < (first_after-1); old_cell++) {
                    old_left = node_position[TIME_ORDER][old_cell];
                    old_right = node_position[TIME_ORDER][old_cell + 1];
                    old_centre = 0.5*(old_left + old_right);
                    integrate(reconstruction_polynomials[old_cell], old_left, old_centre, old_right, contribution);
                    for (int var = 0; var < NUM_VARIABLES; var++) {
                        field_in_cell[0][new_cell][var] += contribution[var];
                    }
                }
                old_left = node_position[TIME_ORDER][first_after-1];
                old_centre = 0.5*(old_left + node_position[TIME_ORDER][first_after]);
                integrate(reconstruction_polynomials[first_after-1], old_left, old_centre, new_right, contribution);
                for (int var = 0; var < NUM_VARIABLES; var++) {
                    field_in_cell[0][new_cell][var] += contribution[var];
                }
            }
        }
        for (int cell = 1; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS -1); cell++) {
            double length = cell_width(node_position[0], cell);
            for (int var = 0; var < NUM_VARIABLES; var++) {
                field_in_cell[0][cell][var] /= length;
            }
        }

        int last_cell = NUM_CELLS + 2*NUM_GHOST_CELLS - 1;
        for (int var = 0; var < NUM_VARIABLES; var++) {
            field_in_cell[0][0][var] = field_in_cell[TIME_ORDER][0][var];
            field_in_cell[0][last_cell][var] = field_in_cell[TIME_ORDER][last_cell][var];
        }
    }
}


void single_step(IN state field_old[NUM_CELLS + 2*NUM_GHOST_CELLS], IN double node_old[NUM_NODES + 2*NUM_GHOST_CELLS], IN double dt, OUT state field_new[NUM_CELLS + 2*NUM_GHOST_CELLS], OUT double node_new[NUM_NODES + 2*NUM_GHOST_CELLS]) {
    find_reconstruction_polynomial(field_old, node_old);
    reconstruct(field_old, node_old);
    // printf("Reconstruction done\n");
    find_node_volocities();
    compute_fluxes();
    // printf("Fluxes done\n");
    for (int cell = 0; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS); cell++) {
        widths[cell] = node_old[cell+1] - node_old[cell];
    }
    single_update(field_old, node_old, widths, dt, field_new, node_new);
    // printf("Single update done\n");
    project_back();
}

void time_advance(double dt) {
    if (TIME_ORDER == 1) {
        single_step(field_in_cell[0], node_position[0], dt, field_in_cell[1], node_position[1]);
    } else if (TIME_ORDER == 2) {
        // RK2
        find_reconstruction_polynomial(field_in_cell[0], node_position[0]);
        reconstruct(field_in_cell[0], node_position[0]);
        // printf("1st reconstruction done\n");
        find_node_volocities();
        compute_fluxes();
        // printf("1st fluxes done\n");
        for (int cell = 0; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS); cell++) {
            widths[cell] = node_position[0][cell+1] - node_position[0][cell];
        }
        single_update(field_in_cell[0], node_position[0], widths, 0.5*dt, field_in_cell[1], node_position[1]);
        // printf("1st single update done\n");

        find_reconstruction_polynomial(field_in_cell[1], node_position[1]);
        reconstruct(field_in_cell[1], node_position[1]);
        // printf("2nd reconstruction done\n");
        find_node_volocities();
        compute_fluxes();
        // printf("2nd fluxes done\n");
        for (int cell = 0; cell < (NUM_CELLS + 2*NUM_GHOST_CELLS); cell++) {
            widths[cell] = node_position[1][cell+1] - node_position[1][cell];
        }
        single_update(field_in_cell[0], node_position[0], widths, dt, field_in_cell[2], node_position[2]);
        // printf("2nd single update done\n");

        project_back();
    } else {
        printf("TBD: higher order in time\n");
        abort();
    }
}

void initialize() {
    const double initial_domain_size = DOMAIN_MAX - DOMAIN_MIN;
    const double initial_dx = initial_domain_size / NUM_CELLS;
    for (int node = NUM_GHOST_CELLS; node < (NUM_NODES + NUM_GHOST_CELLS); node++) {
        node_position[0][node] = DOMAIN_MIN + (initial_dx * ((double)(node - NUM_GHOST_CELLS)));
    }
    for (int i = 0; i < NUM_GHOST_CELLS; i++) {
        int left_node = NUM_GHOST_CELLS - 1 - i;
        int right_node = NUM_NODES + NUM_GHOST_CELLS + i;
        node_position[0][left_node]  = node_position[0][left_node + 1]  - 0.5*initial_domain_size;
        node_position[0][right_node] = node_position[0][right_node - 1] + 0.5*initial_domain_size;
    }

    double rhoLeft, rhoRight, vLeft, vRight, pLeft, pRight, discontinuity;
    #if TEST == 0
        rhoLeft = 1.0; vLeft = 1.0; pLeft = 1.0; rhoRight = 1.0; vRight = 1.0; pRight = 1.0; discontinuity = 0.5; Tend = 1.0;
    #elif TEST == 1
        rhoLeft = 1.0; vLeft = 0.75; pLeft = 1.0; rhoRight = 0.125; vRight = 0.0; pRight = 0.1; discontinuity = 0.3; Tend = 0.2;
    #elif TEST == 2
        rhoLeft = 1.0; vLeft = -2.0; pLeft = 0.4; rhoRight = 1.0; vRight = 2.0; pRight = 0.4; discontinuity = 0.5; Tend = 0.15;
    #elif TEST == 3
        rhoLeft = 1.0; vLeft = 0.0; pLeft = 1000.0; rhoRight = 1.0; vRight = 0.0; pRight = 0.01; discontinuity = 0.5; Tend = 0.012;
    #elif TEST == 4
        rhoLefty = 5.99924; vLeft = 19.5975; pLeft = 460.894; rhoRight = 5.99242; vRight = -6.19633; pRight = 46.0950; discontinuity = 0.4; Tend = 0.035;
    #elif TEST == 5
         rhoLeft = 1.0; vLeft = -19.59745; pLeft = 1000.0; rhoRight = 1.0; vLeft = -19.59745; pRight = 0.01; discontinuity = 0.8; Tend = 0.012;
    #endif

    double rhovLeft, rhovRight, ELeft, ERight;
    rhovLeft  = rhoLeft  * vLeft;
    rhovRight = rhoRight * vRight;
    ELeft  = pLeft / (GAMMA - 1.0) + 0.5*rhoLeft  * vLeft  * vLeft;
    ERight = pRight/ (GAMMA - 1.0) + 0.5*rhoRight * vRight * vRight;

    for (int cell = 0; cell < (NUM_CELLS + (NUM_GHOST_CELLS*2)); cell++) {
        if (node_position[0][cell+1] <= discontinuity) { // entirely left half
            field_in_cell[0][cell][0] = rhoLeft;
            field_in_cell[0][cell][1] = rhovLeft;
            field_in_cell[0][cell][2] = ELeft;
        } else if (node_position[0][cell] >= discontinuity) { // entirely right half
            field_in_cell[0][cell][0] = rhoRight;
            field_in_cell[0][cell][1] = rhovRight;
            field_in_cell[0][cell][2] = ERight;
        } else { // mixed
            double LeftLength = discontinuity - node_position[0][cell];
            double RightLength = node_position[0][cell + 1] - discontinuity;
            double Length = LeftLength + RightLength;
            field_in_cell[0][cell][0] = (rhoLeft *LeftLength + rhoRight *RightLength) / Length;
            field_in_cell[0][cell][1] = (rhovLeft*LeftLength + rhovRight*RightLength) / Length;
            field_in_cell[0][cell][2] = (ELeft   *LeftLength + ERight   *RightLength) / Length;
        }
    }
}

void print_result(state values[NUM_CELLS + 2*NUM_GHOST_CELLS], double points[NUM_NODES + 2*NUM_GHOST_CELLS]) {
    printf("nodes positions:\n");
    for (int i = 0; i < NUM_NODES + 2*NUM_GHOST_CELLS; i++) {
        printf("%lf ", points[i]);
    }
    printf("\n");
    printf("rho:\n");
    for (int i = 0; i < NUM_CELLS + 2*NUM_GHOST_CELLS; i++) {
        printf("%lf ", values[i][0]);
    }
    printf("\n");
    printf("rho*v:\n");
    for (int i = 0; i < NUM_CELLS + 2*NUM_GHOST_CELLS; i++) {
        printf("%lf ", values[i][1]);
    }
    printf("\n");
    printf("E:\n");
    for (int i = 0; i < NUM_CELLS + 2*NUM_GHOST_CELLS; i++) {
        printf("%lf ", values[i][NUM_VARIABLES-1]);
    }
    printf("\n");
}

int main(void) {
    initialize();
    print_result(field_in_cell[0], node_position[0]);
    double TimeStep;
    double currentTime = Tstart;
    int step = 1;
    while ((currentTime < Tend) && (step < 100000)) {
        TimeStep = find_dt(field_in_cell[0], node_position[0]);
        if ((currentTime + TimeStep) > Tend) {
            TimeStep = Tend - currentTime;
        }
        printf("step %d, dt = %lf\n", step, TimeStep);
        time_advance(TimeStep);
        // printf("Time advancement done in step %d\n", step);
        print_result(field_in_cell[0], node_position[0]);
        currentTime += TimeStep;
        step++;
    }
    print_result(field_in_cell[0], node_position[0]);

    printf("Done in %d steps\nFinal time %lf\n", step-1, currentTime);
    return 0;
}
