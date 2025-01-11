#define NUM_VARIABLES 3

#define NUM_CELLS 100
#define NUM_NODES (NUM_CELLS + 1)
#define NUM_GHOST_CELLS (SPACE_ORDER + TIME_ORDER)

#define SPACE_ORDER 4
#define TIME_ORDER 2

#define DOMAIN_MIN 0.0
#define DOMAIN_MAX 1.0

typedef double state[NUM_VARIABLES];
