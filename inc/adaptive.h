
#ifndef ADAPTIVE_H
#define ADAPTIVE_H


#include "kernels.h"
#include "operator.h"
#include "fmm.h"
#include "utility.h"
#include "toptree.h"
#include "remotes.h"

void fmm_solver_adaptive( double halfdkick );

void update_local( int adaptive_level );

void update_p2p_local(int im, int jm, int adaptive_level);

void active_particle (double ai, double af);

void kdk_level (double dkh, double dd, int level);

#endif // ADAPTIVE_H
