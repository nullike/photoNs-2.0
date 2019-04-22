#ifndef INIT_H
#define INIT_H

void initialize_nbody(char parameterfile[]);
void finalize_nbody();
void simple_snapshot(char fbase[]);
void find_boundary();

double drift_loga(double loga_i, double loga_f) ;
double kick_loga(double loga_i, double loga_f) ;

double t_flat_lcdm_a(double a);


#endif

