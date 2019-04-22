#ifndef TOPTREE_H
#define TOPTREE_H

void use_sample_estimate();
void use_time_estimate();

void domaintree_initialize() ;

void domaintree_finalize() ;

void construct_toptree();
void clean_toptree();
void deconstruct_toptree();
void decompose(int *NPART);

void walk_toptree_m2l(int iNode, double L[]) ;
void walk_toptree_m2m(int iNode) ;

void center_toptree(int direct, int iNode, double left[3], double right[3]);
void connect_local_toptree();


#endif 
