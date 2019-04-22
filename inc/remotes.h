#ifndef REMOTE_H
#define REMOTE_H


void p2p_kernel_remote(int inode, int jnode);
void walk_m2l_remote(int im, int jm);
void prepare_sendtree_remote(int isend, int ilocal, int tNode, int D);
void prepare_sendtree(int isend, int ilocal, int tNode, int D);

void p2p_kernel_mirror(int inode, int jnode);
void walk_m2l_mirror(int im, int jm);
void prepare_sendtree_mirror(int isend, int ilocal, int tNode, int D, double displace[3]);

void prepare_sendtree_mirror2(int isend, int ilocal, int tNode, int D, double displace[3]);

void fmm_remote(int idx, double shift[3]);
#endif
