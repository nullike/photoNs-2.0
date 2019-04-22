/*
 * photoNs-2
 * 
 * Cosmological Simulation 
 * on flat LambdaCDM model
 *
 *      2017 - 8 - 21
 *	qwang@nao.cas.cn
 */	 
#ifndef PHOTONS_H
#define PHOTONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utility.h"
#include <mpi.h>
#include "typesdef.h"

double open_angle;
double fixed_stepping;
int adaptive_level_maximum;
int first_node, last_node;
int first_leaf, last_leaf;
int MAXLEAF;

int PROC_RANK, PROC_SIZE;
double BOXSIZE;
double InitialTime;
double BoxMinimum;
double BoxCenter;
double BoxMaximum;

int NPART;
int NLEAF;
int NNODE;
int NPART_MEAN;

long NPART_TOTAL;
long data_length;
double MASSPART;
double GravConst;
double SoftSquare;
double SoftenScale;
double  *data;

double OmegaM0;
double OmegaX0;
double Hubble0;
double Redshift_Time;

////////////////// PM ////////////////
typedef struct {
	int x, y, z;
	double v;
} MKey;

int NSIDE;
int NPARTSIDE;
int local_xstart[3];
int local_xend[3];
int local_xsize[3];
int direct_local_start;

//MPI_Comm PM_COMM_WORLD;
////////////////// PM ////////////////


int NumThread;
int vproc[2];
int pside[2];

int *MSIZE1;
int *MSIZE0;

int SnapNumber;
int SnapFormat;

int p_side;

char OutputPath[128];
char OutputName[128];
char InputPath[128];

char CodeProj[120];

//////////////////// TIMER ////////////////////


double dtime();
double DTIME_THIS_DOMAIN;
double DTIME_FRACTION;
double DTIME_PARTMESH;

double dtime_p2p;
double dtime_m2l;
double dtime_p2p_mirror;
double dtime_p2p_remote;
double dtime_p2p_adlc;
double dtime_p2p_adex;

double dtime_fmm;
double dtime_fmm_remote;
double dtime_fmm_adlc;
double dtime_fmm_adex;

double dtime_pm;
double dtime_prep;
double dtime_task;
double dtime_ext;
double dtime_prep_ext;

///////////////////////////////////////////////

#define MAXSNAP 128
int number_snaptime;
double snaptime[MAXSNAP];


//////////////////////////////

int *ExtDomain;  // free
int nExtDomain;
double width_this_domain;
double cutoffRadius;
double splitRadius;
double *gravfunc;
int NumGravFunc;
double invGravFunc;

double SPLIT_THIS_DOMAIN;

INT64 accept_count;
INT64 open_count;
INT64 active_count;
INT64 p2p_count;
INT64 p2p_count_remote;
INT64 p2p_count_inleaf;
INT64 pack2pack_count;
INT64 pack2pack_inleaf_count;
INT64 walk_m2l_count;
INT64 m2l_count;

INT64 idxM2L;
INT64 idxP2P;

INT64 idxM2Lext;
INT64 idxP2Pext;

int *taskM2Ls;
int *taskM2Lt;

int *taskP2Ps;
int *taskP2Pt;




int limiter_domain;

int count_crushing_step ;
int loop_step;




//////////////////////////////////////////////

MPI_Datatype strMKey;
MPI_Datatype strReNode;
MPI_Datatype strReBody;
MPI_Datatype strBody;
//MPI_Datatype strSample;

///////////////////////////////////////////////////

typedef struct {
	int npart;
	int  son[NSON];
	double width[3];
	double center[3];
	double M[NMULTI];
} RemoteNode;

typedef struct {
	double pos[3];
	//double mass;
	double replenish;
} RemoteBody;
typedef struct{
	   RemoteBody *slave_remote;
	   Body *slave_part;
	   double rs;
	   double coeff;
           double SoftenScale;	 
	  double mass;
	   int *s_i_part;
	   int *s_j_part;
	   int idxp2p;
}slave_data;

int *taskM2Lexts;
int *taskM2Lextt;

int *taskP2Pexts;
int *taskP2Pextt;

RemoteNode *sendtree;
RemoteBody *sendbody;

RemoteNode *recvtree;
RemoteBody *recvbody;

RemoteNode *exstree;
RemoteBody *exsbody;

RemoteNode *exrtree;
RemoteBody *exrbody;

int last_sendtree;
int last_sendbody;

void *bk_sendtree, *bk_sendbody;
void *bk_recvtree, *bk_recvbody;

int *sendcnt_node ;
int *sendisp_node ;
int *sendcnt_body ;
int *sendisp_body ;
int *recvcnt_node ;
int *recdisp_node ;

int *recvcnt_body ;
int *recdisp_body ;

RemoteNode *sendtreem;
RemoteBody *sendbodym;

RemoteNode *recvtreem;
RemoteBody *recvbodym;

int last_sendtreem;
int last_sendbodym;

void *bk_sendtreem, *bk_sendbodym;
void *bk_recvtreem, *bk_recvbodym;

int *sendcnt_nodem ;
int *sendisp_nodem ;
int *sendcnt_bodym ;
int *sendisp_bodym ;
int *recvcnt_nodem ;
int *recdisp_nodem ;

int *recvcnt_bodym ;
int *recdisp_bodym ;

RemoteNode *sendtree_remote;
RemoteBody *sendbody_remote;

RemoteNode *recvtree_remote;
RemoteBody *recvbody_remote;

int last_sendtree_remote;
int last_sendbody_remote;

void *bk_sendtree_remote, *bk_sendbody_remote;
void *bk_recvtree_remote, *bk_recvbody_remote;

int *sendcnt_node_remote ;
int *sendisp_node_remote ;
int *sendcnt_body_remote ;
int *sendisp_body_remote ;
int *recvcnt_node_remote ;
int *recdisp_node_remote ;

int *recvcnt_body_remote ;
int *recdisp_body_remote ;


/////////////////////////////////////////////
//    (toptree)              (domain)
//  0  ... ...  first_domain ... ...  last_domain

typedef struct {
	int  son[NSON];
	double split;
	double M[NMULTI];
	double center[3];
	double width[3];
} TopNode;

TopNode *toptree;

int first_domain;
int last_domain;
int this_domain;
int mostleft;

typedef struct {
	int  son[NSON];
	double time_node;
	double time_left;
	double time_right;
	double split;
	double split_previous;
} DomainNode;

DomainNode *domtree;

#ifdef PMTHREAD
MPI_Comm PM_COMM_WORLD;
#endif
///////////////////////////////////////////////////////

double SamplingRate;
long int seed;
float ran3(long *idum);

///////////////////////////////


int DEBUG_RUNTIME;

int acceptance(double wi[3], double wj[3], double dist[3]) ;
void m2m(double dx, double dy, double dz, double M[], double tM[]);
void p2m(int iPart, int npart, double center[3], double M[]);
void m2l(double dx, double dy, double dz, double M[], double toL[]) ;
void walk_m2m(int iNode);
void walk_toptree_m2m(int iNode);
void initialize_nbody(char argv[1]);
/////////////////////////////////////////////////////


#endif /// PHOTONS_H ///
