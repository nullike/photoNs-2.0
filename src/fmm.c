/*
 *     photoNs-2
 *
 *     2017 - 8 - 21
 *	qwang@nao.cas.cn
 *
 *     modified by Q. Wang on 19 Apr 2019
 */	 

#include <omp.h>
#include <math.h>
#include "photoNs.h"
#include "fmm.h"
#include <pthread.h>


static int *task_s[2];
static int *task_t[2];
static double *task_accx[2];
static double *task_accy[2];
static double *task_accz[2];
static int ntask[2];

static int stask;
static int ttask;
static int idxtask;
static int LEN_TASK;


void bksort_inplace(int D, Body *p, int length, int npart[2], double *split)
{
	npart[0] = npart[1] = 0;
	if (length <2 ) {
		npart[1] = length;
		return;
	}
	Body temp;
	if (length == 2 )
	{
		npart[0] = npart[1] = 1;
		*split = 0.5*(p[0].pos[D] + p[1].pos[D]);
		if (p[0].pos[D] > p[1].pos[D]) {
			temp = p[0];
			p[0] = p[1];
			p[1] = temp;
		}
		return;
	}
	int n;

	double mean = 0.0;

	for (n=0; n<length; n++) {
		mean += p[n].pos[D];
	}
	mean /= (double)length;

	int but = length - 1;

	n=0;
	while ( n < but ) {
		if ( p[n].pos[D] > mean ) {
			temp = p[n];

			while ( p[but].pos[D] > mean && but > n )
				but--;

			p[n] = p[but];
			p[but] = temp;
		}
		n++;
	}

	npart[0] = but;
	npart[1] = length - npart[0];

	*split = mean;
}


void build_kdtree(int direct, int iPart, int length, int iNode)
{
	if (length == 0)
		return;


	btree[iNode].npart = length;
	Body *p = part + iPart;

	int npart[2];
	double split;

	bksort_inplace(direct, p, length, npart, &split);

	btree[iNode].split = split;

	int n, ip;
	ip = iPart;
	int new_direct = (direct + 1)%3;

	for (n=0; n<NSON; n++) {

		if (npart[n]  <= MAXLEAF) {
			leaf[last_leaf].npart = npart[n];
			leaf[last_leaf].ipart = ip;

			btree[iNode].son[n] = last_leaf;

			last_leaf++;
		}
		else {
			last_node++;
			btree[iNode].son[n] = last_node;
			build_kdtree(new_direct, ip, npart[n],  last_node);
		}

		ip += npart[n];
	}

}


void center_kdtree(int direct, int iNode, double left[3], double right[3])
{
	int n, idx;
	btree[iNode].width[0] = right[0]-left[0];
	btree[iNode].width[1] = right[1]-left[1];
	btree[iNode].width[2] = right[2]-left[2];
	btree[iNode].center[0] = 0.5*(right[0]+left[0]);
	btree[iNode].center[1] = 0.5*(right[1]+left[1]);
	btree[iNode].center[2] = 0.5*(right[2]+left[2]);

	int newd = (direct + 1)%3;
	double tmp;
	for (n=0; n<NSON; n++) {
		idx = btree[iNode].son[n];

		if (idx < last_leaf ) {

			leaf[idx].width[0] = btree[iNode].width[0];
			leaf[idx].width[1] = btree[iNode].width[1];
			leaf[idx].width[2] = btree[iNode].width[2];

			leaf[idx].center[0] = btree[iNode].center[0];
			leaf[idx].center[1] = btree[iNode].center[1];
			leaf[idx].center[2] = btree[iNode].center[2];

			if (0==n){
				leaf[idx].width[direct] = btree[iNode].split-left[direct];
				leaf[idx].center[direct] = 0.5*(left[direct]+btree[iNode].split);
			}

			if (1==n) {
				leaf[idx].width[direct] = right[direct] - btree[iNode].split;
				leaf[idx].center[direct] = 0.5*(right[direct]+btree[iNode].split);
			}


		} else {

			if (0==n) {
				tmp = right[direct];
				right[direct] = btree[iNode].split;
				center_kdtree(newd, idx, left, right);
				right[direct] =tmp;
			}
			if (1==n) {
				tmp = left[direct];
				left[direct] = btree[iNode].split;
				center_kdtree(newd, idx, left, right);
				left[direct] =tmp;

			}
		}
	}

}


void build_localtree() {

	dtime_p2p = 0.0;
	dtime_p2p_remote = 0.0;
	dtime_p2p_mirror = 0.0;

	dtime_p2p_adlc = 0.0;
	dtime_p2p_adex = 0.0;

	double bdleft[3];
	double bdright[3];
	int d;
	int direct = direct_local_start;

	for (d=0; d<3; d++) {
		bdleft[d] = toptree[this_domain].center[d] - 0.5*toptree[this_domain].width[d];
		bdright[d] = toptree[this_domain].center[d] + 0.5*toptree[this_domain].width[d];
	}

	/*
	// 0...NPART, fist_leaf...last_leaf, first_node...last_node
	*/

	NNODE = (int)(2.0*((double)NPART)/((double)MAXLEAF));
	NLEAF = (int)(2.0*((double)NPART)/((double)MAXLEAF));

	if (NNODE > NPART)
		NNODE = NPART +1 ;
	if (NLEAF > NPART)
		NLEAF = NPART +1;

	first_leaf = last_leaf = NPART;
	first_node = last_node = NPART+NLEAF;

	// 0...NPART, fist_leaf...last_leaf, first_node...last_node

	leaf  = (Pack*)pmalloc(sizeof(Pack)*NLEAF, 30);
	btree = (Node*)pmalloc(sizeof(Node)*NNODE, 31);

	int n, m;
	for (n=0; n<NNODE; n++) {
		btree[n].npart = 0;
		btree[n].son[0] = btree[n].son[1] = -1;
		btree[n].split = 0.0;
		btree[n].width[0] = 0.0;
		btree[n].width[1] = 0.0;
		btree[n].width[2] = 0.0;

		btree[n].center[0] = 0.0;
		btree[n].center[1] = 0.0;
		btree[n].center[2] = 0.0;

		for (m=0; m<NMULTI; m++) {
			btree[n].M[m] = 0.0;
			btree[n].L[m] = 0.0;
		}
	}


	for (n=0; n<NLEAF; n++) {
		leaf[n].npart = 0;
		leaf[n].ipart = 0;

		leaf[n].width[0] =0.0;
		leaf[n].width[1] =0.0;
		leaf[n].width[2] =0.0;

		leaf[n].center[0] = 0.0;
		leaf[n].center[1] = 0.0;
		leaf[n].center[2] = 0.0;

		for (m=0; m<NMULTI; m++) {
			leaf[n].M[m] = 0.0;
			leaf[n].L[m] = 0.0;
		}

	}

	btree -= first_node;
	leaf  -= first_leaf;

	build_kdtree(direct,0, NPART,  first_node);

	center_kdtree(direct, first_node, bdleft, bdright);
}

// 0==open, 1==accept, -1==abort
inline int acceptance(double wi[3], double wj[3], double dist[3]) 
{
	double min[3];
	double w[3];
    	double dd2, dm2, comp = 1.0;

	w[0] = ( wi[0] + wj[0] ) * 0.5;	
	w[1] = ( wi[1] + wj[1] ) * 0.5;	
	w[2] = ( wi[2] + wj[2] ) * 0.5;	

	min[0] = dist[0];
	min[1] = dist[1];
	min[2] = dist[2];

	dd2 = dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];

	if (min[0] < 0.0)
		min[0] = - min[0];
	if (min[1] < 0.0)
		min[1] = - min[1];
	if (min[2] < 0.0)
		min[2] = - min[2];

	min[0] -= w[0];
	min[1] -= w[1];
	min[2] -= w[2];

	if (min[0] <= 0.0)
		min[0] = 0.0;
	if (min[1] <= 0.0)
		min[1] = 0.0;
	if (min[2] <= 0.0)
		min[2] = 0.0;

	// neighbour
	if (min[0] + min[1] + min[2] < 0.0001)
		return 0;

        dm2 =  min[0]*min[0] + min[1]*min[1] + min[2]*min[2] ;

#ifdef LONGSHORT
	double c2 = cutoffRadius * cutoffRadius;
        if ( dm2 >= c2) 
                return -1;
        
	if ( dd2 > 1.0*c2 )
		return 0;

#endif
	double wmax = w[0];
	if (w[1]> wmax)wmax = w[1];
	if (w[2]> wmax)wmax = w[2];
	wmax *= 2;
       	 if (  wmax*wmax < open_angle * open_angle * dd2 )  {
		return 1;
	}
        else
                return 0;

}


void fmm_construct( ) {

	int rank, d;
	int n, mi, mj, mk;
	double t0, t1, t2, t3, t4, t5, t6, t7, t8;

	t0 = dtime();

	construct_toptree(PROC_SIZE);

	for (n=0; n<=last_domain; n++) {
		toptree[n].split = domtree[n].split;
	}

#ifdef PERIODIC_CONDITION

	double bdl[3] = {     0.0,     0.0,     0.0 };
	double bdr[3] = { BOXSIZE, BOXSIZE, BOXSIZE };

#else

	double bdl[3] = { BoxMinimum, BoxMinimum, BoxMinimum };
	double bdr[3] = { BoxMaximum, BoxMaximum, BoxMaximum };

#endif

	int direct = 0;


	center_toptree(direct, 0, bdl, bdr);

}



////////////////////////////////////////////////////////////////

///// need stort im, jm and for m2l or p2p list.

int p1st = 1;
int alt;
int argt[2];

int *ts;
int *tt;

pthread_t tid;
void *status;

void *task_compute_p2p(void *arg);


void turn2compute_p2p(){
	if (p1st == 1) {
		p1st = 0;
	} 
	else {
		pthread_join(tid, &status);
	}
	ntask[alt] = idxtask;

	argt[0] = alt;
	argt[1] = ntask[alt];


	pthread_create(&tid, NULL, task_compute_p2p, (void*)argt);

	alt = (alt+1)%2;

	ts = task_s[alt];
	tt = task_t[alt];

	idxP2P +=idxtask;
	idxtask = 0;

}

void walk_task_p2p(int im, int jm)
{
	if ( -1 == im  || -1 == jm)
		return;

	if ( im < first_leaf || jm < first_leaf ) {
		printf("error\n");
		exit(0);
	}

	if ( im == jm ) {
		if ( im <  last_leaf ) {	
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;

			if ( idxtask == LEN_TASK ) {
				turn2compute_p2p();
			}
		}
		if ( im >= first_node ) {
			walk_task_p2p(btree[im].son[0], btree[jm].son[0]);
			walk_task_p2p(btree[im].son[0], btree[jm].son[1]);
			walk_task_p2p(btree[im].son[1], btree[jm].son[0]);
			walk_task_p2p(btree[im].son[1], btree[jm].son[1]);
		}
		return;
	}
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < last_leaf && jm < last_leaf ) {
		*(ts + idxtask) = jm;
		*(tt + idxtask) = im;
		idxtask ++;

		if ( idxtask == LEN_TASK ) {

			turn2compute_p2p();

		}
		return;
	}

	if ( im < first_node && jm >= first_node ) 
	{
		dx = leaf[im].center[0] - btree[jm].center[0];
		dy = leaf[im].center[1] - btree[jm].center[1];
		dz = leaf[im].center[2] - btree[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {
		}
		else if (0 == flag) {
			walk_task_p2p(im, btree[jm].son[0]);
			walk_task_p2p(im, btree[jm].son[1]);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}	

	if ( im >= first_node && jm < first_node ) {
		dx = btree[im].center[0] - leaf[jm].center[0];
		dy = btree[im].center[1] - leaf[jm].center[1];
		dz = btree[im].center[2] - leaf[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[jm].width, btree[im].width, dist);

		if ( 1 == flag ) {
		}
		else if (0 == flag) {
			walk_task_p2p(btree[im].son[0], jm);
			walk_task_p2p(btree[im].son[1], jm);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}

	if ( im >= first_node && jm >= first_node ) {
		dx = btree[im].center[0] - btree[jm].center[0];
		dy = btree[im].center[1] - btree[jm].center[1];
		dz = btree[im].center[2] - btree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(btree[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {

		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
					> btree[jm].width[0]+btree[jm].width[1] + btree[jm].width[2] ) 
			{
				walk_task_p2p(btree[im].son[0], jm);
				walk_task_p2p(btree[im].son[1], jm);
			}
			else {
				walk_task_p2p(im, btree[jm].son[0]);
				walk_task_p2p(im, btree[jm].son[1]);
			}
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}
}


void *task_compute_m2l(void *arg);

void turn2compute_m2l(){
	if (p1st == 1) {
		p1st = 0;
	} 
	else {
		pthread_join(tid, &status);
	}
	ntask[alt] = idxtask;

	argt[0] = alt;
	argt[1] = ntask[alt];

	pthread_create(&tid, NULL, task_compute_m2l, (void*)argt);

	alt = (alt+1)%2;

	ts = task_s[alt];
	tt = task_t[alt];

	idxM2L +=idxtask;
	idxtask = 0;

}



void walk_task_m2l(int im, int jm)
{
	if ( -1 == im  || -1 == jm)
		return;

	if ( im < first_leaf || jm < first_leaf ) {
		printf("error\n");
		exit(0);
	}

	if ( im == jm ) {

//		if ( im < last_leaf ) // leaf 
//		{}

		if ( im >= first_node ) 
		{
			walk_task_m2l(btree[im].son[0], btree[jm].son[0]);
			walk_task_m2l(btree[im].son[0], btree[jm].son[1]);
			walk_task_m2l(btree[im].son[1], btree[jm].son[0]);
			walk_task_m2l(btree[im].son[1], btree[jm].son[1]);
		}

		return ;
	}
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < last_leaf && jm < last_leaf ) {
		return;
	}

	if ( im < first_node && jm >= first_node ) 
	{
		dx = leaf[im].center[0] - btree[jm].center[0];
		dy = leaf[im].center[1] - btree[jm].center[1];
		dz = leaf[im].center[2] - btree[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
			if ( idxtask == LEN_TASK ) {
				turn2compute_m2l();
			}
		}
		else if (0 == flag) {
			walk_task_m2l(im, btree[jm].son[0]);
			walk_task_m2l(im, btree[jm].son[1]);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}

		return ;


	}	

	if ( im >= first_node && jm < first_node ) {
		dx = btree[im].center[0] - leaf[jm].center[0];
		dy = btree[im].center[1] - leaf[jm].center[1];
		dz = btree[im].center[2] - leaf[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[jm].width, btree[im].width, dist);

		if ( 1 == flag ) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
			if ( idxtask == LEN_TASK ) {
				turn2compute_m2l();
			}
		}
		else if (0 == flag) {
			walk_task_m2l(btree[im].son[0], jm);
			walk_task_m2l(btree[im].son[1], jm);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;

	}


	if ( im >= first_node && jm >= first_node ) {
		dx = btree[im].center[0] - btree[jm].center[0];
		dy = btree[im].center[1] - btree[jm].center[1];
		dz = btree[im].center[2] - btree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(btree[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
			if ( idxtask == LEN_TASK ) {
				turn2compute_m2l();
			}

		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
					> btree[jm].width[0]+btree[jm].width[1] + btree[jm].width[2] ) 
			{
				walk_task_m2l(btree[im].son[0], jm);
				walk_task_m2l(btree[im].son[1], jm);
			}
			else {
				walk_task_m2l(im, btree[jm].son[0]);
				walk_task_m2l(im, btree[jm].son[1]);
			}
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return;
	}
}


static double dTlocal, dTremote, dTmirror;
static INT64 NTASKP2P;
static INT64 NTASKM2L;

void fmm_prepare() 
{
	int rank, d;
	int n, mi, mj, mk;
	double t0, t1, t2, t3, t4, t5, t6, t7, t8;
	double dTlocal, dTremote, dTmirror;
	dtime_prep = dtime();

	dtime_p2p_remote = 0.0;
	dtime_p2p_mirror = 0.0;

	t0 = dtime();
	m2l_count = 0;
	double bdl[3], bdr[3];

	for (d=0; d<3; d++) {
		bdl[d] = toptree[this_domain].center[d] - 0.5*toptree[this_domain].width[d];
		bdr[d] = toptree[this_domain].center[d] + 0.5*toptree[this_domain].width[d];
	}

	build_localtree();

	for (n=first_leaf; n<last_leaf; n++)
		p2m(leaf[n].ipart, leaf[n].npart, leaf[n].center, leaf[n].M);

	walk_m2m(first_node);

	int idx[8];

	connect_local_toptree();


	width_this_domain = toptree[this_domain].width[0];

	if ( width_this_domain < toptree[this_domain].width[1])
		width_this_domain = toptree[this_domain].width[1];

	if ( width_this_domain  < toptree[this_domain].width[2])
		width_this_domain = toptree[this_domain].width[2];

	ExtDomain = (int*)pmalloc(sizeof(ExtDomain)*PROC_SIZE, 32);

	for (n=0; n<PROC_SIZE; n++)
		ExtDomain[n] = 0;

	dtime_prep = dtime() - dtime_prep;

	t1 = dtime();
	dtime_fmm = t1 - t0;
	idxM2L = 0; 
	idxP2P = 0;
}

void task_prepare_p2p() {
	alt = 0;
	p1st = 1;
	idxtask = 0;

	ts = task_s[alt];
	tt = task_t[alt];

	walk_task_p2p(first_node, first_node);
}


void task_prepare_m2l() {
	alt = 0;
	p1st = 1;
	idxtask = 0;

	ts = task_s[alt];
	tt = task_t[alt];

	walk_task_m2l(first_node, first_node);
}


void *task_compute_p2p(void *argt) {
	int n;
	int *par = argt;
	int c = par[0];
	int nt =  par[1];
	double t0 = dtime();


	double rs = splitRadius;

	double coeff = 2.0/sqrt(M_PI);

//#pragma omp parallel for
	for (n=0; n<nt; n++) {

		int inode = task_t[c][n];
		int jnode = task_s[c][n];

		int p = n*MAXLEAF;

		double dx[3], dr, ir3, ir, mp, x2;
		int ip, jp;

		int n,m, idx;


		int q;
		for (q=0; q<leaf[inode].npart; q++)
		{
			task_accx[c][p+q] = 0.0;
			task_accy[c][p+q] = 0.0;
			task_accz[c][p+q] = 0.0;
			ip = leaf[inode].ipart + q;
			for (jp=leaf[jnode].ipart; jp<leaf[jnode].ipart+leaf[jnode].npart; jp++)
			{
				if (jp == ip)
					continue;

				dx[0] = part[jp].pos[0] - part[ip].pos[0];
				dx[1] = part[jp].pos[1] - part[ip].pos[1];
				dx[2] = part[jp].pos[2] - part[ip].pos[2];

				dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
				if (dr < SoftenScale)
					ir3 = MASSPART/(SoftenScale*SoftenScale*SoftenScale);
				else
					ir3 = MASSPART/(dr*dr*dr);


#ifdef LONGSHORT
				double drs = 0.5*dr/rs;
				ir3 *= (erfc(drs) + coeff*drs*exp(-drs*drs));
#endif

				task_accx[c][p+q] +=  dx[0] * ir3;
				task_accy[c][p+q] +=  dx[1] * ir3;
				task_accz[c][p+q] +=  dx[2] * ir3;

			}
		}

	}

	for (n=0; n<nt; n++) {
		int inode = task_t[c][n];
		int q, ip = leaf[inode].ipart;
		for (q=0; q<leaf[inode].npart; q++) {

			part[ip+q].acc[0] += task_accx[c][n*MAXLEAF+q];
			part[ip+q].acc[1] += task_accy[c][n*MAXLEAF+q];
			part[ip+q].acc[2] += task_accz[c][n*MAXLEAF+q];
		}
	}


	dtime_p2p += dtime() - t0;
}


void *task_compute_m2l(void *argt) {
	int n;
	int *par = argt;
	int c = par[0];
	int nt =  par[1];
	double t0 = dtime();

	for (n=0; n<nt; n++) {
		int im = task_t[c][n];
		int jm = task_s[c][n];
		double dx, dy, dz;
		//   omp_set_lock(&writelock);
		if ( im < first_node ) {
			dx = leaf[im].center[0] - btree[jm].center[0];
			dy = leaf[im].center[1] - btree[jm].center[1];
			dz = leaf[im].center[2] - btree[jm].center[2];
			m2l(dx, dy ,dz, btree[jm].M, leaf[im].L);
		}
		else if ( jm < first_node ) {
			dx = btree[im].center[0] - leaf[jm].center[0];
			dy = btree[im].center[1] - leaf[jm].center[1];
			dz = btree[im].center[2] - leaf[jm].center[2];
			m2l(dx, dy ,dz, leaf[jm].M, btree[im].L);
		}
		else {
			dx = btree[im].center[0] - btree[jm].center[0];
			dy = btree[im].center[1] - btree[jm].center[1];
			dz = btree[im].center[2] - btree[jm].center[2];
			m2l(dx, dy ,dz, btree[jm].M, btree[im].L);
		}	
	}
	dtime_m2l += dtime() - t0;
}

void fmm_task() {
	int n;
	double t1  = dtime();
	LEN_TASK = 16384;

	task_s[0] = (int*)pmalloc(sizeof(int)*LEN_TASK,71);
	task_t[0] = (int*)pmalloc(sizeof(int)*LEN_TASK,72);

	task_s[1] = (int*)pmalloc(sizeof(int)*LEN_TASK,73);
	task_t[1] = (int*)pmalloc(sizeof(int)*LEN_TASK,74);

	task_accx[0] = (double*)pmalloc(sizeof(double)*LEN_TASK*MAXLEAF, 75);
	task_accy[0] = (double*)pmalloc(sizeof(double)*LEN_TASK*MAXLEAF, 76);
	task_accz[0] = (double*)pmalloc(sizeof(double)*LEN_TASK*MAXLEAF, 77);

	task_accx[1] = (double*)pmalloc(sizeof(double)*LEN_TASK*MAXLEAF, 75);
	task_accy[1] = (double*)pmalloc(sizeof(double)*LEN_TASK*MAXLEAF, 76);
	task_accz[1] = (double*)pmalloc(sizeof(double)*LEN_TASK*MAXLEAF, 77);



	task_prepare_p2p();

	if (p1st != 1) {

		pthread_join(tid, &status);

	}
	idxP2P +=idxtask;

	ntask[alt] = idxtask;

	argt[0] = alt;
	argt[1] = ntask[alt];


	pthread_create(&tid, NULL, task_compute_p2p, argt);
	pthread_join(tid, &status);

	double t0 = dtime();
	task_prepare_m2l();

	if (p1st != 1) {
		pthread_join(tid, &status);
	}
	idxM2L +=idxtask;

	ntask[alt] = idxtask;

	argt[0] = alt;
	argt[1] = ntask[alt];


	pthread_create(&tid, NULL, task_compute_m2l, argt);
	pthread_join(tid, &status);



	pfree(task_s[0],71);
	pfree(task_t[0],72);

	pfree(task_s[1],73);
	pfree(task_t[1],74);

	pfree(task_accx[0], 75);
	pfree(task_accy[0], 76);
	pfree(task_accz[0], 77);

	pfree(task_accx[1], 75);
	pfree(task_accy[1], 76);
	pfree(task_accz[1], 77);

	dtime_task=dtime() - t1;

	if (0==PROC_RANK) printf(" fmm time = %lf,task_m2l=%lf \n", dtime() - t1,dtime() - t0);
}


void fmm_ext(){
	int rank, d;
	int n, mi, mj, mk;
	double t0, t1, t2, t3, t4, t5, t6, t7, t8;
	double dTlocal, dTremote, dTmirror;

	t0 = dtime();

	dtime_ext = dtime();

	double lenExTree = 3.0/MAXLEAF;
	double lenExBody = 4.5;

	lenExTree *= ((double)NPART_TOTAL/(double)PROC_SIZE) * sizeof(RemoteNode) ;
	lenExBody *= ((double)NPART_TOTAL/(double)PROC_SIZE) * sizeof(RemoteBody) ;


	lenExTree = (NNODE*1.1 * sizeof(RemoteNode) );
	lenExBody = (NPART*1.1 * sizeof(RemoteBody) );


	exstree = (RemoteNode*)pmalloc((size_t)lenExTree,71);
	exsbody = (RemoteBody*)pmalloc((size_t)lenExBody,72);
	exrtree = (RemoteNode*)pmalloc((size_t)lenExTree,73);
	exrbody = (RemoteBody*)pmalloc((size_t)lenExBody,74);

	t4 = dtime();

	MPI_Barrier(MPI_COMM_WORLD);

	double shift[3] = {0.0, 0.0, 0.0};

	t5 = dtime();

	for (n=1; n<PROC_SIZE; n++) {
		fmm_remote(n, shift);
	}

	t6 = dtime();


#ifdef PERIODIC_CONDITION
	for (mi=-1; mi<=1; mi++) {	
		for (mj=-1; mj<=1; mj++) {	
			for (mk=-1; mk<=1; mk++) {
				if ( 0== mi && 0== mj && 0 == mk)
					continue;

				shift[0] = (double)mi*BOXSIZE;	
				shift[1] = (double)mj*BOXSIZE;	
				shift[2] = (double)mk*BOXSIZE;	

				for (n=0; n<PROC_SIZE; n++) {
					fmm_remote(n, shift);
				}
			}
		}
	}
#endif

	t7 = dtime();

	pfree( exstree , 71);
	pfree( exsbody , 72);
	pfree( exrtree , 73);
	pfree( exrbody , 74);
	MPI_Barrier(MPI_COMM_WORLD);
	walk_l2l(first_node);

	for (n=first_leaf; n<last_leaf; n++)
		l2p(n);

	t8 = dtime();

	dtime_fmm += (t8 - t7);

	double dTm2l = (t3-t2);
	double dTp2p = (t4-t3);
	dTremote = (t6-t5);
	dTmirror = (t7-t6);


	DTIME_THIS_DOMAIN = idxP2P + idxM2L;

	dtime_fmm_remote = t7 - t5;

	dtime_ext = dtime() - dtime_ext;

}

void fmm_deconstruct() 
{
	btree += first_node;
	pfree(btree, 31);

	leaf  += first_leaf;
	pfree(leaf, 30);

	pfree(ExtDomain, 32);
	deconstruct_toptree();
}

