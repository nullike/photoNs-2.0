/*
 * photoNs-2
 *
 *      2018 - 8 - 12
 *	qwang@nao.cas.cn
 */	 
#include <pthread.h>
#include "photoNs.h"
#include "remotes.h"

static int numbody;
static int numnode;

void p2p_kernel_ex(int inode, int jnode)
{
	double dx[3], ir3, ir, dr,x2;
	int ip, jp;

	double rs = splitRadius;	
	double coeff = 2.0/sqrt(M_PI);

	int n,m, idx;
	int jstart = exrtree[jnode].son[0];
	int jend   = exrtree[jnode].son[0] + exrtree[jnode].npart;

	for (ip=leaf[inode].ipart; ip<leaf[inode].ipart+leaf[inode].npart; ip++)
	{
		for (jp=jstart; jp<jend; jp++)
		{
			dx[0] = exrbody[jp].pos[0] - part[ip].pos[0];
			dx[1] = exrbody[jp].pos[1] - part[ip].pos[1];
			dx[2] = exrbody[jp].pos[2] - part[ip].pos[2];


			dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
			if (dr < SoftenScale)
				ir3 = MASSPART/(SoftenScale*SoftenScale*SoftenScale);
			else
				ir3 = MASSPART/(dr*dr*dr);


#ifdef LONGSHORT
			double drs = 0.5*dr/rs;
			ir3 *= (erfc(drs) + coeff*drs*exp(-drs*drs));
#endif


			part[ip].acc[0] +=  dx[0] * ir3;
			part[ip].acc[1] +=  dx[1] * ir3;
			part[ip].acc[2] +=  dx[2] * ir3;

			p2p_count_remote++;

		}

	}
}


void prepare_sendtree2(int isend, int ilocal, int tNode, int D, double displace[3])
{
	int n, idx, p;
	double dx, dy, dz, dr;

	if (ilocal < first_node && ilocal >= first_leaf ) {
		int ileaf = ilocal;

		exstree[isend].npart = leaf[ileaf].npart;
		exstree[isend].center[0] = leaf[ileaf].center[0] + displace[0];
		exstree[isend].center[1] = leaf[ileaf].center[1] + displace[1];
		exstree[isend].center[2] = leaf[ileaf].center[2] + displace[2];

		exstree[isend].width[0] = leaf[ileaf].width[0];
		exstree[isend].width[1] = leaf[ileaf].width[1];
		exstree[isend].width[2] = leaf[ileaf].width[2];

		for (p=0; p<NMULTI; p++) {
			exstree[isend].M[p] = leaf[ileaf].M[p];
		}

		exstree[isend].son[0] = numbody;

		for (p=leaf[ileaf].ipart; p<leaf[ileaf].ipart+leaf[ileaf].npart; p++) {
			exsbody[numbody].pos[0] = part[p].pos[0] + displace[0];
			exsbody[numbody].pos[1] = part[p].pos[1] + displace[1];
			exsbody[numbody].pos[2] = part[p].pos[2] + displace[2];
			//exsbody[numbody].mass = part[p].mass;
			numbody++;
		}

		exstree[isend].son[1] = numbody;
		numnode ++;

		return;
	}

	dx = toptree[tNode].center[0] - btree[ilocal].center[0] - displace[0];
	dy = toptree[tNode].center[1] - btree[ilocal].center[1] - displace[1];
	dz = toptree[tNode].center[2] - btree[ilocal].center[2] - displace[2];

	if (dx < 0.0)
		dx = -dx;
	if (dy < 0.0)
		dy = -dy;
	if (dz < 0.0)
		dz = -dz;

	dx -= (toptree[tNode].width[0] + btree[ilocal].width[0])*0.5;
	dy -= (toptree[tNode].width[1] + btree[ilocal].width[1])*0.5;
	dz -= (toptree[tNode].width[2] + btree[ilocal].width[2])*0.5;

	dr = 0.0;
	if (dx > 0.0)
		dr += dx*dx;
	if (dy > 0.0)
		dr += dy*dy;
	if (dz > 0.0)
		dr += dz*dz;
	dr = sqrt(dr);

	exstree[isend].npart     = btree[ilocal].npart;

	exstree[isend].center[0] = btree[ilocal].center[0] + displace[0];
	exstree[isend].center[1] = btree[ilocal].center[1] + displace[1];
	exstree[isend].center[2] = btree[ilocal].center[2] + displace[2];

	exstree[isend].width[0]  = btree[ilocal].width[0];
	exstree[isend].width[1]  = btree[ilocal].width[1];
	exstree[isend].width[2]  = btree[ilocal].width[2];

	for (p=0; p<NMULTI; p++) {
		exstree[isend].M[p] = btree[ilocal].M[p];
	}

	numnode ++;

	double width_max;
	width_max = btree[ilocal].width[0];
	if (width_max < btree[ilocal].width[1])    
		width_max = btree[ilocal].width[1];
	if (width_max < btree[ilocal].width[2])    
		width_max = btree[ilocal].width[2];


#ifdef LONGSHORT
	if (dr >= cutoffRadius) {
		exstree[isend].son[0] = -1;
		exstree[isend].son[1] = -1;
		return;
	}
#endif


	if ( width_max < 0.95*open_angle*dr ) {
		exstree[isend].son[0] = -1;
		exstree[isend].son[1] = -1;
		return;
	}
	else {
		for (n=0; n<NSON; n++) {
			idx = btree[ilocal].son[n] ;
			if (idx >= first_leaf) {
				exstree[isend].son[n] = numnode;
				prepare_sendtree2(numnode, idx, tNode, (D+1)%3, displace);
			} 
		}
	}

}


static int *task_s_ex[2];
static int *task_t_ex[2];
static int ntask_ex[2];

static int idxtask_ex;
static int LEN_TASK;

int *ts_ex ;
int *tt_ex ;
pthread_t tid_ex ;
int p1st_ex = 1;
void *status_ex ;
int alt_ex ;
int argt_ex [2];

void *task_compute_p2p_ext(void *arg);

void turn2compute_p2p_ext(){
	if (p1st_ex == 1) {
		p1st_ex = 0;
	} 
	else {
		pthread_join(tid_ex , &status_ex );
	}
	ntask_ex[alt_ex] = idxtask_ex;

	argt_ex[0] = alt_ex ;
	argt_ex[1] = ntask_ex[alt_ex ];

	pthread_create(&tid_ex , NULL, task_compute_p2p_ext, (void*)argt_ex );

	alt_ex  = (alt_ex +1)%2;

	ts_ex  = task_s_ex[alt_ex ];
	tt_ex  = task_t_ex[alt_ex ];

	idxtask_ex = 0;

}


void walk_task_p2p_ext(int im, int jm)
{
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < first_leaf ) {
		printf(" error1 \n");
		exit(0);
	}
	if ( jm < 0 ) {
		printf(" error3 im = %d jm = %d\n", im, jm);
		exit(0);
	}

	// copy leaf
	if ( im < first_node && exrtree[jm].npart <= MAXLEAF ) 
	{

		*(ts_ex  + idxtask_ex) = jm;
		*(tt_ex  + idxtask_ex) = im;
		idxtask_ex ++;

		if ( idxtask_ex == LEN_TASK ) {
			turn2compute_p2p_ext();
		}

		return;
	}


	if ( im < first_node && exrtree[jm].npart > MAXLEAF ) 
	{
		dx = leaf[im].center[0] - exrtree[jm].center[0];
		dy = leaf[im].center[1] - exrtree[jm].center[1];
		dz = leaf[im].center[2] - exrtree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag || exrtree[jm].son[0]<0 || exrtree[jm].son[1] < 0) {

			//			*(ts + idxtask) = jm;
			//			*(tt + idxtask) = im;
			//			idxtask ++;

			//	taskM2Lexts[idxM2Lext] = jm;
			//	taskM2Lextt[idxM2Lext] = im;
			//	idxM2Lext ++;

			//        	m2l(dx, dy, dz, exrtree[jm].M, leaf[im].L);
		}
		else if (0 == flag) {
			walk_task_p2p_ext(im, exrtree[jm].son[0]);
			walk_task_p2p_ext(im, exrtree[jm].son[1]);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}

		return ;

	}

	if ( im >= first_node && exrtree[jm].npart <= MAXLEAF ) {
		dx = btree[im].center[0] - exrtree[jm].center[0];
		dy = btree[im].center[1] - exrtree[jm].center[1];
		dz = btree[im].center[2] - exrtree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(btree[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag) {
			//	taskM2Lexts[idxM2Lext] = jm;
			//	taskM2Lextt[idxM2Lext] = im;
			//	idxM2Lext ++;

			//			*(ts + idxtask) = jm;
			//			*(tt + idxtask) = im;
			//			idxtask ++;

			//			m2l(dx, dy, dz, exrtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			walk_task_p2p_ext(btree[im].son[0], jm);
			walk_task_p2p_ext(btree[im].son[1], jm);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}


		return;
	}

	//printf(" %d %lf\n", exrtree[jm].npart, exrtree[jm].M[2]);
	if ( im >= first_node && exrtree[jm].npart > MAXLEAF ) 
	{
		dx = btree[im].center[0] - exrtree[jm].center[0];
		dy = btree[im].center[1] - exrtree[jm].center[1];
		dz = btree[im].center[2] - exrtree[jm].center[2];

		r2 = dx*dx + dy*dy + dz*dz ;

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag ) {

			//			*(ts + idxtask) = jm;
			//			*(tt + idxtask) = im;
			//			idxtask ++;
			//	taskM2Lexts[idxM2Lext] = jm;
			//	taskM2Lextt[idxM2Lext] = im;
			//	idxM2Lext ++;


			//			m2l(dx, dy, dz, exrtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1]+btree[im].width[2] 
					> exrtree[jm].width[0]+exrtree[jm].width[1]+exrtree[jm].width[2]
					|| exrtree[jm].son[0] < 0 || exrtree[jm].son[1] < 0 ) 
			{
				walk_task_p2p_ext(btree[im].son[0], jm);
				walk_task_p2p_ext(btree[im].son[1], jm);
			}
			else {
				walk_task_p2p_ext(im, exrtree[jm].son[0]);
				walk_task_p2p_ext(im, exrtree[jm].son[1]);
			}


		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return;

	}
}

void *task_compute_m2l_ext(void *arg);

void turn2compute_m2l_ext(){
	if (p1st_ex == 1) {
		p1st_ex = 0;
	} 
	else {
		pthread_join(tid_ex , &status_ex );
	}
	ntask_ex[alt_ex ] = idxtask_ex;

	argt_ex[0] = alt_ex ;
	argt_ex[1] = ntask_ex[alt_ex ];
	pthread_create(&tid_ex , NULL, task_compute_m2l_ext, (void*)argt_ex );
	//last join (if not the first)
	//fork current
	//exchange & reset

	//	task_compute_p2p(NULL);
	// printf(" turn %d, idxtas = %d\n", alt, idxtask);
	alt_ex  = (alt_ex +1)%2;

	ts_ex  = task_s_ex[alt_ex ];
	tt_ex  = task_t_ex[alt_ex ];

	idxtask_ex = 0;

}



void walk_task_m2l_ext(int im, int jm)
{
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < first_leaf ) {
		printf(" error1 \n");
		exit(0);
	}
	if ( jm < 0 ) {
		printf(" error3 im = %d jm = %d\n", im, jm);
		exit(0);
	}

	// copy leaf
	if ( im < first_node && exrtree[jm].npart <= MAXLEAF ) 
	{

		return;
	}


	if ( im < first_node && exrtree[jm].npart > MAXLEAF ) 
	{
		dx = leaf[im].center[0] - exrtree[jm].center[0];
		dy = leaf[im].center[1] - exrtree[jm].center[1];
		dz = leaf[im].center[2] - exrtree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag || exrtree[jm].son[0]<0 || exrtree[jm].son[1] < 0) {

			*(ts_ex  + idxtask_ex) = jm;
			*(tt_ex  + idxtask_ex) = im;
			idxtask_ex ++;

			if ( idxtask_ex == LEN_TASK ) {
				turn2compute_m2l_ext();
			}

		}
		else if (0 == flag) {
			walk_task_m2l_ext(im, exrtree[jm].son[0]);
			walk_task_m2l_ext(im, exrtree[jm].son[1]);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}

		return ;

	}

	if ( im >= first_node && exrtree[jm].npart <= MAXLEAF ) {
		dx = btree[im].center[0] - exrtree[jm].center[0];
		dy = btree[im].center[1] - exrtree[jm].center[1];
		dz = btree[im].center[2] - exrtree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag) {

			*(ts_ex  + idxtask_ex) = jm;
			*(tt_ex  + idxtask_ex) = im;
			idxtask_ex ++;

			if ( idxtask_ex == LEN_TASK ) {
					turn2compute_m2l_ext();
			}
		}
		else if (0 == flag) {
			walk_task_m2l_ext(btree[im].son[0], jm);
			walk_task_m2l_ext(btree[im].son[1], jm);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}


		return;
	}


	if ( im >= first_node && exrtree[jm].npart > MAXLEAF ) 
	{
		dx = btree[im].center[0] - exrtree[jm].center[0];
		dy = btree[im].center[1] - exrtree[jm].center[1];
		dz = btree[im].center[2] - exrtree[jm].center[2];

		r2 = dx*dx + dy*dy + dz*dz ;

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag ) {

			*(ts_ex  + idxtask_ex) = jm;
			*(tt_ex  + idxtask_ex) = im;
			idxtask_ex ++;
			if ( idxtask_ex == LEN_TASK ) {
		
				turn2compute_m2l_ext();
			}

		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1]+btree[im].width[2] 
					> exrtree[jm].width[0]+exrtree[jm].width[1]+exrtree[jm].width[2]
					|| exrtree[jm].son[0] < 0 || exrtree[jm].son[1] < 0 ) 
			{
				walk_task_m2l_ext(btree[im].son[0], jm);
				walk_task_m2l_ext(btree[im].son[1], jm);
			}
			else {
				walk_task_m2l_ext(im, exrtree[jm].son[0]);
				walk_task_m2l_ext(im, exrtree[jm].son[1]);
			}


		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return;

	}
}



void task_prepare_p2p_ext() {
	alt_ex  = 0;
	p1st_ex = 1;
	idxtask_ex = 0;

	ts_ex  = task_s_ex[alt_ex];
	tt_ex  = task_t_ex[alt_ex];

	walk_task_p2p_ext(first_node, 0);


}


void task_prepare_m2l_ext() {
	alt_ex  = 0;
	p1st_ex = 1;
	idxtask_ex = 0;

	ts_ex  = task_s_ex[alt_ex ];
	tt_ex  = task_t_ex[alt_ex ];

	walk_task_m2l_ext(first_node, 0);


}

void *task_compute_p2p_ext(void *argt_ex ) {
	int n;
	int *par = argt_ex ;
	int c = par[0];
	int nt =  par[1];
	double t0 = dtime();

	for (n=0; n<nt; n++) {
		int im = task_t_ex[c][n];
		int jm = task_s_ex[c][n];
		p2p_kernel_ex(im, jm);
	}
	dtime_p2p += dtime() - t0;
}

void *task_compute_m2l_ext(void *argt_ex ) {
	int n;
	int *par = argt_ex ;
	int c = par[0];
	int nt =  par[1];
	double t0 = dtime();

	for (n=0; n<nt; n++) {
		int im = task_t_ex[c][n];
		int jm = task_s_ex[c][n];
		double dx, dy, dz;
	
		if ( im < first_node) {
			dx = leaf[im].center[0] - exrtree[jm].center[0];
			dy = leaf[im].center[1] - exrtree[jm].center[1];
			dz = leaf[im].center[2] - exrtree[jm].center[2];

			m2l(dx, dy ,dz, exrtree[jm].M, leaf[im].L);
		}
		if ( im >= first_node ) {

			dx = btree[im].center[0] - exrtree[jm].center[0];
			dy = btree[im].center[1] - exrtree[jm].center[1];
			dz = btree[im].center[2] - exrtree[jm].center[2];

			m2l(dx, dy ,dz, exrtree[jm].M, btree[im].L);
		}

	}
	dtime_m2l += dtime() - t0;
}


void fmm_remote_task(int numrecv) {
	int n;
	double t1  = dtime();
	LEN_TASK = 100000;

	task_s_ex[0] = (int*)pmalloc(sizeof(int)*LEN_TASK, 81);
	task_t_ex[0] = (int*)pmalloc(sizeof(int)*LEN_TASK, 82);

	task_s_ex[1] = (int*)pmalloc(sizeof(int)*LEN_TASK, 83);
	task_t_ex[1] = (int*)pmalloc(sizeof(int)*LEN_TASK, 84);

	task_prepare_p2p_ext();

	if (p1st_ex != 1) {
		pthread_join(tid_ex , &status_ex );
	}
	ntask_ex[alt_ex ] = idxtask_ex;

	argt_ex[0] = alt_ex ;
	argt_ex[1] = ntask_ex[alt_ex ];

	pthread_create(&tid_ex , NULL, task_compute_p2p_ext, argt_ex);

	pthread_join(tid_ex, &status_ex);


	ntask_ex[0]= ntask_ex[1]= 0;
	task_prepare_m2l_ext();

	if (p1st_ex != 1) {
		pthread_join(tid_ex, &status_ex);

	}
	ntask_ex[alt_ex] = idxtask_ex;
	argt_ex[0] = alt_ex;
	argt_ex[1] = ntask_ex[alt_ex];

	pthread_create(&tid_ex, NULL, task_compute_m2l_ext, argt_ex);

	pthread_join(tid_ex, &status_ex);

	if (NULL != task_s_ex[0])
		pfree(task_s_ex[0] ,81);
	if (NULL != task_t_ex[0])
		pfree(task_t_ex[0], 82);

	if (NULL != task_s_ex[1])
		pfree(task_s_ex[1], 83);
	if (NULL != task_t_ex[1])
		pfree(task_t_ex[1], 84);

}

void fmm_remote(int idx, double displace[3])
{
	int srank, d, rrank;
	int n, mi, mj, mk;


	srank = (PROC_RANK + idx ) % PROC_SIZE;
	rrank = (PROC_RANK - idx + PROC_SIZE) % PROC_SIZE;

	

	double bdl[3], bdr[3];

	for (d=0; d<3; d++) {
		bdl[d] = toptree[this_domain].center[d] - 0.5*toptree[this_domain].width[d];
		bdr[d] = toptree[this_domain].center[d] + 0.5*toptree[this_domain].width[d];
	}

	width_this_domain = toptree[this_domain].width[0];

	if ( width_this_domain < toptree[this_domain].width[1])
		width_this_domain = toptree[this_domain].width[1];

	if ( width_this_domain < toptree[this_domain].width[2])
		width_this_domain = toptree[this_domain].width[2];

	int sendnumnode, sendnumbody;
	int recvnumnode, recvnumbody;

	int tNode = srank + mostleft ;
	if (tNode > last_domain)
		tNode -= PROC_SIZE;

	numbody = 0;
	numnode = 0;

	prepare_sendtree2(numnode, first_node, tNode, direct_local_start, displace);

	sendnumnode = numnode ;
	sendnumbody = numbody ;

	MPI_Request request, request2;
	MPI_Status status, status2;

	MPI_Isend(&sendnumnode, 1, MPI_INT, srank, 101, MPI_COMM_WORLD, &request);
	MPI_Recv( &recvnumnode, 1, MPI_INT, rrank, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Wait(&request, &status);

	MPI_Isend(&sendnumbody, 1, MPI_INT, srank, 102, MPI_COMM_WORLD, &request2);
	MPI_Recv( &recvnumbody, 1, MPI_INT, rrank, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Wait(&request2, &status2);


	if ( 0 == recvnumbody && 0 == recvnumnode)
		return;

	MPI_Isend(exstree, sendnumnode, strReNode, srank, 111, MPI_COMM_WORLD, &request);
	MPI_Recv( exrtree, recvnumnode, strReNode, rrank, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Wait(&request, &status);    

	MPI_Isend(exsbody, sendnumbody, strReBody, srank, 112, MPI_COMM_WORLD, &request2);
	MPI_Recv( exrbody, recvnumbody, strReBody, rrank, 112, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Wait(&request2, &status2);    

	fmm_remote_task ( recvnumnode );
	//	walk_m2l_remote2(first_node, 0);

}




