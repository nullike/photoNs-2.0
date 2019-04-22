/*
 * photoNs-2
 *
 *      2017 - 12 - 8
 *	qwang@nao.cas.cn
 */	 

#include "photoNs.h"


void connect_local_toptree()
{
	int n;
	int  idx = PROC_RANK + mostleft;
	if (idx > last_domain)
		idx -= PROC_SIZE;

	idx = this_domain;

	// optional
	toptree[idx].width[0] = btree[first_node].width[0];
	toptree[idx].width[1] = btree[first_node].width[1];
	toptree[idx].width[2] = btree[first_node].width[2];

	toptree[idx].center[0] = btree[first_node].center[0];
	toptree[idx].center[1] = btree[first_node].center[1];
	toptree[idx].center[2] = btree[first_node].center[2];
	// optional

	for (n=0; n<NMULTI; n++) {
		toptree[idx].M[n] = btree[first_node].M[n];
	}

	TopNode *recvbuf = (TopNode*)pmalloc(sizeof(TopNode)*PROC_SIZE, 26);

	MPI_Allgather( toptree+idx, sizeof(TopNode), MPI_BYTE, recvbuf, sizeof(TopNode), MPI_BYTE, MPI_COMM_WORLD);

	for (n=0; n<PROC_SIZE; n++) {
		idx = ( n +  mostleft );
		if (idx > last_domain)
			idx -= PROC_SIZE;


		toptree[idx] = recvbuf[n];
	}

	pfree(recvbuf, 26);

	walk_toptree_m2m(0);
}



void construct_toptree()
{

	int n, m,  length = 2*PROC_SIZE - 1;
	toptree = (TopNode*)pmalloc(sizeof(TopNode) * length, 27);

	for (n=0; n<length; n++) {
		toptree[n].width[0] = 0.0;
		toptree[n].width[1] = 0.0;
		toptree[n].width[2] = 0.0;

		toptree[n].center[0] = 0.0;
		toptree[n].center[1] = 0.0;
		toptree[n].center[2] = 0.0;

		for (m=0; m<NMULTI; m++)
			toptree[n].M[m] = 0.0;

		toptree[n].split = 0.0;
	}

	for (n=0; n<first_domain; n++) {
		toptree[n].son[0] = 2*n+1;
		toptree[n].son[1] = 2*n+2;

	}

	for (n=first_domain; n<=last_domain; n++) {
		toptree[n].son[0] = -1;
		toptree[n].son[1] = -1;
	}

}


void clean_toptree()  {
	int n, m,  length = 2*PROC_SIZE - 1;

	for (n=0; n<length; n++)  {
		toptree[n].width[0] = 0.0;
		toptree[n].width[1] = 0.0;
		toptree[n].width[2] = 0.0;

		toptree[n].center[0] = 0.0;
		toptree[n].center[1] = 0.0;
		toptree[n].center[2] = 0.0;

		for (m=0; m<NMULTI; m++)
			toptree[n].M[m] = 0.0;

		toptree[n].split = 0.0;
	}

	for (n=0; n<first_domain; n++) {
		toptree[n].son[0] = 2*n+1;
		toptree[n].son[1] = 2*n+2;
	}
	for (n=first_domain; n<=last_domain; n++) {
		toptree[n].son[0] = -1;
		toptree[n].son[1] = -1;
	}

}

void deconstruct_toptree() {
	pfree(toptree, 27);
}



#ifdef CHECK_TOPTREE
int walk_part(int D, int iNode, double x[3]) {
	if (iNode >=first_domain)
		return iNode - first_domain;

	if (x[D] > toptree[iNode].split &&toptree[iNode].son[1] > -1 )
		return walk_part((D+1)%3, toptree[iNode].son[1],  x) ;
	if (x[D] <= toptree[iNode].split && toptree[iNode].son[0] > -1 )
		return walk_part((D+1)%3, toptree[iNode].son[0],  x) ;
}

void check(int npart) {
	int n;

	int count[3] = { 0, 0, 0};

	for (n=0; n<npart; n++) {
		count[walk_part(0, 0, part[n].pos)] ++ ;
	}
	for (n=0; n<PROC_SIZE; n++)
		printf(" check count [%d] = %d at rank %d\n", n, count[n], PROC_RANK);

}
#endif


void center_toptree(int direct, int iNode, double left[3], double right[3])
{

	toptree[iNode].width[0] = right[0]-left[0];
	toptree[iNode].width[1] = right[1]-left[1];
	toptree[iNode].width[2] = right[2]-left[2];

	toptree[iNode].center[0] = 0.5*(right[0]+left[0]);
	toptree[iNode].center[1] = 0.5*(right[1]+left[1]);
	toptree[iNode].center[2] = 0.5*(right[2]+left[2]);

	if (iNode == this_domain )
		direct_local_start = direct;

	int newd = (direct + 1)%3;
	double tmp;


	if (toptree[iNode].son[0] > -1) {
		tmp = right[direct];
		right[direct] = toptree[iNode].split;
		center_toptree(newd, toptree[iNode].son[0], left, right);
		right[direct] =tmp;
	}

	if (toptree[iNode].son[1] > -1) {
		tmp = left[direct];
		left[direct] = toptree[iNode].split;
		center_toptree(newd, toptree[iNode].son[1], left, right);
		left[direct] =tmp;
	}

}


void walk_toptree_m2m(int iNode) {
	int n, idx;
	double dx, dy, dz;

	for (n=0; n<NSON; n++) {
		idx = toptree[iNode].son[n];

		if (idx > -1) {
			walk_toptree_m2m(idx);

			dx =  toptree[iNode].center[0] - toptree[idx].center[0];
			dy =  toptree[iNode].center[1] - toptree[idx].center[1];
			dz =  toptree[iNode].center[2] - toptree[idx].center[2];

			m2m(dx, dy, dz, toptree[idx].M, toptree[iNode].M);
		}
	}
}

void walk_toptree_m2l(int iNode, double L[]) 
{
	int n, idx;
	double dx, dy, dz, dr;
	double max;

	if (iNode == this_domain)
		return;

	dx = toptree[this_domain].center[0] - toptree[iNode].center[0];
	dy = toptree[this_domain].center[1] - toptree[iNode].center[1];
	dz = toptree[this_domain].center[2] - toptree[iNode].center[2];

	dr = sqrt(dx*dx + dy*dy +dz*dz);

	max = toptree[iNode].width[0];

	if (max < toptree[iNode].width[1])
		max = toptree[iNode].width[1];

	if (max < toptree[iNode].width[2])
		max = toptree[iNode].width[2];

	if (max < open_angle * dr) {
		m2l(dx, dy, dz, toptree[iNode].M, L);

		return ;
	}
	else
		for (n=0; n<NSON; n++) {
			if (toptree[iNode].son[n] > -1) {
				if (toptree[iNode].son[n] < first_domain)
					walk_toptree_m2l(toptree[iNode].son[n], L);
				else if (toptree[iNode].son[n] >= first_domain ) {
					idx = (toptree[iNode].son[n]-mostleft+PROC_SIZE)%PROC_SIZE;
					if (idx != PROC_RANK)
						ExtDomain[idx] = 1;
				}
			}
		}
}




