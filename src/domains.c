#include "photoNs.h"
#include "domains.h"

////////////////////////////////////////////////////////////////////////////////
double fill_time_domtree(int iNode)
{
	if (iNode >= first_domain)
		return domtree[iNode].time_node;

	if (domtree[iNode].son[0] > 0)
		domtree[iNode].time_left  = fill_time_domtree(domtree[iNode].son[0]);
	if (domtree[iNode].son[1] > 0)
		domtree[iNode].time_right = fill_time_domtree(domtree[iNode].son[1]);

	domtree[iNode].time_node =  domtree[iNode].time_left +  domtree[iNode].time_right;

	return domtree[iNode].time_node;
}

void measure_domain_runtime(double dTimeFrac) {
	int n;
	double *time_rank = (double*)malloc(sizeof(double)*PROC_SIZE);

	MPI_Allgather(&dTimeFrac, 1, MPI_DOUBLE, 
			time_rank, 1, MPI_DOUBLE, MPI_COMM_WORLD);

	for (n=0; n<PROC_SIZE; n++) {
		int idom = n + mostleft;
		if (idom  > 2*PROC_SIZE - 2)
			idom  -= PROC_SIZE;

		domtree[idom].time_node = time_rank[n];
	}

	fill_time_domtree(0);

	free(time_rank);
}

////////////////////////////////////////////////////////////////////////////////

double fraction(int size, int *l, int *r) {
	int left, right;
	double frac;

	if (size == 1) {
		left  = 1;
		right = 0;
		frac = 1.0;
	}
	else if (size == 2 ) {
		left  = 1;
		right = 1;
		frac = 0.5;
	}
	else if (size == 3 ) {
		left  = 2;
		right = 1;
		frac =  0.666666667;
	}
	else {
		left = 1;
		right = 2;

		while ( size - left >= right -size ) {
			left *=2;
			right *= 2;
		}


		left >>= 1;
		right = size - left;
		if (left < right) {
			left = right;
			right = size-left;
		}
		frac = ((double)left)/((double)size);
	}


	*l = left;
	*r = right;
	return  frac;
}

void determine_split_node(int D, int nproc, int iNode, DomainNode *tree, double bl[], double br[])
{
	int nleft, nright;
	int npart[2];
	if (iNode >= first_domain)
		return;

	double frac_split =  fraction(nproc, &nleft, &nright);

	double relax = 0.3;
	double shift;

	double new_split, previous_split = domtree[iNode].split;
	double t1 = domtree[iNode].time_left/nleft;
	double t2 = domtree[iNode].time_right/nright;

	double w0l =  domtree[iNode].split - bl[D];
	double w0r =  br[D] - domtree[iNode].split;

	if ( t2 > t1) {
		shift = w0r*(1.0-(t1+t2 )/(2*t2));
		if (shift > relax * w0r) {
			shift = relax * w0r;
		}

	}
	else {
		shift =  w0l*(1.0-(t1+t2 )/(2*t1));
		if (shift > relax * w0l) {
			shift = relax * w0l;
		}
		shift *= -1;
	}

	shift = 0.5*relax*(t2-t1)/(t1*nleft/w0l+t2*nright/w0r);

	new_split = domtree[iNode].split + shift;

	double bound_left[3], bound_right[3];

	bound_left[0] = bl[0];
	bound_left[1] = bl[1];
	bound_left[2] = bl[2];

	bound_right[0] = br[0];
	bound_right[1] = br[1];
	bound_right[2] = br[2];

	bound_right[D] = tree[iNode].split;

	determine_split_node((D+1)%3, nleft,  tree[iNode].son[0], tree, bound_left, bound_right);

	bound_left[D] = tree[iNode].split;
	bound_right[D] = br[D];

	determine_split_node((D+1)%3, nright, tree[iNode].son[1], tree, bound_left, bound_right);

	tree[iNode].split = new_split;
}

void determine_split_domtree(int size, int npart, DomainNode *tree)
{

	int  n,iNode, length = 2*PROC_SIZE - 1;

	int D = 0;
	double bl[3] = {    0.0, 0.0, 0.0};
	double br[3] = {BOXSIZE, BOXSIZE, BOXSIZE};

	determine_split_node(D, PROC_SIZE,0, tree, bl, br);

}


///////////////////////////////////////////////////////////////////////////////


void bksort_body_inplace(int D, Body *p, int length, int npart[2], double split)
{
	int n, count=0;
	for (n=0; n<length; n++) {
		if ( p[n].pos[D] > split )
			count++;

	}

	if (length == 0 ) {
		npart[0] = 0;
		npart[1] = 0;
		return;
	}

	if (length == 1 ) {
		if (p[0].pos[D] > split) {
			npart[0] = 0;
			npart[1] = 1;
		}
		else {
			npart[0] = 1;
			npart[1] = 0;
		}
		return ;
	}


	Body temp;

	if (length == 2 )
	{
		if (p[0].pos[D] > p[1].pos[D]) {
			temp = p[0];
			p[0] = p[1];
			p[1] = temp;
		}
		if (p[0].pos[D] > split) {
			npart[0] = 0;
			npart[1] = 2;
		}
		else if (p[1].pos[D] <= split) {
			npart[0] = 2;
			npart[1] = 0;
		}
		else {
			npart[0] = npart[1] = 1;
		}
		return;
	}

	int top = 0;
	while (top < length && p[top].pos[D] <= split)
		top++;

	int but = length - 1;
	while (but >= 0 && p[but].pos[D] > split )
		but--;

	//printf(" top = %d  but = %d  [%d] \n", top, but, PROC_RANK);

	if (top == length) {
		npart[0] = length;
		npart[1] = 0;
		return;
	}


	if (but == -1) {
		npart[1] = length;
		npart[0] = 0;
		return;
	}   

	n=top;
	for (n=top; n<=but; n++) {
		if ( p[n].pos[D] > split ) {
			temp = p[n];
			p[n] = p[but];
			p[but] = temp;

			while ( p[but].pos[D] > split) {
				but--;
			}
		}    

	}

	if (n==but)
		npart[0] = but+1;
	else    
		npart[0] = n;

	npart[1] = length - npart[0];



	for (n=0; n<npart[0]; n++) {
		if ( p[n].pos[D] > split )
			printf("error %d [%d]\n", n, PROC_RANK);
	}

	for (n=npart[0]; n<length; n++) {
		if ( p[n].pos[D] <= split )
			printf("error %d [%d]\n", n, PROC_RANK);
	}

}

void prepare_body_inOrderOf_domain(int D,  Body *p, int length, int iNode, int send[]) 
{

	int npart[2];
	int rank;

	if (iNode >= first_domain) {
		rank = (iNode - mostleft + PROC_SIZE ) % PROC_SIZE;
		send[rank] = length;
		return;
	}


	//    npart[0] = npart[1] = 0;
	bksort_body_inplace(D, p, length, npart, domtree[iNode].split);

	if (iNode <= last_domain) {
		prepare_body_inOrderOf_domain((D+1)%3,  p, npart[0], domtree[iNode].son[0], send);
	} 

	if (iNode <= last_domain) {
		prepare_body_inOrderOf_domain((D+1)%3,p+npart[0],npart[1], domtree[iNode].son[1],send);
	}

}

void prepare_deliver_realloc_body(int *npart)
{
	int length = *npart;
	int n, new_npart;
	int  *sendcount = (int*)pmalloc(sizeof(int)*PROC_SIZE, 19);
	int  *recvcount = (int*)pmalloc(sizeof(int)*PROC_SIZE, 20);
	int  *senddisp  = (int*)pmalloc(sizeof(int)*PROC_SIZE, 21);
	int  *recvdisp  = (int*)pmalloc(sizeof(int)*PROC_SIZE, 22);

	for (n=0; n<PROC_SIZE; n++) {
		sendcount[n] = 0;
		recvcount[n] = 0;
		senddisp[n]  = 0;
		recvdisp[n]  = 0;
	}

	prepare_body_inOrderOf_domain(0, part, length, 0, sendcount);

	MPI_Alltoall(sendcount, 1, MPI_INT, recvcount, 1, MPI_INT, MPI_COMM_WORLD);

	new_npart = 0;

	for (n=0; n<PROC_SIZE; n++)
		new_npart += recvcount[n];

	recvdisp[0] = 0;
	senddisp[0] = 0;

	for (n=1; n<PROC_SIZE; n++) {
		recvdisp[n] = recvcount[n-1] + recvdisp[n-1];
		senddisp[n] = sendcount[n-1] + senddisp[n-1];
	}

	Body *part_b = (Body*)pmalloc(sizeof(Body)*new_npart, 1);
	//	if(PROC_RANK==0) printf("NPART_MEAN=%d,new_npart=%d\n",NPART_MEAN,new_npart);

#ifdef MYALLTOALLV
	int rank,nproc;
	rank=PROC_RANK;
	nproc=PROC_SIZE;
	MPI_Status status;
	MPI_Status *vstatus = (MPI_Status*)pmalloc(sizeof(MPI_Status)  * nproc, 57);
	MPI_Request *vreq  = (MPI_Request*) pmalloc(sizeof(MPI_Request) * nproc, 58);
	for (n=0; n<sendcount[rank]; n++) {
		*(part_b+recvdisp[rank] + n) = *(part+senddisp[rank] + n);
	}
	for (n=0;n<nproc;n++) {
		if (n != rank && sendcount[n] > 0) {
			MPI_Isend(part+senddisp[n], sendcount[n], strBody, n,n, MPI_COMM_WORLD, &vreq[n]);
		}
	}
	for (n=0;n<nproc;n++) {
		if (n != rank && recvcount[n] > 0)
			MPI_Recv(part_b+recvdisp[n], recvcount[n], strBody, n,rank, MPI_COMM_WORLD, &status);
	}	
	for (n=0;n<nproc;n++) {
		if (n != rank && sendcount[n] > 0)
			MPI_Wait(&vreq[n], &vstatus[n]);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	pfree(vstatus, 57);
	pfree(vreq, 58);
#else
	MPI_Alltoallv(part, sendcount, senddisp, strBody, part_b, recvcount, recvdisp, strBody, MPI_COMM_WORLD);
#endif

	pfree(part, 0);
	part = part_b;

	*npart = new_npart;
	mem_shift(1, 0);

	pfree(sendcount, 19);
	pfree(recvcount, 20);
	pfree(senddisp, 21);
	pfree(recvdisp, 22);

}






/////////////////////////////////////////////////////////////////////////////////

void domain_decomposition() 
{
	MPI_Barrier(MPI_COMM_WORLD);

	measure_domain_runtime(DTIME_FRACTION);

	determine_split_domtree(PROC_SIZE, NPART, domtree);


	prepare_deliver_realloc_body(&NPART);
}

/////////////////////////////////////////////////////////////////////////////////


void domain_volume_part(int iNode, double boxl[], double boxr[], int dim) {
	if (iNode >= first_domain)
		return;


	double bl[3], br[3];

	bl[0] = boxl[0];
	bl[1] = boxl[1];
	bl[2] = boxl[2];

	br[0] = boxr[0];
	br[1] = boxr[1];
	br[2] = boxr[2];

	double norm = domtree[iNode].time_left + domtree[iNode].time_right ;
	double frac   = bl[dim] +  (br[dim] - bl[dim]) * domtree[iNode].time_left / norm ;

	double brd = br[dim];
	br[dim] = frac;
	domtree[iNode].split = frac;


	domain_volume_part(domtree[iNode].son[0] , bl, br, (dim+1)%3) ;

	bl[dim] = frac;
	br[dim] = brd;

	domain_volume_part(domtree[iNode].son[1] , bl, br, (dim+1)%3) ;
}


void domain_initialize() {
	int n, length = 2*PROC_SIZE - 1;
	domtree = (DomainNode*)pmalloc(sizeof(DomainNode) * length, 24);

	for (n=0; n<first_domain; n++) {
		domtree[n].son[0] = 2*n+1;
		domtree[n].son[1] = 2*n+2;


		domtree[n].split = 0.0;
		domtree[n].split_previous = 0.0;

		domtree[n].time_node = 0.0;
		domtree[n].time_left = 0.0;
		domtree[n].time_right =0.0;

	}

	for (n=first_domain; n<=last_domain; n++) {
		domtree[n].son[0] = -1;
		domtree[n].son[1] = -1;

		domtree[n].time_node = 1.0;
		domtree[n].time_left = 1.0;
		domtree[n].time_right =1.0;
	}

	DTIME_FRACTION = 1.0;

	fill_time_domtree(0);

	double bl[3] = {    0.0, 0.0, 0.0};
	double br[3] = {BOXSIZE, BOXSIZE, BOXSIZE};
	int D = 0;

	domain_volume_part(0, bl, br, 0) ;
}


void domain_finalize() {
	pfree(domtree, 24);
}

