#include "photoNs.h"
#include <sys/time.h>
#include <stdio.h>
#include <string.h>

static FILE *flog;

void LogMessage(int loop, double a, double time_short, double time_pm, double time_total, double imbalance) {
	if (0== PROC_RANK) {
		if (0 == flog) {
			printf(" error log files!\n");
			exit(0);
		}

		fprintf(flog, "%5d %3d  %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", loop_step, adaptive_level_maximum,  imbalance, a, time_pm, time_short, dtime_p2p, dtime_m2l, dtime_fmm, dtime_fmm_remote , time_total );
		//    fprintf(flog, "\n");
	}
}

void Logfile_flush(){
	if (0== PROC_RANK) {
		fflush(flog);
	}
}

void initializeLogfile() {
	if (0== PROC_RANK) {
		char fname[256];
		sprintf(fname, "%s/LOG%s.TXT", OutputPath, CodeProj);

		flog = fopen(fname, "w");
		if (0 == flog) {
			printf(" error log files!\n");
			exit(0);
		}
		fprintf(flog, "### n lvl  imbalance a_t      dTpm      dTshort  dTp2p    dTm2l  dTfmm    dText     Ttot \n");

	}
}



void finalizeLogfile() {

	if (0== PROC_RANK) {
		if (0 != flog) {
			fclose(flog);
		}
	}
}

void reset_mem() {
	int n;
	for (n=0; n<MEMBLOCK; n++) {
        if (n!= 6)
		memused[n] = 0;
	}

}

size_t total_mem_used() {
	int n;
	size_t total = 0;
	for (n=0; n<MEMBLOCK; n++)
		total += memused[n];	

	TotalMemory = total;

	if (total > MaxMemory)
		MaxMemory = total;

	return total;
} 

void* pmalloc (size_t size, int idx) {
	void *alloc = malloc(size);

	if (alloc == 0) {
		printf("[%d]out of memory (malloc) idx=%d\n", PROC_RANK, idx);
		exit(1000);
	}
	memset(alloc, 0, size);

	if (idx >= MEMBLOCK) {
		printf("uncorrect idx in xmallox!\n");
		exit(0);
	}

	memused[idx] = size ;
	TotalMemory += memused[idx];

	if (TotalMemory > MaxMemory)
		MaxMemory = TotalMemory;

	return alloc;
}

void pfree(void* ptr, int idx) {

	if (idx >= MEMBLOCK) {
		printf("uncorrect idx in xmallox!\n");
		exit(0);
	}

	free(ptr);

	TotalMemory -= memused[idx];

	memused[idx] = 0 ;
}

void mem_shift(int source, int target){
	memused[target] = memused[source];
}


double dtime() // Return double style timestamp using gettimeofday().
{
	double tseconds;
	struct timeval mytime;

	gettimeofday(&mytime, NULL);
	tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);

	return (tseconds);

} /* dtime() */


////////////////// NR ran3 ///////////////

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1; i<=54; i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1; k<=4; k++)
			for (i=1; i<=55; i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


