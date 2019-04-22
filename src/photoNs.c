#include "photoNs.h"
#include "toptree.h"
#include "remotes.h"
#include "domains.h"
#include "fmm.h"
#include "snapshot.h"
#include "partmesh.h"
#include "initial.h"
#include <mpi.h>
#include <pthread.h>

void make_title(){
	if (0 == PROC_RANK) {
		printf(" \n\n\n");	
		printf("                    / photoNs 2 /\n\n");

	}
}

void make_appendix(double total_time) {
	if (0 == PROC_RANK) {
		printf(" \n\n");	
		printf("                      complete !\n\n");
		printf(" Elapsed (wall-clock) Time %.1lf[sec]\n", total_time);
		printf(" Project: %s\n", CodeProj);
		printf(" output: %s\n", OutputPath);
		printf(" \n\n");	
	}
}

double total_time_start; 
double total_time;


void driver(double ai, double af, int snap_idx, int nstep_fix) {
	int n;
	double memproc;  
	double time0, time1, DTIME_DOMAIN;
	double DTIME_TOTAL_DOMAIN ;
	double DTIME_MAXIMUM;

	double time_loop;
	double maxmem;  
	double maxgb;

	int loop;
	int Nstep = nstep_fix;


	double dKick_half;
	size_t mem;
	double dloga =  (log(af) - log(ai) ) /Nstep;
	double loga_i, loga_f, dk, dkh, dd;


	double imbalance ;


	for (n=0; n<NPART; n++) {
		part[n].acc[0] = 0.0;
		part[n].acc[1] = 0.0;
		part[n].acc[2] = 0.0;

		part[n].acc_pm[0] = 0.0;
		part[n].acc_pm[1] = 0.0;
		part[n].acc_pm[2] = 0.0;
	}
	DTIME_FRACTION = 1.0;

#ifdef PMONLY
	DTIME_THIS_DOMAIN = 1.0;
#endif
	reset_mem();

	//  use_sample_estimate();
	// use_time_estimate();

	MPI_Barrier(MPI_COMM_WORLD);
	if (0==PROC_RANK) {
		printf("begin domain_initialize\n");
	}

	domain_initialize();


	if (0==PROC_RANK) {
		printf("begin domain_decomposition\n");
	}
	domain_decomposition() ;


	MPI_Barrier(MPI_COMM_WORLD);


#ifndef PMONLY

	fmm_construct() ;


	fmm_prepare();

#ifdef PMTHREAD
		pthread_t tpm;
		pthread_create(&tpm, NULL, pm_thread, NULL);
#endif


	fmm_task();

#ifdef PMTHREAD
		pthread_join(tpm, NULL);
#else
		partmesh();
#endif

	fmm_ext();

#endif


	if (0==PROC_RANK) {
		printf("\n dtime_pm = %lf\n dtime_m2l = %lf, dtime_p2p = %lf [sec]\n", dtime_pm, dtime_m2l, dtime_p2p);
		printf(" dtime_fmm/task = %lf, task = %lf [sec]\n", dtime_fmm, dtime_task);
		if (dtime_pm > dtime_task)
			printf("   dTime PM > dTime Task\n");

		printf(" dtime_prep = %lf, ext = %lf [sec]\n", dtime_prep, dtime_ext); 

	}


	MPI_Barrier(MPI_COMM_WORLD);

	imbalance =  0.0;

	for (loop = 0; loop < Nstep; loop++) 
	{
		time_loop = dtime();
		loop_step ++;

		loga_i = loop * dloga + log( ai);
		loga_f = (loop+1) * dloga + log( ai);

		dk = kick_loga(loga_i, loga_f);
		dd = drift_loga(loga_i, loga_f);

		if (PROC_RANK  == 0)
			printf("\n\nLOOP         a=( %lf to %lf )        %5d\n\n",  exp(loga_i), exp(loga_f), loop_step);

		dkh = 0.5*dk*GravConst;

		dKick_half = dkh;
		pack2pack_count = 0;
		pack2pack_inleaf_count = 0;



		for (n=0; n<NPART; n++) {
			part[n].vel[0] += part[n].acc_pm[0]*dkh;
			part[n].vel[1] += part[n].acc_pm[1]*dkh;
			part[n].vel[2] += part[n].acc_pm[2]*dkh;
		}


		for (n=0; n<NPART; n++) {
			part[n].vel[0] += part[n].acc[0]*dkh;
			part[n].vel[1] += part[n].acc[1]*dkh;
			part[n].vel[2] += part[n].acc[2]*dkh;
		}

		for (n=0; n<NPART; n++) {
			part[n].pos[0] += part[n].vel[0]*dd;
			part[n].pos[1] += part[n].vel[1]*dd;
			part[n].pos[2] += part[n].vel[2]*dd;
		}


		for (n=0; n<NPART; n++) {
			while (part[n].pos[0]  < 0.0 )
				part[n].pos[0] += BOXSIZE;

			while (part[n].pos[1]  < 0.0 )
				part[n].pos[1] += BOXSIZE;

			while (part[n].pos[2]  < 0.0 )
				part[n].pos[2] += BOXSIZE;

			while (part[n].pos[0]  >= BOXSIZE )
				part[n].pos[0] -= BOXSIZE;

			while (part[n].pos[1]  >= BOXSIZE )
				part[n].pos[1] -= BOXSIZE;

			while (part[n].pos[2]  >= BOXSIZE )
				part[n].pos[2] -= BOXSIZE;
		}

#ifndef PMONLY
		fmm_deconstruct() ;
#endif

		time0 = dtime();

		//use_time_estimate();
		domain_decomposition() ;


		DTIME_DOMAIN = dtime() - time0;

		if (0 ==PROC_RANK)
			printf(" * DOM reconstruct     (dTime = %lf sec)\n", DTIME_DOMAIN);

	

		for (n=0; n<NPART; n++) {
			part[n].acc_pm[0] = 0.0;
			part[n].acc_pm[1] = 0.0;
			part[n].acc_pm[2] = 0.0;
		}


		for (n=0; n<NPART; n++) {
			part[n].acc[0] = 0.0;
			part[n].acc[1] = 0.0;
			part[n].acc[2] = 0.0;
		}	


#ifndef PMONLY

		fmm_construct() ;

		fmm_prepare();


#ifdef PMTHREAD
		pthread_create(&tpm, NULL, pm_thread, NULL);
#endif


		fmm_task();


#ifdef PMTHREAD
		pthread_join(tpm, NULL);
#else
	
		partmesh();

#endif

		fmm_ext();

#endif


		for (n=0; n<NPART; n++) {
			part[n].vel[0] += part[n].acc[0]*dkh;
			part[n].vel[1] += part[n].acc[1]*dkh;
			part[n].vel[2] += part[n].acc[2]*dkh;
		}


		for (n=0; n<NPART; n++) {
			part[n].vel[0] += part[n].acc_pm[0]*dkh;
			part[n].vel[1] += part[n].acc_pm[1]*dkh;
			part[n].vel[2] += part[n].acc_pm[2]*dkh;
		}


#ifdef PMONLY
		DTIME_THIS_DOMAIN = 1.0;
#endif
		DTIME_TOTAL_DOMAIN = 0.0;
		DTIME_MAXIMUM= 0.0;

		MPI_Allreduce( &DTIME_THIS_DOMAIN, &DTIME_TOTAL_DOMAIN, 
				1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

		MPI_Allreduce( &DTIME_THIS_DOMAIN, &DTIME_MAXIMUM, 
				1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);

		DTIME_FRACTION = DTIME_THIS_DOMAIN*PROC_SIZE/(DTIME_TOTAL_DOMAIN+0.0001);


		imbalance =  1.0 - DTIME_TOTAL_DOMAIN /((double) PROC_SIZE * DTIME_MAXIMUM ); 

		if (0 ==PROC_RANK)
			printf("\n work-load imbalance = %1.3f %%\n", imbalance*100.0);


		total_time = dtime() - total_time_start;


		Logfile_flush();

		if (0==PROC_RANK) {
			printf("\n dtime_pm = %lf\n dtime_m2l = %lf, dtime_p2p = %lf [sec]\n", dtime_pm, dtime_m2l, dtime_p2p);
			printf(" dtime_fmm/task = %lf, task = %lf [sec]\n", dtime_fmm, dtime_task);
			if (dtime_pm > dtime_task)
				printf("    dTime PM > dTime Task\n");

			printf(" dtime_prep = %lf, ext = %lf [sec]\n", dtime_prep, dtime_ext); 

		}

		time_loop = dtime() - time_loop;
		if (0==PROC_RANK)
			printf(" dtime_loop = %lf, total_time= %lf[sec]\n", time_loop, total_time ); 

		LogMessage(loop_step, 0.5*(exp(loga_f)+exp(loga_i)), time_loop, dtime_pm, total_time, imbalance ) ;

	}
	//////////////////////// LOOP //////////////////////// 

#ifndef PMONLY
	fmm_deconstruct() ;
#endif

	domain_finalize();

	Redshift_Time = 1.0/af - 1.0;

	//	simple_snapshot("snap");
	write_snapshot(SnapFormat, snap_idx, PROC_RANK );

}

int main(int argc, char* argv[]) 
{

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &PROC_SIZE);
	MPI_Comm_rank(MPI_COMM_WORLD, &PROC_RANK);

#ifdef PMTHREAD
	MPI_Comm_dup(MPI_COMM_WORLD, &PM_COMM_WORLD);
#endif

	make_title();

	total_time_start = dtime(); 

	initialize_nbody(argv[1]);

	loop_step = 0;
	count_crushing_step = 0;

	driver(1.0/(1.0+InitialTime), 1.0, 3, 512);

	finalize_nbody();

	total_time = dtime() - total_time_start;
	make_appendix(total_time);


	MPI_Finalize();
	return 0;
}


