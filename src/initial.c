#include "photoNs.h"
#include "initial.h"
#include "snapshot.h"
#include "domains.h"
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <mpi.h>

/* snapshot parameters */

static double prmBOXSIZE = -1.0;
static double prmOMEGAM = -1.0;
static double prmOMEGAX = -1.0;
static double prmHUBBLE = -1.0;
static double prmSPLITSCALE  = -1.0;
static double prmINITIALTIME  = -1.0;
static double prmSOFTENSCALE  = -1.0;
static unsigned int prmNUMPART  = 0;

void ic_uniform(double box, int npart) ;

static void print_usage () {}

void readin_parameter_file(char parameterfile[])
{
	int n, sharp, space, first,second, length;
	char line[256];
	char pname[32];
	char value[256];

	FILE *fd;
	if ( !( fd = fopen(parameterfile, "r") ) ) {
		printf(" Cannot open parameterfile `%s'\n", parameterfile);
		printf(" Check the name or path of parameter file! \n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}



	while(1) {
		memset(line,0,sizeof(line));

		if ( fgets(line, 256, fd) != NULL ) {
			sharp = 0;
			for (n=0; n<256; n++) {
				if ('#' == line[n] && 0 == sharp  )
					sharp = 1;
				if (1 == sharp )
					line[n] = 0;
			}

			memset(pname,0,sizeof(pname));
			memset(value,0,sizeof(value));
			sscanf(line, "%s %s", pname, value);

			if ( 0 == (int)pname[0])
				continue;
			else if ( 0 == (int)value[0] ) {
				printf("wrong value of parameter '%s' !\n", pname);
				print_usage(0);
				exit(0);
			}



			if ( 0 == strcmp(pname, "OutputPath") ) {
				sprintf(OutputPath, "%s", value);
			}
			else if ( 0 == strcmp(pname, "OutputName") ) {
				sprintf(OutputName, "%s", value);
			}
			else if ( 0 == strcmp(pname, "InputPath") ) {
				sprintf(InputPath, "%s", value);
			}
			else if ( 0 == strcmp(pname, "CodeProj") ) {
				sprintf(CodeProj, "%s", value);
			}
			else if ( 0 == strcmp(pname, "OPENANGLE") ) {
				open_angle = atof(value);
				if (open_angle < 0.0) {
					printf(" error ! check OPENANGLE \n");
					exit(0);
				}
			}
			else if ( 0 == strcmp(pname, "SAMPLINGRATE") ) {
				SamplingRate = atof(value);
				if (SamplingRate < 0.0) {
					printf(" error ! check SAMPLINGRATE \n");
					exit(0);
				}
			}
			else if ( 0 == strcmp(pname, "BOXSIZE") ) {
				prmBOXSIZE = atof(value);
				if (prmBOXSIZE < 0.0) {
					printf(" error ! check BOXSIZE \n");
					exit(0);
				}
			}
			else if ( 0 == strcmp(pname, "OMEGAM") ) {
				prmOMEGAM = atof(value);
				if (prmOMEGAM < 0.0) {
					printf(" error ! check OMEGAM\n");
					exit(0);
				}
			}
			else if ( 0 == strcmp(pname, "OMEGAX") ) {
				prmOMEGAX = atof(value);
				if (prmOMEGAX < 0.0) {
					printf(" error ! check OMEGAX \n");
					exit(0);
				}
			}
			else if ( 0 == strcmp(pname, "HUBBLE") ) {
				prmHUBBLE = atof(value);
				if (prmHUBBLE < 0.0) {
					printf(" error ! check HUBBLE \n");
					exit(0);
				}
			}
			else if ( 0 == strcmp(pname, "SOFTENING") ) {
				prmSOFTENSCALE = atof(value);
				if (prmSOFTENSCALE < 0.0) {
					printf(" error ! check SOFTENING \n");
					exit(0);
				}
			}
			else if ( 0 == strcmp(pname, "SPLITSCALE") ) {
				prmSPLITSCALE = atof(value);
				if (prmSPLITSCALE < 0.0) {
					printf(" error ! check SPLITSCALE \n");
					exit(0);
				}
			}
			else if ( 0 == strcmp(pname, "GravConstant") ) {
				GravConst = atof(value);
			}
			else if ( 0 == strcmp(pname, "INITIALTIME") ) {
				prmINITIALTIME = atof(value);
			}
			else if ( 0 == strcmp(pname, "NUMPART") ) {
				prmNUMPART = atoi(value);
			}
			else if ( 0 == strcmp(pname, "SnapTime") ) {
				snaptime[number_snaptime] = atof(value);
				if (number_snaptime < MAXSNAP ) {
					number_snaptime ++;
				}
				else {
					printf(" error ! too many snaptime \n");
					exit(0);
				}
			}
			else if ( 0 == strcmp(pname, "SnapFormat") ) {
				SnapFormat = atoi(value);
			}	
			else if ( 0 == strcmp(pname, "NumMeshSide") ) {
				NSIDE = atoi(value);
			}
			else if ( 0 == strcmp(pname, "NumPartSide") ) {
				NPARTSIDE = atoi(value);
			}
			else if ( 0 == strcmp(pname, "NumThread") ) {
				NumThread = atoi(value);
			}
			else if ( 0 == strcmp(pname, "NprocVertical") ) {
				vproc[0] = atoi(value);
			}
			else if ( 0 == strcmp(pname, "NprocHorizon") ) {
				vproc[1] = atoi(value);
			}
			else if ( 0 == strcmp(pname, "MaxPackage") ) {
				MAXLEAF = atoi(value);
			}
			else if ( 0 == strcmp(pname, "SnapNumber") ) {
				SnapNumber = atoi(value);
			}
			else {
				printf("check the parameters %s, %s\n", pname, value);
				print_usage(0);
	
				exit(0);
			}

		}
		else {
			fclose(fd);
			break;
		}

	}
	if ( 1 == number_snaptime ) {
		printf(" miss SNAPTIME item !\n");
		exit(0);
	}

}

void setup_domain_index() {

	mostleft = 1;

	while ( mostleft < 2*PROC_SIZE-1 )
		mostleft *= 2;

	mostleft /=2;
	mostleft -= 1;

	if (PROC_SIZE == 1)
		mostleft = 0;


	this_domain = PROC_RANK + mostleft;

	if (this_domain > 2*PROC_SIZE - 2)
		this_domain -= PROC_SIZE;

	first_domain  = PROC_SIZE - 1;

	last_domain = 2*PROC_SIZE - 2;


}


#define decomp_2d_init get_local_size_

void initialize_nbody(char parameterfile[])
{
	int n, m;
	long long  n_start, n_end;

	MPI_Type_contiguous(sizeof(MKey), MPI_CHAR, &strMKey );
	MPI_Type_commit(&strMKey);

	MPI_Type_contiguous(sizeof(RemoteNode), MPI_CHAR, &strReNode );
	MPI_Type_commit(&strReNode);

	MPI_Type_contiguous(sizeof(RemoteBody), MPI_CHAR, &strReBody );
	MPI_Type_commit(&strReBody);

	MPI_Type_contiguous(sizeof(Body), MPI_CHAR, &strBody );
	MPI_Type_commit(&strBody);

	vproc[0] = PROC_SIZE;
	vproc[1] = 1; 

	////////////////
	number_snaptime = 1;

	NumThread = 1;
	open_angle = 0.3;
	SamplingRate = 0.3;

	InitialTime = 0.01;

	//////////////////


	readin_parameter_file(parameterfile);

	initializeLogfile();



	if (2 == SnapFormat) {
		if (SnapNumber == 1) {
			read_GadgetHeader(InputPath);
		} else {
			char fhead[200];
			sprintf(fhead,"%s.0", InputPath);		

			read_GadgetHeader(fhead);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	///////////////// FORCE PARAMETER //////////////////

	if (prmNUMPART > 0) {
		NPART = prmNUMPART;
		NPART_TOTAL = (long)(NPART * PROC_SIZE);
		if (0==PROC_RANK)
			printf("\n ATTENTION !  cast NPART = %u, NPART_TOTAL = %lu\n\n", NPART, NPART_TOTAL);          
	}
	if (prmINITIALTIME >= 0.0) {
		InitialTime = prmINITIALTIME;
		if (0==PROC_RANK)
			printf("\n ATTENTION !  cast Initial Time= %lf\n\n", InitialTime);          
	}
	if (prmOMEGAM >= 0.0) {
		OmegaM0 = prmOMEGAM ;
		if (0==PROC_RANK)
			printf("\n ATTENTION !  cast OmegaM0     = %lf\n\n", OmegaM0);          
	}
	if (prmOMEGAX >= 0.0) {
		OmegaX0 = prmOMEGAX ;
		if (0==PROC_RANK)
			printf("\n ATTENTION !  cast OmegaX0     = %lf\n\n", OmegaX0 );           
	}
	if (prmHUBBLE >= 0.0) {
		Hubble0 = prmHUBBLE ;
		if (0==PROC_RANK) 
			printf("\n ATTENTION !  cast Hubble0     = %lf\n\n", Hubble0);       
	}
	if (prmBOXSIZE >= 0.0) {
		BOXSIZE = prmBOXSIZE ;
		if (0==PROC_RANK) 
			printf("\n ATTENTION !  cast BOXSIZE     = %lf\n\n", BOXSIZE );   
	}
	///////////////// FORCE PARAMETER //////////////////


	snaptime[0] = 1.0/(1.0 + InitialTime);	

	double invside = BOXSIZE/((double)NSIDE);
	splitRadius = 1.25 * invside;
	SoftenScale = 0.03 * BOXSIZE/pow( ((double)NPART_TOTAL), 0.3333333 );
	if (-1 == SnapFormat) {
		SoftenScale =  0.03 * BOXSIZE/((double)NPARTSIDE);
	}


	if (0==PROC_RANK) { 
		printf(" DERIVED    + splitRadius = %lf \n", splitRadius );
		printf("            + SoftenScale = %lf \n\n", SoftenScale );
	}
	///////////////// FORCE PARAMETER //////////////////


	if (prmSPLITSCALE > 0.0) {
		splitRadius = prmSPLITSCALE;
		if (0==PROC_RANK) 
			printf("\n ATTENTION ! cast splitRadius = %lf \n\n", splitRadius );
	}

	cutoffRadius = 4.5*splitRadius;


	if (prmSOFTENSCALE >= 0.0) {
		SoftenScale = prmSOFTENSCALE;

		if (0==PROC_RANK) 
			printf("\n ATTENTION !  cast SoftenScale = %lf\n\n", SoftenScale );
	}

	nsample_rank = SamplingRate * NPART_MEAN;

	DTIME_FRACTION = 1.0;


	if (0 <  SnapFormat) {

		if (PROC_RANK == PROC_SIZE - 1)
			n_end = NPART_TOTAL;

		NPART  = n_end - n_start;

	}

	setup_domain_index();

	part  = (Body*)pmalloc(sizeof(Body)*NPART, 0);


	if (-2 == SnapFormat) {
		ic_uniform(BOXSIZE, NPART) ;
	}

	if (0 == SnapFormat) {
		if (SnapNumber == 1) {
			read_Particle_text(InputPath, n_start, NPART);
		} else {
			printf(" SnapNumber ! = 1 \n");
			exit(0);
		}
	}

	if (2 == SnapFormat) {
		if (SnapNumber == 1) {
			read_Particle_Gadget2(InputPath, n_start, NPART);
		} else {
			//        printf(" reading ... [%d]\n", PROC_RANK);
			int *np = (int*)pmalloc(sizeof(int)*SnapNumber, 3);
			int id, np_file[6];
			for (id=0; id<SnapNumber; id++) {
				npart_infile(InputPath, id, np_file);
				np[id] = np_file[1];
			}

			for (id=1; id<SnapNumber; id++) {
				np[id] = np[id-1] + np[id];
			}

			int file_start, file_end;

			for( id=0; id<SnapNumber; id++) {
				if(n_start<np[id])
					break;
			}
			file_start = id;

			for( id=SnapNumber-1; id>=0; id--) {
				if (n_start+NPART > np[id]) 
					break;				
			}
			file_end = id+1;

			char fname[200];
			int ip = 0;
			for ( id=file_start; id<=file_end; id++) {
				sprintf(fname,"%s.%d", InputPath, id);		
				int ns,ne;
				int low, high;
				if (id == 0 )
					low = 0 ;
				else
					low = np[id-1];

				high = np[id];

				if ( n_start > low )
					ns = n_start - low;
				else
					ns =  0;


				if ( n_start + NPART < high )
					ne = n_start + NPART - low;
				else
					ne = high - low;

				read_Particle_Gadget2_mfile(fname, ns, ne, ip );

				ip += ne - ns;
			}


			pfree(np, 3);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (NPART > 100000000) {
		printf(" NPART(%d) is too large over range of  int!\n", NPART);
		exit(0);
	}

	////////////////////////////////////////////////////////////////////////////////////

	pside[0] = NSIDE/vproc[0];
	pside[1] = NSIDE/vproc[1];

	if ( PROC_SIZE != vproc[0] * vproc[1] ) {
		if (0==PROC_RANK) {

			printf(" check number of proc !\n");
			printf(" vproc = %d(Y), %d(Z)\n", vproc[0], vproc[1]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		exit(0);
	}

	int nside[3] = { NSIDE, NSIDE, NSIDE };

	decomp_2d_init(nside, vproc, local_xstart, local_xend, local_xsize);



	MSIZE0 = (int*)malloc(sizeof(int) * PROC_SIZE);
	MSIZE1 = (int*)malloc(sizeof(int) * PROC_SIZE);

	MPI_Allgather( &(local_xend[1]), 1, MPI_INT, MSIZE0, 1, MPI_INT, MPI_COMM_WORLD  );
	MPI_Allgather( &(local_xend[2]), 1, MPI_INT, MSIZE1, 1, MPI_INT, MPI_COMM_WORLD  );


	double p_min[3] = {0.0, 0.0, 0.0};
	double p_max[3] = {BOXSIZE, BOXSIZE, BOXSIZE};
	double p_size[3]= {BOXSIZE, BOXSIZE, BOXSIZE};


	p_min[0] = (local_xstart[0])*invside;
	p_min[1] = (local_xstart[1])*invside;
	p_min[2] = (local_xstart[2])*invside;

	p_size[0] = (local_xsize[0])*invside;
	p_size[1] = (local_xsize[1])*invside;
	p_size[2] = (local_xsize[2])*invside;

	//    printf("[%d] %d %d %d\n", PROC_RANK, local_xsize[0],  local_xsize[1], local_xsize[2]);
	data_length = local_xsize[0]*local_xsize[1]*local_xsize[2];

	BoxMinimum = 0.0;
	BoxCenter  = 0.5*BOXSIZE;
	BoxMaximum = BOXSIZE;


	if (0==PROC_RANK) {
		printf("\n      #Proc = %d ( %d x %d )\n", PROC_SIZE, vproc[0], vproc[1]);
		printf("\n      NSIDE = %d x %d x %d\n", nside[0], nside[1], nside[2]);
		printf("\n      Theta = %lf\n", open_angle);
	}
	if (0==PROC_RANK)
		printf("\n\n      PROJECT    / %s /\n", CodeProj);



	if (0 == PROC_RANK) {
		printf("\n              initialize complete !\n\n");
	}

	return ;
}

void finalize_nbody()
{

	MPI_Type_free(&strMKey);
	MPI_Type_free(&strReNode);
	MPI_Type_free(&strReBody);
	MPI_Type_free(&strBody);

	pfree(part, 0); 

	finalizeLogfile() ;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////


void simple_snapshot(char fbase[]) {
	int n;
	FILE *fd;
	char fname[128];
	sprintf(fname, "%s.%d-%d", fbase, PROC_SIZE, PROC_RANK);

	if ( !(fd = fopen(fname,"w")) ) {
		printf("\n\nERROR! Cannot open/write snapshot file\n"); 
		printf("  - '%s' - (simple_snapshot)\n\n", fname);
		exit(0);
	}
	for (n=0; n<NPART; n++) {

	if (part[n].pos[2] < 0.1*BOXSIZE) {
		fprintf(fd, "%lf %lf %lf ", part[n].pos[0], part[n].pos[1], part[n].pos[2]);
//		fprintf(fd, "%lf %lf %lf ", part[n].vel[0], part[n].vel[1], part[n].vel[2]);
//		fprintf(fd, "%e %e %e ", part[n].acc[0], part[n].acc[1], part[n].acc[2]);
//		fprintf(fd, "%e %e %e", part[n].acc_pm[0], part[n].acc_pm[1], part[n].acc_pm[2]);
		fprintf(fd, "\n");
	}
	}
	fclose(fd);
}

void ic_uniform(double box, int npart) {
	int proc = this_domain;
	double xl[3], xr[3];
	xl[0] = xl[1] = xl[2] = 0;
	xr[0] = xr[1] = xr[2] = box;

	int left, right;
	int mid = PROC_SIZE;
	double frac;
	int d=0;

	proc = PROC_RANK;

	while (mid>1)
	{
		left = mid/2 + mid%2;
		right= mid - left;
		frac = (double)left/((double)mid);


		if (proc < left){        
			xr[d] = frac * (xr[d]-xl[d]) + xl[d];
			mid = left;
			//    printf(" [%d] left = %d \n", PROC_RANK, left);
		}
		else{
			xl[d] = frac * (xr[d]-xl[d]) + xl[d];
			proc-=left;
			mid = right;
			//    printf(" [%d] right= %d \n", PROC_RANK, right);
		}
		d = (d+1)%3;
	}


	long n;
	seed = 378412+PROC_RANK;

	double pmass=(OmegaM0*3*0.01)/(8.0*M_PI*GravConst)*(BOXSIZE*BOXSIZE*BOXSIZE)/(double)NPART_TOTAL;
	MASSPART = pmass;

	for (n=0; n<npart; n++) {
		part[n].pos[0] = ran3(&seed)*(xr[0] - xl[0]) + xl[0];
		part[n].pos[1] = ran3(&seed)*(xr[1] - xl[1]) + xl[1];
		part[n].pos[2] = ran3(&seed)*(xr[2] - xl[2]) + xl[2];

		part[n].vel[0] = 0.0;//10* (ran3(&seed) - 0.5);
		part[n].vel[1] = 0.0;//10* (ran3(&seed) - 0.5);
		part[n].vel[2] = 0.0;//10* (ran3(&seed) - 0.5);

		part[n].acc[0] = 0.0;
		part[n].acc[1] = 0.0;
		part[n].acc[2] = 0.0;


		part[n].acc_pm[0] = 0.0;
		part[n].acc_pm[1] = 0.0;
		part[n].acc_pm[2] = 0.0;
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////

double a_flat_lcdm_t(double time)
{
	double t_star = 3.0*sqrt(OmegaX0)/20.0;
	double kernel = sinh(t_star * time);
	double a = pow(kernel*kernel*OmegaM0/OmegaX0, 0.33333333f);
	return a;
}

double t_flat_lcdm_a(double a)
{
	double t_star = 3.0*sqrt(OmegaX0)/20.0;
	double a3 = a*a*a;
	double f = OmegaX0/OmegaM0;
	double time = log( sqrt(f*a3) + sqrt(1.0+f*a3) )/t_star;
	return time;
}

double kick_loga(double loga_i, double loga_f) {

	int n;
	int Nblock = 128;
	double kick_time = 0.0;
	double dloga = (loga_f - loga_i)/Nblock;
	double a_f = exp(loga_f);
	double a_i = exp(loga_i);
	double z1 = 1.0/(a_i);
	double h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
	kick_time = dloga*z1/h;
	for (n=1; n<Nblock; n++) {
		z1 = 1.0/(exp(loga_i+dloga*n));
		h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
		kick_time += 2.0*(1+n%2)*dloga*z1/h;
	}
	z1 = 1.0/(a_f);
	h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
	kick_time += dloga*z1/h;
	kick_time /= (3.0);
	return kick_time;
}

double drift_loga(double loga_i, double loga_f)
{
	int n;
	int Nblock = 128;
	double kick_time = 0.0;
	double dloga = (loga_f - loga_i)/Nblock;
	double a_f = exp(loga_f);
	double a_i = exp(loga_i);
	double z1 = 1.0/(a_i);
	double h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
	kick_time = dloga*z1*z1/h;
	for (n=1; n<Nblock; n++) {
		z1 = 1.0/(exp(loga_i+dloga*n));
		h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
		kick_time += 2.0*(1+n%2)*dloga*z1*z1/h;
	}
	z1 = 1.0/(a_f);
	h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
	kick_time += dloga*z1*z1/h;
	kick_time /= (3.0);
	return kick_time;
}


