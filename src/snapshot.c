#include "photoNs.h"
#include <stdlib.h>
#include <stdio.h>

typedef struct
{
	int npart[6];
	double mass[6];
	double time;
	double redshift;
	int flag_sfr;
	int flag_feedback;
	int npartTotal[6];
	int flag_cooling;
	int num_files;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];	
	/* fills to 256 Bytes */
} GadgetHeader;


static double partmass[6];
static double Omega0;
static double OmegaLambda;
static double HubbleParam;


void npart_infile(const char filename[], int ith, int np[])
{
	int dummy, byte;
	FILE* fsnap;
	GadgetHeader head;

	char fname[200];
	sprintf(fname, "%s.%d",filename, ith);


	if ( !(fsnap=fopen(fname,"r")) )
	{
		printf("[read_GadgetHeader] cannot open file %s \n", filename);
		exit(0);
	}

	byte = fread(&dummy, sizeof(dummy),        1, fsnap);
	byte = fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte = fread(&dummy, sizeof(dummy),        1, fsnap);

	BOXSIZE = head.BoxSize;

	np[0] = head.npart[0];
	np[1] = head.npart[1];
	np[2] = head.npart[2];
	np[3] = head.npart[3];
	np[4] = head.npart[4];
	np[5] = head.npart[5];

	fclose(fsnap);

}

void read_GadgetHeader(const char filename[])
{
	int dummy, byte;
	FILE* fsnap;
	GadgetHeader head;

	if ( !(fsnap=fopen(filename,"r")) )
	{
		printf("[read_GadgetHeader] cannot open file %s \n", filename);
		exit(0);
	}

	byte = fread(&dummy, sizeof(dummy),        1, fsnap);
	byte = fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte = fread(&dummy, sizeof(dummy),        1, fsnap);

	BOXSIZE = head.BoxSize;

	NPART_TOTAL = head.npartTotal[0];
	NPART_TOTAL+= head.npartTotal[1];
	NPART_TOTAL+= head.npartTotal[2];
	NPART_TOTAL+= head.npartTotal[3];
	NPART_TOTAL+= head.npartTotal[4];
	NPART_TOTAL+= head.npartTotal[5];

	MASSPART  = head.mass[1];

	partmass[0] = head.mass[0];
	partmass[1] = head.mass[1];
	partmass[2] = head.mass[2];
	partmass[3] = head.mass[3];
	partmass[4] = head.mass[4];
	partmass[5] = head.mass[5];

	BOXSIZE = head.BoxSize;
	OmegaM0 = head.Omega0;
	OmegaX0 = head.OmegaLambda;
	Hubble0 = head.HubbleParam;
	InitialTime = head.redshift;

	if (0==PROC_RANK) {
		printf("\n PARAMETERS from IC snapshots : \n");
		printf("            + Omega Matter =  %lf\n", OmegaM0 );
		printf("            + Omega Lambda =  %lf\n", OmegaX0 );
		printf("            + Hubble Const =  %lf\n", Hubble0 );
		printf("            + Box Size = %lf \n", BOXSIZE);   
		printf("            + Num Particle =  %ld\n", NPART_TOTAL);   
		printf("            + IC Red-Shift =  %lf\n", InitialTime);   
		printf("\n");   
		fflush(stdout);
	}
	dummy = byte;

	fclose(fsnap);

}

void read_Particle_Gadget2_mfile(char filename[], int n_start, int n_end, int ip)
{
	int n, m, c, p;
	int n_count;
	int dummy, byte, ID;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	int ID_SIZE = sizeof(int);

	float data[3];

	n_count = n_start - n_end;

	FILE* fsnap;
	GadgetHeader head;
	if ( !(fsnap=fopen(filename,"r")) ) {
		printf("[read_Particle] cannot open mfile %s \n", filename);
		exit(0);
	}

	printf(" < <  reading `%s'\n", filename);

	byte = 0;
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);
	byte += fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);


	Omega0 = head.Omega0;
	OmegaLambda = head.OmegaLambda;
	HubbleParam = head.HubbleParam;

	//    printf(" npart = %d %d %d \n", head.npart[0], head.npart[1], head.npart[2]);
	//    printf(" mass = %lf %lf %lf \n", head.mass[0], head.mass[1], head.mass[2]);
	//    printf(" time = %lf , redshift = %lf\n", head.time, head.redshift);
	//    printf(" OmegaM = %lf, %lf, Hubble= %lf \n", head.Omega0, head.OmegaLambda, head.HubbleParam);


	byte += fread(&dummy, sizeof(int), 1, fsnap);   
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(data, DATA_SIZE, 3, fsnap);

			if (c>=n_start && c<n_end) {
				part[ip+p].pos[0] = (double)data[0];
				part[ip+p].pos[1] = (double)data[1];
				part[ip+p].pos[2] = (double)data[2];

				p++;
			}
			c++;
		}
	}
	byte += fread(&dummy, sizeof(int), 1, fsnap);

	double gdt2unit = pow( 1.0/(1.0+head.redshift), 1.5); //for gadget2 format
	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(data, DATA_SIZE, 3, fsnap);

			if (c>=n_start && c<n_end) {
				part[ip+p].vel[0] = (double)data[0] * gdt2unit;
				part[ip+p].vel[1] = (double)data[1] * gdt2unit;
				part[ip+p].vel[2] = (double)data[2] * gdt2unit;
				p++;
			}
			c++;
		}
	}
	byte += fread(&dummy, sizeof(int), 1, fsnap);
	/*
	   byte += fread(&dummy, sizeof(int), 1, fsnap);
	   for (c=0, p=0, m=0; m<6; m++) {
	   for (n=0; n<head.npart[m]; n++) {
	   byte += fread(&ID, ID_SIZE, 1, fsnap);

	   if (c>=n_start && c<n_end) {
	   part[ip+p].id = (long)ID;
	   p++;
	   }
	   c++;
	   }
	   }
	   byte += fread(&dummy, sizeof(int), 1, fsnap);
	   */
	fclose(fsnap);
}

void read_Particle_Gadget2(char filename[], int n_start, int n_count)
{
	int n, m, c, p;
	int n_end;
	int dummy, byte, ID;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	int ID_SIZE = sizeof(int);

	float data[3];

	n_end = n_start + n_count;

	FILE* fsnap;
	GadgetHeader head;
	if ( !(fsnap=fopen(filename,"r")) ) {
		printf("[read_Particle_Gadget2] cannot open ic %s \n", filename);
		exit(0);
	}

	printf(" > reading `%s'\n", filename);

	byte = 0;
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);
	byte += fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);

	Omega0 = head.Omega0;
	OmegaLambda = head.OmegaLambda;
	HubbleParam = head.HubbleParam;

	byte += fread(&dummy, sizeof(int), 1, fsnap);   
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(data, DATA_SIZE, 3, fsnap);

			if (c>=n_start && c<n_end) {
				part[p].pos[0] = (double)data[0];
				part[p].pos[1] = (double)data[1];
				part[p].pos[2] = (double)data[2];
				//part[p].mass   =  head.mass[m];

				p++;
			}
			c++;
		}
	}
	byte += fread(&dummy, sizeof(int), 1, fsnap);

	double gdt2unit = pow( 1.0/(1.0+head.redshift), 1.5); //for gadget2 format
	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(data, DATA_SIZE, 3, fsnap);

			if (c>=n_start && c<n_end) {
				part[p].vel[0] = (double)data[0] * gdt2unit;
				part[p].vel[1] = (double)data[1] * gdt2unit;
				part[p].vel[2] = (double)data[2] * gdt2unit;
				p++;
			}
			c++;
		}
	}
	byte += fread(&dummy, sizeof(int), 1, fsnap);
	/*
	   byte += fread(&dummy, sizeof(int), 1, fsnap);
	   for (c=0, p=0, m=0; m<6; m++) {
	   for (n=0; n<head.npart[m]; n++) {
	   byte += fread(&ID, ID_SIZE, 1, fsnap);

	   if (c>=n_start && c<n_end) {
	   part[p].id = (long)ID;
	   p++;
	   }
	   c++;
	   }
	   }
	   byte += fread(&dummy, sizeof(int), 1, fsnap);
	   */
	fclose(fsnap);
}


void read_Particle(char filename[], int n_start, int n_count)
{
	int n;
	int dummy, byte;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	float data[3];

	FILE* fsnap;
	GadgetHeader head;
	if ( !(fsnap=fopen(filename,"r")) ) {
		printf("[read_Particle] cannot open file %s \n", filename);
		exit(0);
	}

	printf(" > reading file `%s' \n", filename);

	byte = 0;
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);
	byte += fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);

	num_part = NPART_TOTAL;

	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (n=0; n<n_start; n++)
		byte += fread(data, DATA_SIZE, 3, fsnap);

	for (n=0; n<n_count; n++) {
		byte += fread(data, DATA_SIZE, 3, fsnap);
		part[n].pos[0] = (double)data[0];
		part[n].pos[1] = (double)data[1];
		part[n].pos[2] = (double)data[2];
	}

	//    printf(" [%d] pos byte = %d\n", PROC_RANK, byte);
	for (n=0; n<num_part-n_count-n_start; n++)
		byte += fread(data, DATA_SIZE , 3, fsnap);

	byte += fread(&dummy, sizeof(int), 1, fsnap);

	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (n=0; n<n_start; n++)
		byte += fread(data, DATA_SIZE, 3, fsnap);

	for (n=0; n<n_count; n++) {
		byte += fread(&data[0], DATA_SIZE, 3, fsnap);
		part[n].vel[0] = (double)data[0];
		part[n].vel[1] = (double)data[1];
		part[n].vel[2] = (double)data[2];
	}

	/*
	   disp   = (2*sizeof(dummy) + sizeof(GadgetHeader));
	   disp  += (2*sizeof(dummy) + 2*num_part*3*DATA_SIZE  );
	   disp  += (2*sizeof(dummy) + n_start*sizeof(int) );
	   fseek(fsnap, disp,0);
	   for (n=0; n<n_count; n++)
	   {
	   fread(&id, sizeof(int), 1, fsnap);
	   part[n].ID = (long)id;
	   }
	   */

	fclose(fsnap);

}

void read_Particle_text(char filename[], int n_start, int n_count) 
{

	int n, cnt;
	FILE *fsnap;
	if ( !(fsnap=fopen(filename,"r")) ) {
		printf("[read_Particle] cannot open file %s \n", filename);
		exit(0);
	}

	float data[6];

	for (n=0; n<n_start; n++)
		cnt=fscanf(fsnap, "%f %f %f %f %f %f", data, data+1, data+2, data+3, data+4, data+5);

	for (n=0; n<n_count; n++) {

		cnt=fscanf(fsnap, "%f %f %f %f %f %f", data, data+1, data+2, data+3, data+4, data+5);
		part[n].pos[0] = data[0];
		part[n].pos[1] = data[1];
		part[n].pos[2] = data[2];

		part[n].vel[0] = data[3];
		part[n].vel[1] = data[4];
		part[n].vel[2] = data[5];

	}

	fclose(fsnap);
}


void write_Particle_Gadget2(char filename[], int n_start, int n_count)
{
	int n, m, c, p;
	int n_end;
	int dummy, byte, ID;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	int ID_SIZE = sizeof(int);

	float data[3];

	n_end = n_start + n_count;

	FILE* fsnap;
	GadgetHeader head;
	if ( !(fsnap=fopen(filename,"w")) ) {
		printf("[read_Particle] cannot open file %s \n", filename);
		exit(0);
	}

	printf(" > writing `%s' \n", filename);

	head.num_files = 1;
	head.BoxSize = BOXSIZE;

	head.Omega0 = Omega0 ;
	head.OmegaLambda = OmegaLambda ;
	head.HubbleParam = HubbleParam ;

	head.npart[0] = 0; 
	head.npart[1] = n_count;
	head.npart[2] = 0; 
	head.npart[3] = 0; 
	head.npart[4] = 0; 
	head.npart[5] = 0; 

	head.mass[0]  = partmass[0];
	head.mass[1]  = partmass[1];
	head.mass[2]  = partmass[2];
	head.mass[3]  = partmass[3];
	head.mass[4]  = partmass[4];
	head.mass[5]  = partmass[5];


	head.npartTotal[0] = 0; 
	head.npartTotal[1] = NPART_TOTAL; 
	head.npartTotal[2] = 0; 
	head.npartTotal[3] = 0; 
	head.npartTotal[4] = 0; 
	head.npartTotal[5] = 0; 
	head.time = 1.0/(Redshift_Time+1);
	head.redshift = Redshift_Time;

	byte = 0;
	byte += fwrite(&dummy, sizeof(dummy),        1, fsnap);
	byte += fwrite(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte += fwrite(&dummy, sizeof(dummy),        1, fsnap);



	byte += fwrite(&dummy, sizeof(int), 1, fsnap);   
	for (c=n_start; c<n_end; c++) {

		data[0] = (float)part[c].pos[0];
		data[1] = (float)part[c].pos[1];
		data[2] = (float)part[c].pos[2];

		byte += fwrite(data, DATA_SIZE, 3, fsnap);
	}
	byte += fwrite(&dummy, sizeof(int), 1, fsnap);

	double gdt2unit = pow( 1.0/(1.0+head.redshift), 1.5); 
	byte += fwrite(&dummy, sizeof(int), 1, fsnap);   
	for (c=n_start; c<n_end; c++) {

		data[0] = (float)part[c].vel[0] / gdt2unit;
		data[1] = (float)part[c].vel[1] / gdt2unit;
		data[2] = (float)part[c].vel[2] / gdt2unit;

		byte += fwrite(data, DATA_SIZE, 3, fsnap);
	}
	byte += fwrite(&dummy, sizeof(int), 1, fsnap);   
	/*
	   byte += fwrite(&dummy, sizeof(int), 1, fsnap);   
	   for (c=n_start; c<n_end; c++) {
	   ID = part[c].id;
	   byte += fwrite(&ID, ID_SIZE, 1, fsnap);
	   }
	   byte += fwrite(&dummy, sizeof(int), 1, fsnap);   
	   */
	fclose(fsnap);

	MPI_Barrier(MPI_COMM_WORLD);

/*
	if (0 == PROC_RANK) {
		printf(" [%d] head byte = %d\n", PROC_RANK, byte);

		printf(" npart = %d %d %d \n", head.npart[0], head.npart[1], head.npart[2]);
		printf(" mass = %lf %lf %lf \n", head.mass[0], head.mass[1], head.mass[2]);
		printf(" > Scale factor = %lf , redshift = %lf\n", head.time, head.redshift);
		printf(" OmegaM = %lf, %lf, Hubble= %lf \n", head.Omega0, head.OmegaLambda, head.HubbleParam);
	}
*/

}

void write_Particle_text(char filename[], int n_start, int n_count)
{
	int n, m, c, p;
	int n_end;
	int dummy, byte, ID;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	int ID_SIZE = sizeof(int);

	float data[3];

	n_end = n_start + n_count;

	FILE* fsnap;
	if ( !(fsnap=fopen(filename,"w")) ) {
		printf("[read_Particle] cannot open file %s \n", filename);
		exit(0);
	}

	printf(" > writing `%s' \n", filename);

	for (c=n_start; c<n_end; c++) {
		fprintf(fsnap, "%lf %lf %lf\n", part[c].pos[0], part[c].pos[1], part[c].pos[2]);
	}

	fclose(fsnap);

	MPI_Barrier(MPI_COMM_WORLD);

}


void write_snapshot(int type, int timestamp, int rank ) 
{
	char fout[64];

	sprintf(fout, "%s/%s_%d.%d", OutputPath, OutputName, timestamp, rank);

	//		write_Particle_text(fout, 0, NPART);

	if (type == 2)
		write_Particle_Gadget2(fout, 0, NPART);

}


