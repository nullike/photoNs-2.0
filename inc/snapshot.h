#ifndef SNAPSHOT_H
#define SNAPSHOT_H

void read_GadgetHeader(const char filename[]);
void read_Particle_Gadget2(char filename[], int n_start, int n_count);
void read_Particle(char filename[], int n_start, int n_count);
void read_Particle_text(char filename[], int n_start, int n_count);
void read_Particle_Gadget2_mfile(char filename[], int n_start, int n_end, int ip);

void npart_infile(const char filename[], int ith, int np[]);
void write_snapshot(int type, int timestamp, int rank ); 
void write_Particle_Gadget2(char filename[], int n_start, int n_count);

#endif
