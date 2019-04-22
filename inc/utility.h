#ifndef UTILITY_H
#define UTILITY_H





#define MEMBLOCK 100
size_t memused[MEMBLOCK];
size_t TotalMemory;
size_t MaxMemory;



size_t total_mem_used();

void reset_mem() ;
void* pmalloc (size_t size, int idx) ;
void mem_shift(int source, int target);
void pfree(void* ptr, int idx);


double  dtime();
float ran3(long *idum);
void* xmalloc (size_t size);

void LogMessage(int loop,double a,double time_fmm,double time_pm,double time_total,double inbalance);
void initializeLogfile();
void Logfile_flush();

void finalizeLogfile();

#endif
