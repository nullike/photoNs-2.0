/*
 * photoNs-2
 *
 *      2017 - 8 - 21
 *	qwang@nao.cas.cn
 */	 

#include "photoNs.h"
#include "operator.h"
#include <stdio.h>
#include <stdlib.h>

void p2m(int iPart, int npart, double center[3], double M[]) 
{
	int n, idx;
	double dx[3];
	double m0, m1, m2, m3, m4, m5;
	for (n=0; n<NMULTI; n++) 
		M[n] = 0.0;


	for (idx=iPart; idx<npart+iPart;  idx++)
	{
		dx[0] = part[idx].pos[0] - center[0];
		dx[1] = part[idx].pos[1] - center[1];	
		dx[2] = part[idx].pos[2] - center[2];	

		m0 = MASSPART;
		m1 = -m0;

		M[0] += m0;

		M[X] += m1*dx[0];
		M[Y] += m1*dx[1];
		M[Z] += m1*dx[2];

#ifdef QUADRUPOLE
		m2 = m0;

		// 2nd order

		M[XX] += m2*dx[0]*dx[0]/2;
		M[XY] += m2*dx[0]*dx[1];
		M[XZ] += m2*dx[0]*dx[2];
		M[YY] += m2*dx[1]*dx[1]/2;
		M[YZ] += m2*dx[1]*dx[2];
		M[ZZ] += m2*dx[2]*dx[2]/2;

#endif

#ifdef  OCTUPOLE
		m3 = -m0;

		// 3rd order
		M[XXX] += m3 * dx[0] * dx[0] * dx[0]/6;
		M[XXY] += m3 * dx[1] * dx[0] * dx[0]/2;
		M[XXZ] += m3 * dx[2] * dx[0] * dx[0]/2;
		M[XYY] += m3 * dx[1] * dx[1] * dx[0]/2;
		M[XYZ] += m3 * dx[2] * dx[1] * dx[0];
		M[XZZ] += m3 * dx[2] * dx[2] * dx[0]/2;
		M[YYY] += m3 * dx[1] * dx[1] * dx[1]/6;
		M[YYZ] += m3 * dx[1] * dx[2] * dx[1]/2;
		M[YZZ] += m3 * dx[2] * dx[2] * dx[1]/2;
		M[ZZZ] += m3 * dx[2] * dx[2] * dx[2]/6;	

#endif

#ifdef  HEXADECAPOLE

		m4 = m0;  

		M[XXXX] += m4 * dx[0] * dx[0] * dx[0] * dx[0]/24;
		M[XXXY] += m4 * dx[0] * dx[0] * dx[0] * dx[1]/6;
		M[XXXZ] += m4 * dx[0] * dx[0] * dx[0] * dx[2]/6;
		M[XXYY] += m4 * dx[0] * dx[0] * dx[1] * dx[1]/4;
		M[XXYZ] += m4 * dx[0] * dx[0] * dx[1] * dx[2]/2;
		M[XXZZ] += m4 * dx[0] * dx[0] * dx[2] * dx[2]/4;
		M[XYYY] += m4 * dx[0] * dx[1] * dx[1] * dx[1]/6;
		M[XYYZ] += m4 * dx[0] * dx[1] * dx[1] * dx[2]/2;
		M[XYZZ] += m4 * dx[0] * dx[1] * dx[2] * dx[2]/2;
		M[XZZZ] += m4 * dx[0] * dx[2] * dx[2] * dx[2]/6;
		M[YYYY] += m4 * dx[1] * dx[1] * dx[1] * dx[1]/24;
		M[YYYZ] += m4 * dx[1] * dx[1] * dx[1] * dx[2]/6;
		M[YYZZ] += m4 * dx[1] * dx[1] * dx[2] * dx[2]/4;
		M[YZZZ] += m4 * dx[1] * dx[2] * dx[2] * dx[2]/6;
		M[ZZZZ] += m4 * dx[2] * dx[2] * dx[2] * dx[2]/24;

#endif


	}

}


void m2m(double dx, double dy, double dz, double M[], double tM[]) 
{
	double toM[NMULTI];

	toM[0] = M[0];

	toM[X] = M[X]+M[0]*dx;
	toM[Y] = M[Y]+M[0]*dy;
	toM[Z] = M[Z]+M[0]*dz;

#ifdef QUADRUPOLE

	toM[XX] = M[XX] + M[X]*dx + M[0]*dx*dx/2;
	toM[XY] = M[XY] + M[X]*dy + M[Y]*dx + M[0]*dx*dy;
	toM[XZ] = M[XZ] + M[X]*dz + M[Z]*dx + M[0]*dx*dz;
	toM[YY] = M[YY] + M[Y]*dy + M[0]*dy*dy/2;
	toM[YZ] = M[YZ] + M[Y]*dz + M[Z]*dy + M[0]*dy*dz; 
	toM[ZZ] = M[ZZ] + M[Z]*dz + M[0]*dz*dz/2;

#endif

#ifdef  OCTUPOLE

	toM[XXX] = M[XXX] + M[XX]*dx + M[X]*dx*dx/2 + M[0]*dx*dx*dx/6;
	toM[YYY] = M[YYY] + M[YY]*dy + M[Y]*dy*dy/2 + M[0]*dy*dy*dy/6;
	toM[ZZZ] = M[ZZZ] + M[ZZ]*dz + M[Z]*dz*dz/2 + M[0]*dz*dz*dz/6;
	toM[XXY] = M[XXY] + M[XX]*dy + M[XY]*dx + M[X]*dx*dy + M[Y]*dx*dx/2 + M[0]*dx*dx*dy/2;
	toM[XXZ] = M[XXZ] + M[XX]*dz + M[XZ]*dx + M[X]*dx*dz + M[Z]*dx*dx/2 + M[0]*dx*dx*dz/2;
	toM[XYY] = M[XYY] + M[YY]*dx + M[XY]*dy + M[Y]*dx*dy + M[X]*dy*dy/2 + M[0]*dx*dy*dy/2;
	toM[XZZ] = M[XZZ] + M[ZZ]*dx + M[XZ]*dz + M[Z]*dx*dz + M[X]*dz*dz/2 + M[0]*dx*dz*dz/2;
	toM[YZZ] = M[YZZ] + M[ZZ]*dy + M[YZ]*dz + M[Z]*dy*dz + M[Y]*dz*dz/2 + M[0]*dy*dz*dz/2;
	toM[YYZ] = M[YYZ] + M[YY]*dz + M[YZ]*dy + M[Y]*dy*dz + M[Z]*dy*dy/2 + M[0]*dy*dy*dz/2;

	toM[XYZ] = M[XYZ] + M[XY]*dz + M[YZ]*dx + M[XZ]*dy + M[X]*dy*dz + M[Y]*dx*dz + M[Z]*dy*dx + M[0]*dx*dy*dz;

#endif 


#ifdef  HEXADECAPOLE
	toM[XXXX] = M[XXXX] + M[XXX]*dx + M[XX]*dx*dx/2 + M[X]*dx*dx*dx/6 + M[0]*dx*dx*dx*dx/24;
	toM[YYYY] = M[YYYY] + M[YYY]*dy + M[YY]*dy*dy/2 + M[Y]*dy*dy*dy/6 + M[0]*dy*dy*dy*dy/24;
	toM[ZZZZ] = M[ZZZZ] + M[ZZZ]*dz + M[ZZ]*dz*dz/2 + M[Z]*dz*dz*dz/6 + M[0]*dz*dz*dz*dz/24;

	toM[XXXY] = M[XXXY] + M[XXX]*dy + M[XXY]*dx + M[XX]*dx*dy + M[XY]*dx*dx/2 +  M[X]*dx*dx*dy/2 + M[Y]*dx*dx*dx/6 +  M[0]*dx*dx*dx*dy/6 ;
	toM[XXXZ] = M[XXXZ] + M[XXX]*dz + M[XXZ]*dx + M[XX]*dx*dz + M[XZ]*dx*dx/2 +  M[X]*dx*dx*dz/2 + M[Z]*dx*dx*dx/6 +  M[0]*dx*dx*dx*dz/6 ;
	toM[XYYY] = M[YYYX] + M[YYY]*dx + M[YYX]*dy + M[YY]*dy*dx + M[XY]*dy*dy/2 +  M[Y]*dy*dy*dx/2 + M[X]*dy*dy*dy/6 +  M[0]*dy*dy*dy*dx/6 ;
	toM[YYYZ] = M[YYYZ] + M[YYY]*dz + M[YYZ]*dy + M[YY]*dy*dz + M[ZY]*dy*dy/2 +  M[Y]*dy*dy*dz/2 + M[Z]*dy*dy*dy/6 +  M[0]*dy*dy*dy*dz/6 ;
	toM[XZZZ] = M[ZZZX] + M[ZZZ]*dx + M[ZZX]*dz + M[ZZ]*dz*dx + M[XZ]*dz*dz/2 +  M[Z]*dz*dz*dx/2 + M[X]*dz*dz*dz/6 +  M[0]*dz*dz*dz*dx/6 ;
	toM[YZZZ] = M[ZZZY] + M[ZZZ]*dy + M[ZZY]*dz + M[ZZ]*dz*dy + M[YZ]*dz*dz/2 +  M[Z]*dz*dz*dy/2 + M[Y]*dz*dz*dz/6 +  M[0]*dz*dz*dz*dy/6 ;

	toM[XXYY] = M[XXYY] + M[XXY]*dy + M[XYY]*dx + M[XX]*dy*dy/2 + M[YY]*dx*dx/2 + M[XY]*dx*dy + M[X]*dx*dy*dy/2 + M[Y]*dx*dx*dy/2 + M[0]*dx*dx*dy*dy/4;
	toM[XXZZ] = M[XXZZ] + M[XXZ]*dz + M[XZZ]*dx + M[XX]*dz*dz/2 + M[ZZ]*dx*dx/2 + M[XZ]*dx*dz + M[X]*dx*dz*dz/2 + M[Z]*dx*dx*dz/2 + M[0]*dx*dx*dz*dz/4;
	toM[YYZZ] = M[YYZZ] + M[YYZ]*dz + M[YZZ]*dy + M[YY]*dz*dz/2 + M[ZZ]*dy*dy/2 + M[YZ]*dy*dz + M[Y]*dy*dz*dz/2 + M[Z]*dy*dy*dz/2 + M[0]*dy*dy*dz*dz/4;

	toM[XXYZ] = M[XXYZ] + M[XXY]*dz + M[XXZ]*dy + M[XYZ]*dx + M[XX]*dy*dz + M[XY]*dx*dz + M[XZ]*dx*dy + M[YZ]*dx*dx/2 + M[X]*dx*dy*dz + M[Y]*dx*dx*dz/2 + M[Z]*dx*dx*dy/2 + M[0]*dx*dx*dy*dz/2;
	toM[YYXZ] = M[YYXZ] + M[YYX]*dz + M[YYZ]*dx + M[XYZ]*dy + M[YY]*dx*dz + M[XY]*dy*dz + M[YZ]*dx*dy + M[XZ]*dy*dy/2 + M[Y]*dx*dy*dz + M[X]*dy*dy*dz/2 + M[Z]*dy*dy*dx/2 + M[0]*dy*dy*dx*dz/2;
	toM[ZZXY] = M[ZZXY] + M[ZZX]*dy + M[ZZY]*dx + M[XYZ]*dz + M[ZZ]*dx*dy + M[XZ]*dz*dy + M[YZ]*dx*dz + M[XY]*dz*dz/2 + M[Z]*dx*dz*dy + M[X]*dz*dz*dy/2 + M[Y]*dz*dz*dx/2 + M[0]*dz*dz*dx*dy/2;
#endif


	int n;
	for (n=0; n<NMULTI; n++)
		tM[n] += toM[n];

}




void walk_m2m(int iNode)
{
	int n, idx;
	double *Mptr;
	double dx, dy, dz;

	for (n=0; n<NSON; n++)
	{
		idx = btree[iNode].son[n];
		if (idx < first_node) {
			/// qwang ????
			dx =  btree[iNode].center[0] - leaf[idx].center[0];
			dy =  btree[iNode].center[1] - leaf[idx].center[1];
			dz =  btree[iNode].center[2] - leaf[idx].center[2];


			m2m(dx, dy, dz, leaf[idx].M, btree[iNode].M);
		}
		else {
			walk_m2m(idx);

			dx =  btree[iNode].center[0] - btree[idx].center[0];
			dy =  btree[iNode].center[1] - btree[idx].center[1];
			dz =  btree[iNode].center[2] - btree[idx].center[2];

			m2m(dx, dy, dz, btree[idx].M, btree[iNode].M);
		}
	}

}


void l2p( int ileaf )
{
	double dx, dy, dz;
	int n;
	double *Fn = leaf[ileaf].L;
	double pot, pot0,pot1,pot2, pot3;
	double Fx, Fy, Fz;
	for (n=leaf[ileaf].ipart; n<leaf[ileaf].ipart+leaf[ileaf].npart; n++) 
	{

		dx = part[n].pos[0] - leaf[ileaf].center[0];
		dy = part[n].pos[1] - leaf[ileaf].center[1];
		dz = part[n].pos[2] - leaf[ileaf].center[2];


		pot = Fn[0];


		pot += (Fn[X]*dx+Fn[Y]*dy+Fn[Z]*dz);

		Fx = (Fn[X]);
		Fy = (Fn[Y]);
		Fz = (Fn[Z]);


#ifdef QUADRUPOLE
		pot += 0.5*(Fn[XX]*dx*dx + 2*Fn[XY]*dx*dy+ 2*Fn[XZ]*dx*dz+Fn[YY]*dy*dy + 2*Fn[YZ]*dy*dz+Fn[ZZ]*dz*dz);

		Fx += (Fn[XX]*dx+Fn[XY]*dy+Fn[XZ]*dz);
		Fy += (Fn[XY]*dx+Fn[YY]*dy+Fn[YZ]*dz);
		Fz += (Fn[XZ]*dx+Fn[YZ]*dy+Fn[ZZ]*dz);
#endif 

#ifdef  OCTUPOLE
		pot += (Fn[XXX]*dx*dx*dx + 3*Fn[XXY]*dx*dx*dy + 3*Fn[XXZ]*dx*dx*dz + 6*Fn[XYZ]*dx*dy*dz + 3*Fn[XYY]*dx*dy*dy + 3*Fn[XZZ]*dx*dz*dz + Fn[YYY]*dy*dy*dy + 3*Fn[YYZ]*dy*dy*dz + 3*Fn[YZZ]*dy*dz*dz + Fn[ZZZ]*dz*dz*dz )/6.0;

		Fx +=Fn[XXX]*dx*dx/2 + Fn[XXY]*dx*dy + Fn[XXZ]*dx*dz + Fn[XYZ]*dy*dz + Fn[XYY]*dy*dy/2 + Fn[XZZ]*dz*dz/2;
		Fy +=Fn[YYY]*dy*dy/2 + Fn[YYZ]*dy*dz + Fn[XYY]*dx*dy + Fn[XYZ]*dx*dz + Fn[YZZ]*dz*dz/2 + Fn[XXY]*dx*dx/2;
		Fz +=Fn[ZZZ]*dz*dz/2 + Fn[XZZ]*dx*dz + Fn[YZZ]*dy*dz + Fn[XYZ]*dx*dy + Fn[XXZ]*dx*dx/2 + Fn[YYZ]*dy*dy/2;
#endif

#ifdef  HEXADECAPOLE
		pot += ( Fn[XXXX]*dx*dx*dx*dx + Fn[YYYY]*dy*dy*dy*dy + Fn[ZZZZ]*dz*dz*dz*dz + 4*Fn[XXXY]*dx*dx*dx*dy + 4*Fn[XXXZ]*dx*dx*dx*dz + 4*Fn[ZZZY]*dz*dz*dz*dy + 4*Fn[ZZZX]*dx*dz*dz*dz + 4*Fn[YYYX]*dy*dy*dy*dx + 4*Fn[YYYZ]*dy*dy*dy*dz  + 6*Fn[XXYY]*dx*dx*dy*dy + 6*Fn[XXZZ]*dx*dx*dz*dz + 6*Fn[YYZZ]*dy*dy*dz*dz + 12*Fn[XXYZ]*dx*dx*dy*dz + 12*Fn[XYYZ]*dx*dy*dy*dz + 12*Fn[XYZZ]*dz*dx*dy*dz )/24.0;

		Fx += Fn[XXXX]*dx*dx*dx/6 + Fn[XXXY]*dx*dx*dy/2 + Fn[XXXZ]*dx*dx*dz/2 + Fn[XXYZ]*dx*dy*dz + Fn[XXYY]*dx*dy*dy/2 + Fn[XXZZ]*dx*dz*dz/2 + Fn[XYYY]*dy*dy*dy/6 + Fn[XYYZ]*dy*dy*dz/2 + Fn[XYZZ]*dy*dz*dz/2 + Fn[XZZZ]*dz*dz*dz/6;
		Fy += Fn[YYYY]*dy*dy*dy/6 + Fn[YYYZ]*dy*dy*dz/2 + Fn[YYYX]*dy*dy*dx/2 + Fn[YYZX]*dy*dz*dx + Fn[YYZZ]*dy*dz*dz/2 + Fn[YYXX]*dy*dx*dx/2 + Fn[YZZZ]*dz*dz*dz/6 + Fn[YZZX]*dz*dz*dx/2 + Fn[YZXX]*dz*dx*dx/2 + Fn[YXXX]*dx*dx*dx/6;
		Fz += Fn[ZZZZ]*dz*dz*dz/6 + Fn[ZZZX]*dz*dz*dx/2 + Fn[ZZZY]*dz*dz*dy/2 + Fn[ZZXY]*dz*dx*dy + Fn[ZZXX]*dz*dx*dx/2 + Fn[ZZYY]*dz*dy*dy/2 + Fn[ZXXX]*dx*dx*dx/6 + Fn[ZXXY]*dx*dx*dy/2 + Fn[ZXYY]*dx*dy*dy/2 + Fn[ZYYY]*dy*dy*dy/6;
#endif

		part[n].acc[0] += Fx;
		part[n].acc[1] += Fy;
		part[n].acc[2] += Fz;
		// part[n].acc[3] -= pot;
	}
}



void m2l(double dx, double dy, double dz, double M[], double toL[]) 
{
	double dr, r2;
	double ir, ir2, ir3, ir4, ir5, ir6, ir7, ir8,ir9;
	double rs2, irs2, irs3, irs, irs5, irs7;
	double Fn[NMULTI];
	double Dn[NMULTI];

	walk_m2l_count++;
	//	printf(" m2l");

	r2 = (dx*dx + dy*dy + dz*dz);
	dr = sqrt(r2);	

	ir  = 1.0/dr;
	ir2 = ir*ir;
	ir3 = ir*ir*ir;
	ir4 = ir3 * ir;
	ir5 = ir*ir*ir3;
	ir6 = ir2*ir4;
	ir7 = ir5*ir2;
	ir8 = ir7 * ir;
	ir9 = ir7*ir2; 

	irs = 1.0/splitRadius;
	rs2 = splitRadius*splitRadius;
	irs2 = irs *irs;
	irs3 = irs2*irs;
	irs5 = irs3*irs2;
	irs7 = irs5*irs2;

	double fac[5];

	fac[0] = ir;
	fac[1] =-ir3;
	fac[2] = 3.0*ir5;
	fac[3] =-3.0*5.0*ir7;
	fac[4] = 3.0*5.0*7.0*ir9;

#ifdef LONGSHORT

	double drs = 0.5*dr/splitRadius;
	double coeff = 1.0/sqrt(M_PI);
	double facExp = exp(-drs*drs) * coeff;
	double facErc = erfc(drs);

	fac[0] =  ir * facErc ; 
	fac[1] = -ir3*(facErc + dr * facExp * irs );
	fac[2] =  3.0*ir5*facErc + (3.0*irs*ir4  + 0.5*ir2*irs3)*facExp;
	fac[3] = -3.0*5.0*ir7*facErc-(15.0*ir6*irs + 2.5*ir4*irs3 +0.25*ir2*irs5)*facExp;
	fac[4] =  3.0*5.0*7.0*ir9*facErc + (3.0*5.0*7.0*ir8*irs + 7.5*ir6*irs3 + 10*ir6*irs3 +1.25*ir4*irs5 + 0.5*ir4* irs5 +0.125*ir2*irs7 )*facExp;

#endif

	Dn[0]  = fac[0];
	Dn[X]  = fac[1] * dx;
	Dn[Y]  = fac[1] * dy;
	Dn[Z]  = fac[1] * dz;

	Fn[0]  = M[0]*Dn[0] + M[X]*Dn[X] + M[Y]*Dn[Y] + M[Z]*Dn[Z] ;
	Fn[X]  = M[0]*Dn[X]  ;
	Fn[Y]  = M[0]*Dn[Y]  ;
	Fn[Z]  = M[0]*Dn[Z]  ;

#ifdef QUADRUPOLE
	Dn[XX] =  fac[2]*dx*dx + fac[1];
	Dn[XY] =  fac[2]*dx*dy;
	Dn[XZ] =  fac[2]*dx*dz;
	Dn[YY] =  fac[2]*dy*dy + fac[1];
	Dn[YZ] =  fac[2]*dy*dz;
	Dn[ZZ] =  fac[2]*dz*dz + fac[1];

	Fn[0]   += M[XX]*Dn[XX] +  M[XY] *Dn[XY] + M[XZ]*Dn[XZ] + M[YY]*Dn[YY] + M[YZ]*Dn[YZ] + M[ZZ]*Dn[ZZ];

	Fn[X]   += M[X]*Dn[XX] + M[Y]*Dn[XY] + M[Z]*Dn[XZ] ;
	Fn[Y]   += M[X]*Dn[YX] + M[Y]*Dn[YY] + M[Z]*Dn[YZ] ;
	Fn[Z]   += M[X]*Dn[ZX] + M[Y]*Dn[ZY] + M[Z]*Dn[ZZ] ;

	Fn[XX]  = M[0]*Dn[XX];
	Fn[XY]  = M[0]*Dn[XY];
	Fn[XZ]  = M[0]*Dn[XZ];
	Fn[YY]  = M[0]*Dn[YY];
	Fn[YZ]  = M[0]*Dn[YZ];
	Fn[ZZ]  = M[0]*Dn[ZZ];

#endif

#ifdef  OCTUPOLE

	Dn[XXX]= fac[3]*dx*dx*dx + fac[2]*dx*3;
	Dn[XYY]= fac[3]*dx*dy*dy + fac[2]*dx;
	Dn[XZZ]= fac[3]*dx*dz*dz + fac[2]*dx;

	Dn[XXY]= fac[3]*dx*dx*dy + fac[2]*dy;
	Dn[YYY]= fac[3]*dy*dy*dy + fac[2]*dy*3;
	Dn[YZZ]= fac[3]*dy*dz*dz + fac[2]*dy;

	Dn[XXZ]= fac[3]*dx*dx*dz + fac[2]*dz;
	Dn[YYZ]= fac[3]*dy*dy*dz + fac[2]*dz;
	Dn[ZZZ]= fac[3]*dz*dz*dz + fac[2]*dz*3;

	Dn[XYZ]= fac[3]*dx*dy*dz;

	Fn[0]   += M[XXX]*Dn[XXX] + M[XXY]*Dn[XXY] + M[XXZ]*Dn[XXZ] + M[XYY]*Dn[XYY] + M[XYZ]*Dn[XYZ] + M[YYY]*Dn[YYY] + M[YYZ]*Dn[YYZ] + M[YZZ]*Dn[YZZ] + M[XZZ]*Dn[XZZ] + M[ZZZ]*Dn[ZZZ];

	Fn[X]   += M[XX]*Dn[XXX] + M[XY]*Dn[XXY] + M[XZ]*Dn[XXZ] + M[YY]*Dn[XYY] + M[YZ]*Dn[XYZ] + M[ZZ]*Dn[XZZ] ;
	Fn[Y]   += M[YY]*Dn[YYY] + M[YZ]*Dn[YYZ] + M[YX]*Dn[YYX] + M[ZZ]*Dn[YZZ] + M[ZX]*Dn[YZX] + M[XX]*Dn[YXX] ;
	Fn[Z]   += M[ZZ]*Dn[ZZZ] + M[ZX]*Dn[ZZX] + M[ZY]*Dn[ZZY] + M[XX]*Dn[ZXX] + M[XY]*Dn[ZXY] + M[YY]*Dn[ZYY] ;

	Fn[XX]  += M[X]*Dn[XXX] + M[Y]*Dn[XXY] + M[Z]*Dn[XXZ];
	Fn[XY]  += M[X]*Dn[XYX] + M[Y]*Dn[XYY] + M[Z]*Dn[XYZ];
	Fn[XZ]  += M[X]*Dn[XZX] + M[Y]*Dn[XZY] + M[Z]*Dn[XZZ];
	Fn[YY]  += M[X]*Dn[YYX] + M[Y]*Dn[YYY] + M[Z]*Dn[YYZ];
	Fn[YZ]  += M[X]*Dn[YZX] + M[Y]*Dn[YZY] + M[Z]*Dn[YZZ];
	Fn[ZZ]  += M[X]*Dn[ZZX] + M[Y]*Dn[ZZY] + M[Z]*Dn[ZZZ];

	Fn[XXX] = M[0]*Dn[XXX];
	Fn[XXY] = M[0]*Dn[XXY];
	Fn[XXZ] = M[0]*Dn[XXZ];
	Fn[XYY] = M[0]*Dn[XYY];
	Fn[XYZ] = M[0]*Dn[XYZ];
	Fn[XZZ] = M[0]*Dn[XZZ];
	Fn[YYY] = M[0]*Dn[YYY];
	Fn[YYZ] = M[0]*Dn[YYZ];
	Fn[YZZ] = M[0]*Dn[YZZ];
	Fn[ZZZ] = M[0]*Dn[ZZZ];

#endif



	int n;
	for (n=0; n<NMULTI; n++) {
		toL[n] += Fn[n];
	}


}


void l2l(double dx, double dy, double dz, double L[], double toL[]) 
{
	toL[0] += L[0]+L[X]*dx+L[Y]*dy+L[Z]*dz;

	toL[X] += L[X];
	toL[Y] += L[Y];
	toL[Z] += L[Z];

#ifdef QUADRUPOLE


	toL[0] += L[XX]*dx*dx/2+L[XY]*dx*dy+L[XZ]*dx*dz+L[YZ]*dy*dz+L[YY]*dy*dy/2+L[ZZ]*dz*dz/2;

	toL[X] += L[XX]*dx+L[XY]*dy+L[XZ]*dz;
	toL[Y] += L[YX]*dx+L[YY]*dy+L[YZ]*dz;
	toL[Z] += L[ZX]*dx+L[ZY]*dy+L[ZZ]*dz;

	toL[XX] += L[XX] ;
	toL[XY] += L[XY] ;
	toL[XZ] += L[XZ] ;
	toL[YY] += L[YY] ;
	toL[YZ] += L[YZ] ;
	toL[ZZ] += L[ZZ] ;
#endif

#ifdef  OCTUPOLE

	toL[0] += L[XXX]*dx*dx*dx/6+L[XXY]*dx*dx*dy/2+L[XXZ]*dx*dx*dz/2+L[XYY]*dx*dy*dy/2+L[XYZ]*dx*dy*dz+L[XZZ]*dx*dz*dz/2+L[YYY]*dy*dy*dy/6+L[YZZ]*dy*dz*dz/2+L[YYZ]*dy*dy*dz/2+L[ZZZ]*dz*dz*dz/6;

	toL[X] += L[XXX]*dx*dx*0.5+L[XXY]*dx*dy+L[XXZ]*dx*dz+L[XYY]*dy*dy*0.5+L[XYZ]*dy*dz+L[XZZ]*dz*dz*0.5;
	toL[Y] += L[YXX]*dx*dx*0.5+L[YXY]*dx*dy+L[YXZ]*dx*dz+L[YYY]*dy*dy*0.5+L[YYZ]*dy*dz+L[YZZ]*dz*dz*0.5;
	toL[Z] += L[ZXX]*dx*dx*0.5+L[ZXY]*dx*dy+L[ZXZ]*dx*dz+L[ZYY]*dy*dy*0.5+L[ZYZ]*dy*dz+L[ZZZ]*dz*dz*0.5;

	toL[XX] += L[XXX]*dx + L[XXY]*dy + L[XXZ]*dz;
	toL[XY] += L[XYX]*dx + L[XYY]*dy + L[XYZ]*dz;
	toL[XZ] += L[XZX]*dx + L[XZY]*dy + L[XZZ]*dz;
	toL[YY] += L[YYX]*dx + L[YYY]*dy + L[YYZ]*dz;
	toL[YZ] += L[YZX]*dx + L[YZY]*dy + L[YZZ]*dz;
	toL[ZZ] += L[ZZX]*dx + L[ZZY]*dy + L[ZZZ]*dz;

	toL[XXX] += L[XXX];
	toL[XZZ] += L[XZZ];
	toL[XYY] += L[XYY];
	toL[XXY] += L[XXY];
	toL[XXZ] += L[XXZ];
	toL[YYY] += L[YYY];
	toL[XYZ] += L[XYZ];
	toL[YYZ] += L[YYZ];
	toL[YZZ] += L[YZZ];
	toL[ZZZ] += L[ZZZ];
#endif

#ifdef  HEXADECAPOLE
	toL[XXXX] += L[XXXX];
	toL[YYYY] += L[YYYY];
	toL[ZZZZ] += L[ZZZZ];

	toL[XXXY] += L[XXXY];
	toL[XXXZ] += L[XXXZ];
	toL[XYYY] += L[XYYY];
	toL[YYYZ] += L[YYYZ];
	toL[XZZZ] += L[XZZZ];
	toL[YZZZ] += L[YZZZ];

	toL[XXYY] += L[XXYY];
	toL[XXZZ] += L[XXZZ];
	toL[YYZZ] += L[YYZZ];

	toL[XXYZ] += L[XXYZ];
	toL[YYXZ] += L[YYXZ];
	toL[ZZXY] += L[ZZXY];

	toL[XXX] += L[XXXX]*dx + L[XXXY]*dy + L[XXXZ]*dz;
	toL[YYY] += L[YYYX]*dx + L[YYYY]*dy + L[YYYZ]*dz;
	toL[ZZZ] += L[ZZZX]*dx + L[ZZZY]*dy + L[ZZZZ]*dz;
	toL[XXY] += L[XXYX]*dx + L[XXYY]*dy + L[XXYZ]*dz;
	toL[XXZ] += L[XXZX]*dx + L[XXZY]*dy + L[XXZZ]*dz;
	toL[XYY] += L[XYYX]*dx + L[XYYY]*dy + L[XYYZ]*dz;
	toL[YYZ] += L[YYZX]*dx + L[YYZY]*dy + L[YYZZ]*dz;
	toL[XZZ] += L[XZZX]*dx + L[XZZY]*dy + L[XZZZ]*dz;
	toL[YZZ] += L[YZZX]*dx + L[YZZY]*dy + L[YZZZ]*dz;
	toL[XYZ] += L[XYZX]*dx + L[XYZY]*dy + L[XYZZ]*dz;

	toL[XX] += L[XXXX]*dx*dx/2 + L[XXYY]*dy*dy/2 + L[XXZZ]*dz*dz/2 + L[XXXY]*dx*dy + L[XXXZ]*dx*dz + L[XXYZ]*dy*dz ;
	toL[YY] += L[YYXX]*dx*dx/2 + L[YYYY]*dy*dy/2 + L[YYZZ]*dz*dz/2 + L[YYYX]*dy*dx + L[YYZX]*dz*dx + L[YYYZ]*dy*dz ;
	toL[ZZ] += L[ZZXX]*dx*dx/2 + L[ZZYY]*dy*dy/2 + L[ZZZZ]*dz*dz/2 + L[ZZXY]*dx*dy + L[ZZZX]*dz*dx + L[ZZZY]*dz*dy ;
	toL[XY] += L[XYXX]*dx*dx/2 + L[XYYY]*dy*dy/2 + L[XYZZ]*dz*dz/2 + L[XYXY]*dx*dy + L[XYXZ]*dx*dz + L[XYYZ]*dy*dz ;
	toL[XZ] += L[XZXX]*dx*dx/2 + L[XZYY]*dy*dy/2 + L[XZZZ]*dz*dz/2 + L[XZXY]*dx*dy + L[XZXZ]*dx*dz + L[XZYZ]*dy*dz ;
	toL[YZ] += L[YZXX]*dx*dx/2 + L[YZYY]*dy*dy/2 + L[YZZZ]*dz*dz/2 + L[YZYX]*dy*dx + L[YZZX]*dz*dx + L[YZYZ]*dy*dz  ;

	toL[X]  += L[XXXX]*dx*dx*dx/6 + L[XXXY]*dx*dx*dy/2 + L[XXXZ]*dx*dx*dz/2 + L[XXYY]*dx*dy*dy/2 + L[XXYZ]*dx*dy*dz + L[XXZZ]*dx*dz*dz/2 + L[XYYY]*dy*dy*dy/6 + L[XYYZ]*dy*dy*dz/2 + L[XYZZ]*dy*dz*dz/2 + L[XZZZ]*dz*dz*dz/6;
	toL[Y]  += L[YXXX]*dx*dx*dx/6 + L[YXXY]*dx*dx*dy/2 + L[YXXZ]*dx*dx*dz/2 + L[YXYY]*dx*dy*dy/2 + L[YXYZ]*dx*dy*dz + L[YXZZ]*dx*dz*dz/2 + L[YYYY]*dy*dy*dy/6 + L[YYYZ]*dy*dy*dz/2 + L[YYZZ]*dy*dz*dz/2 + L[YZZZ]*dz*dz*dz/6;
	toL[Z]  += L[ZXXX]*dx*dx*dx/6 + L[ZXXY]*dx*dx*dy/2 + L[ZXXZ]*dx*dx*dz/2 + L[ZXYY]*dx*dy*dy/2 + L[ZXYZ]*dx*dy*dz + L[ZXZZ]*dx*dz*dz/2 + L[ZYYY]*dy*dy*dy/6 + L[ZYYZ]*dy*dy*dz/2 + L[ZYZZ]*dy*dz*dz/2 + L[ZZZZ]*dz*dz*dz/6;


	toL[0] += L[XXXX]*dx*dx*dx*dx/24 + L[XXXY]*dx*dx*dx*dy/6 + L[XXYY]*dx*dx*dy*dy/4 + L[XXYZ]*dx*dx*dy*dz/2 + L[XXZZ]*dx*dx*dz*dz/4 + L[XYYY]*dx*dy*dy*dy/6 + L[XYYZ]*dx*dy*dy*dz/2 + L[XYZZ]*dx*dy*dz*dz/2 + L[XZZZ]*dx*dz*dz*dz/6 + L[YYYY]*dy*dy*dy*dy/24 + L[YYYZ]*dy*dy*dy*dz/6 + L[YYZZ]*dy*dy*dz*dz/4 + L[YZZZ]*dy*dz*dz*dz/6 + L[ZXXX]*dz*dx*dx*dx/6 + L[ZZZZ]*dz*dz*dz*dz/24;

#endif

}



void walk_l2l(int inode)
{
	if ( inode < first_node ) {
		return;
	}

	int n, idx;
	double dx, dy ,dz;
	for (n=0; n<NSON; n++) {
		idx = btree[inode].son[n];
		if ( idx  >= first_node ) {
			dx = btree[idx].center[0] - btree[inode].center[0];
			dy = btree[idx].center[1] - btree[inode].center[1];
			dz = btree[idx].center[2] - btree[inode].center[2];

			l2l(dx, dy, dz, btree[inode].L, btree[idx].L) ;
			walk_l2l( idx );
		}
		else if (idx >= first_leaf) {

			dx = leaf[idx].center[0] - btree[inode].center[0];
			dy = leaf[idx].center[1] - btree[inode].center[1];
			dz = leaf[idx].center[2] - btree[inode].center[2];

			l2l(dx, dy, dz, btree[inode].L, leaf[idx].L) ;
		}
		else
			return;
	}

}



