#ifndef OPERATOR_H
#define OPERATOR_H

#include "kernels.h"

void p2m(int iPart, int npart, double center[3], double M[]);
void m2m(double dx, double dy, double dz, double M[], double tM[]);
void m2l(double dx, double dy, double dz, double M[], double toL[]);
void l2l(double dx, double dy, double dz, double L[], double toL[]);
void l2p( int ileaf );


void l2p( int ileaf );

void walk_m2m(int iNode);
void walk_m2l(int im, int jm);
void walk_m2l_task(int im, int jm);
void walk_m2l_tasknum(int im, int jm);
void walk_l2l(int inode);

void driver_fmm(int NPART, double BOXSIZE, int max_leaf, double open, int norder);


#define X 	1
#define Y 	2
#define Z 	3

#define XX	4
#define XY	5
#define YX	5
#define XZ	6
#define ZX	6
#define YY	7
#define YZ	8
#define ZY	8
#define ZZ	9


#define XXX	10
#define XYX	11
#define YXX	11
#define XZX	12
#define ZXX	12
#define YYX	13
#define YZX	14
#define ZYX	14
#define ZZX	15

#define XXY	11
#define XYY	13
#define YXY	13
#define XZY	14
#define ZXY	14
#define YYY	16
#define YZY	17
#define ZYY	17
#define ZZY	18

#define XXZ	12
#define XYZ	14
#define YXZ	14
#define XZZ	15
#define ZXZ	15
#define YYZ	17
#define YZZ	18
#define ZYZ	18
#define ZZZ	19


#define XXXX	20

#define XXYX	21
#define YXXX	21
#define XYXX	21
#define XXXY	21

#define XXZX	22
#define ZXXX	22
#define XZXX	22
#define XXXZ	22

#define XYYX	23
#define YXYX	23
#define YYXX	23
#define XXYY	23
#define YXXY	23
#define XYXY	23

#define XYZX	24
#define XZYX	24
#define YXZX	24
#define YZXX	24
#define ZXYX	24
#define ZYXX    24
#define XXZY	24
#define ZXXY	24
#define XZXY	24
#define XXYZ	24
#define YXXZ	24
#define XYXZ	24

#define XZZX	25
#define ZZXX	25
#define ZXZX	25
#define XXZZ	25
#define ZXXZ	25
#define XZXZ	25

#define YYYX	26
#define XYYY	26
#define YXYY	26
#define YYXY	26

#define YYZX	27	
#define YZYX	27	
#define ZYYX	27
#define XYZY	27
#define XZYY	27
#define YXZY	27
#define YZXY	27
#define ZXYY	27
#define ZYXY    27
#define XYYZ	27
#define YXYZ	27
#define YYXZ	27

#define YZZX	28
#define ZZYX	28
#define ZYZX	28
#define XZZY	28
#define ZZXY	28
#define ZXZY	28
#define XYZZ	28
#define XZYZ	28
#define YXZZ	28
#define YZXZ	28
#define ZXYZ	28
#define ZYXZ    28

#define ZZZX    29
#define XZZZ	29
#define ZZXZ	29
#define ZXZZ	29

#define YYYY	30

#define YYZY	31	
#define YZYY	31	
#define ZYYY	31
#define YYYZ	31

#define YYZZ	32
#define YZYZ	32	
#define ZYYZ	32
#define YZZY	32
#define ZZYY	32
#define ZYZY	32

#define ZZZY    33
#define YZZZ	33
#define ZZYZ	33
#define ZYZZ	33

#define ZZZZ    34

#endif

