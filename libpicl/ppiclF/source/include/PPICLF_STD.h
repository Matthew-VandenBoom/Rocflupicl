#include "PPICLF_USER.h"

#ifndef PPICLF_LPART
#define PPICLF_LPART 100000
#endif

#ifndef NDIM
#define NDIM 3
#endif

#ifdef PPICLF_LRP
! #define PPICLF_LRP 1
#undef PPICLF_LRP
#endif

#ifdef PPICLF_LRP2
! #define PPICLF_LRP2 0
#undef PPICLF_LRP2
#endif

#ifndef PPICLF_VU
#define PPICLF_VU 0
! #undef PPICLF_VU
#endif

#ifdef PPICLF_LRP3
! #define PPICLF_LRP3 0 
#undef PPICLF_LRP3
#endif

#ifdef PPICLF_LRP4
! #define PPICLF_LRP4 0
#undef PPICLF_LRP4
#endif

#ifdef PPICLF_LRP5
! #define PPICLF_LRP5 0
#undef PPICLF_LRP5
#endif

#ifdef PPICLF_LIP
! #define PPICLF_LIP 11
#undef PPICLF_LIP
#endif

#ifndef PPICLF_LPART_GP
#define PPICLF_LPART_GP 8*PPICLF_LPART
!#undef PPICLF_LPART_GP
#endif

#ifndef PPICLF_INTERP
#define PPICLF_INTERP 1
! #undef PPICLF_INTERP
#endif

#ifndef PPICLF_LRP_INT
#undef PPICLF_INTERP
#define PPICLF_INTERP 0
#define PPICLF_LRP_INT 1
#endif

#ifndef PPICLF_PROJECT
#define PPICLF_PROJECT 1
! #undef PPICLF_PROJECT
#endif

#ifndef PPICLF_LRP_PRO
#undef PPICLF_PROJECT
#define PPICLF_PROJECT 0
#define PPICLF_LRP_PRO 1
#endif

#ifdef PPICLF_LRP_GP
! #define PPICLF_LRP_GP PPICLF_LRS+PPICLF_LRP+PPICLF_LRP_PRO
#undef PPICLF_LRP_GP
#endif

#ifdef PPICLF_LIP_GP
! #define PPICLF_LIP_GP 8
#undef PPICLF_LIP_GP
#endif

#ifdef PPICLF_LEX
! #define PPICLF_LEX 1
#undef PPICLF_LEX
#endif

#ifdef PPICLF_LEY
! #define PPICLF_LEY 1
#undef PPICLF_LEY
#endif

#ifdef PPICLF_LEZ
! #define PPICLF_LEZ 1
#undef PPICLF_LEZ
#endif

! Maximum number of overlap mesh elements per rank
#ifndef PPICLF_LEE
#define PPICLF_LEE 1
#endif

#ifdef PPICLF_BMAX
! #define PPICLF_BMAX 1
#undef PPICLF_BMAX
#endif

#ifdef PPICLF_BX1
! #define PPICLF_BX1 1
#undef PPICLF_BX1
#endif

#ifdef PPICLF_BY1
! #define PPICLF_BY1 1
#undef PPICLF_BY1
#endif

#ifdef PPICLF_BZ1
! #define PPICLF_BZ1 1
#undef PPICLF_BZ1
#endif

#ifndef PPICLF_LRMAX
#define PPICLF_LRMAX 6
!#undef PPICLF_LRMAX
#endif

! The maximum number of triangular patch boundaries
#ifndef PPICLF_LWALL
#define PPICLF_LWALL 20
!#undef PPICLF_LWALL
#endif


#ifdef __GFORTRAN__
#define PASTE(a) a
#define CAT(a,b) PASTE(a)b
#else
#define PASTE(a) a ## b
#define CAT(a,b) PASTE(a,b)
#endif