//#include <ia32intrin.h>
//#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
//#include <smmintrin.h> //SSE4.1
#include <omp.h>

//#include <nmmintrin.h> //SSE4.2
//****************************************************************************************

#ifndef MAX
	#define MAX( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef MIN
	#define MIN( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#define UP_L8(r1,r2,r3,r4,r5,r6,r7,r8) r8=r7;r7=r6;r6=r5;r5=r4;r4=r3;r3=r2;r2=r1;

//#define COUNT icount+=4.0;
#define COUNT 

#define CONV_L1_INIT(il,ir,ishift,f1,ro1) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	ro1=_mm_mul_ps(x,f1);COUNT;

#define CONV_L1_NEXT(il,ir,ishift,f1,ri1,ro1) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	y=_mm_mul_ps(x,f1);COUNT;ro1=_mm_add_ps(ri1,y);COUNT;

#define CONV_L2_INIT(il,ir,ishift,f1,f2,ro1,ro2) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	ro1=_mm_mul_ps(x,f1);COUNT; \
	ro2=_mm_mul_ps(x,f2);COUNT;

#define CONV_L2_NEXT(il,ir,ishift,f1,f2,ri1,ri2,ro1,ro2) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	y=_mm_mul_ps(x,f1);COUNT;ro1=_mm_add_ps(ri1,y);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT;

#define CONV_L3_INIT(il,ir,ishift,f1,f2,f3,ro1,ro2,ro3) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	ro1=_mm_mul_ps(x,f1);COUNT; \
	ro2=_mm_mul_ps(x,f2);COUNT; \
	ro3=_mm_mul_ps(x,f3);COUNT;

#define CONV_L3_NEXT(il,ir,ishift,f1,f2,f3,ri1,ri2,ri3,ro1,ro2,ro3) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	y=_mm_mul_ps(x,f1);COUNT;ro1=_mm_add_ps(ri1,y);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT;

#define CONV_L4_INIT(il,ir,ishift,f1,f2,f3,f4,ri2,ri3,ri4,ro1,ro2,ro3,ro4) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	ro1=_mm_mul_ps(x,f1);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT; \
	y=_mm_mul_ps(x,f4);COUNT;ro4=_mm_add_ps(ri4,y);COUNT;

#define CONV_L4_NEXT(il,ir,ishift,f1,f2,f3,f4,ri1,ri2,ri3,ri4,ro1,ro2,ro3,ro4) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	y=_mm_mul_ps(x,f1);COUNT;ro1=_mm_add_ps(ri1,y);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT; \
	y=_mm_mul_ps(x,f4);COUNT;ro4=_mm_add_ps(ri4,y);COUNT;

#define CONV_L5_INIT(il,ir,ishift,f1,f2,f3,f4,f5,ri2,ri3,ri4,ri5,ro1,ro2,ro3,ro4,ro5) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	ro1=_mm_mul_ps(x,f1);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT; \
	y=_mm_mul_ps(x,f4);COUNT;ro4=_mm_add_ps(ri4,y);COUNT; \
	y=_mm_mul_ps(x,f5);COUNT;ro5=_mm_add_ps(ri5,y);COUNT;

#define CONV_L5_NEXT(il,ir,ishift,f1,f2,f3,f4,f5,ri1,ri2,ri3,ri4,ri5,ro1,ro2,ro3,ro4,ro5) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	y=_mm_mul_ps(x,f1);COUNT;ro1=_mm_add_ps(ri1,y);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT; \
	y=_mm_mul_ps(x,f4);COUNT;ro4=_mm_add_ps(ri4,y);COUNT; \
	y=_mm_mul_ps(x,f5);COUNT;ro5=_mm_add_ps(ri5,y);COUNT;

#define CONV_L6_INIT(il,ir,ishift,f1,f2,f3,f4,f5,f6,ri2,ri3,ri4,ri5,ri6,ro1,ro2,ro3,ro4,ro5,ro6) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	ro1=_mm_mul_ps(x,f1);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT; \
	y=_mm_mul_ps(x,f4);COUNT;ro4=_mm_add_ps(ri4,y);COUNT; \
	y=_mm_mul_ps(x,f5);COUNT;ro5=_mm_add_ps(ri5,y);COUNT; \
	y=_mm_mul_ps(x,f6);COUNT;ro6=_mm_add_ps(ri6,y);COUNT;

#define CONV_L6_NEXT(il,ir,ishift,f1,f2,f3,f4,f5,f6,ri1,ri2,ri3,ri4,ri5,ri6,ro1,ro2,ro3,ro4,ro5,ro6) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	y=_mm_mul_ps(x,f1);COUNT;ro1=_mm_add_ps(ri1,y);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT; \
	y=_mm_mul_ps(x,f4);COUNT;ro4=_mm_add_ps(ri4,y);COUNT; \
	y=_mm_mul_ps(x,f5);COUNT;ro5=_mm_add_ps(ri5,y);COUNT; \
	y=_mm_mul_ps(x,f6);COUNT;ro6=_mm_add_ps(ri6,y);COUNT;

#define CONV_L7_INIT(il,ir,ishift,f1,f2,f3,f4,f5,f6,f7,ri2,ri3,ri4,ri5,ri6,ri7,ro1,ro2,ro3,ro4,ro5,ro6,ro7) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	ro1=_mm_mul_ps(x,f1);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT; \
	y=_mm_mul_ps(x,f4);COUNT;ro4=_mm_add_ps(ri4,y);COUNT; \
	y=_mm_mul_ps(x,f5);COUNT;ro5=_mm_add_ps(ri5,y);COUNT; \
	y=_mm_mul_ps(x,f6);COUNT;ro6=_mm_add_ps(ri6,y);COUNT; \
	y=_mm_mul_ps(x,f7);COUNT;ro7=_mm_add_ps(ri7,y);COUNT;

#define CONV_L7_NEXT(il,ir,ishift,f1,f2,f3,f4,f5,f6,f7,ri1,ri2,ri3,ri4,ri5,ri6,ri7,ro1,ro2,ro3,ro4,ro5,ro6,ro7) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	y=_mm_mul_ps(x,f1);COUNT;ro1=_mm_add_ps(ri1,y);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT; \
	y=_mm_mul_ps(x,f4);COUNT;ro4=_mm_add_ps(ri4,y);COUNT; \
	y=_mm_mul_ps(x,f5);COUNT;ro5=_mm_add_ps(ri5,y);COUNT; \
	y=_mm_mul_ps(x,f6);COUNT;ro6=_mm_add_ps(ri6,y);COUNT; \
	y=_mm_mul_ps(x,f7);COUNT;ro7=_mm_add_ps(ri7,y);COUNT;

#define CONV_L8_INIT(il,ir,ishift,f1,f2,f3,f4,f5,f6,f7,f8,ri2,ri3,ri4,ri5,ri6,ri7,ri8,ro1,ro2,ro3,ro4,ro5,ro6,ro7,ro8) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	ro1=_mm_mul_ps(x,f1);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT; \
	y=_mm_mul_ps(x,f4);COUNT;ro4=_mm_add_ps(ri4,y);COUNT; \
	y=_mm_mul_ps(x,f5);COUNT;ro5=_mm_add_ps(ri5,y);COUNT; \
	y=_mm_mul_ps(x,f6);COUNT;ro6=_mm_add_ps(ri6,y);COUNT; \
	y=_mm_mul_ps(x,f7);COUNT;ro7=_mm_add_ps(ri7,y);COUNT; \
	y=_mm_mul_ps(x,f8);COUNT;ro8=_mm_add_ps(ri8,y);COUNT;

#define CONV_L8_NEXT(il,ir,ishift,f1,f2,f3,f4,f5,f6,f7,f8,ri1,ri2,ri3,ri4,ri5,ri6,ri7,ri8,ro1,ro2,ro3,ro4,ro5,ro6,ro7,ro8) \
	x=_mm_castsi128_ps(_mm_alignr_epi8(ir,il,ishift)); \
	y=_mm_mul_ps(x,f1);COUNT;ro1=_mm_add_ps(ri1,y);COUNT; \
	y=_mm_mul_ps(x,f2);COUNT;ro2=_mm_add_ps(ri2,y);COUNT; \
	y=_mm_mul_ps(x,f3);COUNT;ro3=_mm_add_ps(ri3,y);COUNT; \
	y=_mm_mul_ps(x,f4);COUNT;ro4=_mm_add_ps(ri4,y);COUNT; \
	y=_mm_mul_ps(x,f5);COUNT;ro5=_mm_add_ps(ri5,y);COUNT; \
	y=_mm_mul_ps(x,f6);COUNT;ro6=_mm_add_ps(ri6,y);COUNT; \
	y=_mm_mul_ps(x,f7);COUNT;ro7=_mm_add_ps(ri7,y);COUNT; \
	y=_mm_mul_ps(x,f8);COUNT;ro8=_mm_add_ps(ri8,y);COUNT;

#define CONV_L1(l1,fm03,fm02,fm01,fp00,fp01,fp02,fp03) \
	CONV_L1_INIT(lz,l1, 4,fm03   ,r0) \
	CONV_L1_NEXT(lz,l1, 8,fm02,r0,r0) \
	CONV_L1_NEXT(lz,l1,12,fm01,r0,r0) \
	CONV_L1_NEXT(l1,lz, 0,fp00,r0,r0) \
	CONV_L1_NEXT(l1,lz, 4,fp01,r0,r0) \
	CONV_L1_NEXT(l1,lz, 8,fp02,r0,r0) \
	CONV_L1_NEXT(l1,lz,12,fp03,r0,r0)

#define CONV_L1C(l1,fm03,fm02,fm01,fp00,fp01,fp02,fp03) \
	CONV_L1_NEXT(lz,l1, 4,fm03,r0,r0) \
	CONV_L1_NEXT(lz,l1, 8,fm02,r0,r0) \
	CONV_L1_NEXT(lz,l1,12,fm01,r0,r0) \
	CONV_L1_NEXT(l1,lz, 0,fp00,r0,r0) \
	CONV_L1_NEXT(l1,lz, 4,fp01,r0,r0) \
	CONV_L1_NEXT(l1,lz, 8,fp02,r0,r0) \
	CONV_L1_NEXT(l1,lz,12,fp03,r0,r0)

#define CONV_L2(l1,l2,fm07,fm06,fm05,fm04,fm03,fm02,fm01,fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07) \
	CONV_L2_INIT(lz,l1, 4,fm07,fm03      ,r0,r1); \
	CONV_L2_NEXT(lz,l1, 8,fm06,fm02,r0,r1,r0,r1); \
	CONV_L2_NEXT(lz,l1,12,fm05,fm01,r0,r1,r0,r1); \
	CONV_L2_NEXT(l1,l2, 0,fm04,fp00,r0,r1,r0,r1); \
	CONV_L2_NEXT(l1,l2, 4,fm03,fp01,r0,r1,r0,r1); \
	CONV_L2_NEXT(l1,l2, 8,fm02,fp02,r0,r1,r0,r1); \
	CONV_L2_NEXT(l1,l2,12,fm01,fp03,r0,r1,r0,r1); \
	CONV_L2_NEXT(l2,lz, 0,fp00,fp04,r0,r1,r0,r1); \
	CONV_L2_NEXT(l2,lz, 4,fp01,fp05,r0,r1,r0,r1); \
	CONV_L2_NEXT(l2,lz, 8,fp02,fp06,r0,r1,r0,r1); \
	CONV_L2_NEXT(l2,lz,12,fp03,fp07,r0,r1,r0,r1);

#define CONV_L2C(l1,l2,fm07,fm06,fm05,fm04,fm03,fm02,fm01,fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07) \
	CONV_L2_NEXT(lz,l1, 4,fm07,fm03,r0,r1,r0,r1); \
	CONV_L2_NEXT(lz,l1, 8,fm06,fm02,r0,r1,r0,r1); \
	CONV_L2_NEXT(lz,l1,12,fm05,fm01,r0,r1,r0,r1); \
	CONV_L2_NEXT(l1,l2, 0,fm04,fp00,r0,r1,r0,r1); \
	CONV_L2_NEXT(l1,l2, 4,fm03,fp01,r0,r1,r0,r1); \
	CONV_L2_NEXT(l1,l2, 8,fm02,fp02,r0,r1,r0,r1); \
	CONV_L2_NEXT(l1,l2,12,fm01,fp03,r0,r1,r0,r1); \
	CONV_L2_NEXT(l2,lz, 0,fp00,fp04,r0,r1,r0,r1); \
	CONV_L2_NEXT(l2,lz, 4,fp01,fp05,r0,r1,r0,r1); \
	CONV_L2_NEXT(l2,lz, 8,fp02,fp06,r0,r1,r0,r1); \
	CONV_L2_NEXT(l2,lz,12,fp03,fp07,r0,r1,r0,r1);

#define CONV_L3(l1,l2,l3,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01,fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11) \
	CONV_L3_INIT(lz,l1, 4,fm11,fm07,fm03         ,r0,r1,r2); \
	CONV_L3_NEXT(lz,l1, 8,fm10,fm06,fm02,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(lz,l1,12,fm09,fm05,fm01,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l1,l2, 0,fm08,fm04,fp00,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l1,l2, 4,fm07,fm03,fp01,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l1,l2, 8,fm06,fm02,fp02,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l1,l2,12,fm05,fm01,fp03,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l2,l3, 0,fm04,fp00,fp04,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l2,l3, 4,fm03,fp01,fp05,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l2,l3, 8,fm02,fp02,fp06,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l2,l3,12,fm01,fp03,fp07,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l3,lz, 0,fp00,fp04,fp08,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l3,lz, 4,fp01,fp05,fp09,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l3,lz, 8,fp02,fp06,fp10,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l3,lz,12,fp03,fp07,fp11,r0,r1,r2,r0,r1,r2);

#define CONV_L3C(l1,l2,l3,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01,fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11) \
	CONV_L3_NEXT(lz,l1, 4,fm11,fm07,fm03,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(lz,l1, 8,fm10,fm06,fm02,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(lz,l1,12,fm09,fm05,fm01,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l1,l2, 0,fm08,fm04,fp00,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l1,l2, 4,fm07,fm03,fp01,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l1,l2, 8,fm06,fm02,fp02,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l1,l2,12,fm05,fm01,fp03,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l2,l3, 0,fm04,fp00,fp04,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l2,l3, 4,fm03,fp01,fp05,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l2,l3, 8,fm02,fp02,fp06,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l2,l3,12,fm01,fp03,fp07,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l3,lz, 0,fp00,fp04,fp08,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l3,lz, 4,fp01,fp05,fp09,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l3,lz, 8,fp02,fp06,fp10,r0,r1,r2,r0,r1,r2); \
	CONV_L3_NEXT(l3,lz,12,fp03,fp07,fp11,r0,r1,r2,r0,r1,r2);

#define CONV_L4(l1,l2,l3,l4,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14) \
	CONV_L3_INIT(lz,l1, 4     ,fm11,fm07,fm03               ,r1,r2,r3); \
	CONV_L4_INIT(lz,l1, 8,fm14,fm10,fm06,fm02   ,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(lz,l1,12,fm13,fm09,fm05,fm01,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l1,l2, 0,fm12,fm08,fm04,fp00,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l1,l2, 4,fm11,fm07,fm03,fp01,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l1,l2, 8,fm10,fm06,fm02,fp02,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l1,l2,12,fm09,fm05,fm01,fp03,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l2,l3, 0,fm08,fm04,fp00,fp04,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l2,l3, 4,fm07,fm03,fp01,fp05,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l2,l3, 8,fm06,fm02,fp02,fp06,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l2,l3,12,fm05,fm01,fp03,fp07,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l3,l4, 0,fm04,fp00,fp04,fp08,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l3,l4, 4,fm03,fp01,fp05,fp09,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l3,l4, 8,fm02,fp02,fp06,fp10,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l3,l4,12,fm01,fp03,fp07,fp11,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l4,lz, 0,fp00,fp04,fp08,fp12,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l4,lz, 4,fp01,fp05,fp09,fp13,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l4,lz, 8,fp02,fp06,fp10,fp14,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L3_NEXT(l4,lz,12,fp03,fp07,fp11     ,r0,r1,r2   ,r0,r1,r2);

#define CONV_L4C(l1,l2,l3,l4,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14) \
	CONV_L3_NEXT(lz,l1, 4     ,fm11,fm07,fm03   ,r1,r2,r3   ,r1,r2,r3); \
	CONV_L4_NEXT(lz,l1, 8,fm14,fm10,fm06,fm02,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(lz,l1,12,fm13,fm09,fm05,fm01,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l1,l2, 0,fm12,fm08,fm04,fp00,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l1,l2, 4,fm11,fm07,fm03,fp01,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l1,l2, 8,fm10,fm06,fm02,fp02,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l1,l2,12,fm09,fm05,fm01,fp03,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l2,l3, 0,fm08,fm04,fp00,fp04,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l2,l3, 4,fm07,fm03,fp01,fp05,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l2,l3, 8,fm06,fm02,fp02,fp06,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l2,l3,12,fm05,fm01,fp03,fp07,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l3,l4, 0,fm04,fp00,fp04,fp08,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l3,l4, 4,fm03,fp01,fp05,fp09,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l3,l4, 8,fm02,fp02,fp06,fp10,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l3,l4,12,fm01,fp03,fp07,fp11,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l4,lz, 0,fp00,fp04,fp08,fp12,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l4,lz, 4,fp01,fp05,fp09,fp13,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L4_NEXT(l4,lz, 8,fp02,fp06,fp10,fp14,r0,r1,r2,r3,r0,r1,r2,r3); \
	CONV_L3_NEXT(l4,lz,12,fp03,fp07,fp11     ,r0,r1,r2   ,r0,r1,r2);

#define CONV_L5(l1,l2,l3,l4,l5,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14) \
	CONV_L3_INIT(lz,l1, 4          ,fm11,fm07,fm03                     ,r2,r3,r4); \
	CONV_L4_INIT(lz,l1, 8     ,fm14,fm10,fm06,fm02      ,r2,r3,r4   ,r1,r2,r3,r4); \
	CONV_L4_NEXT(lz,l1,12     ,fm13,fm09,fm05,fm01   ,r1,r2,r3,r4   ,r1,r2,r3,r4); \
	CONV_L4_NEXT(l1,l2, 0     ,fm12,fm08,fm04,fp00   ,r1,r2,r3,r4   ,r1,r2,r3,r4); \
	CONV_L4_NEXT(l1,l2, 4     ,fm11,fm07,fm03,fp01   ,r1,r2,r3,r4   ,r1,r2,r3,r4); \
	CONV_L5_INIT(l1,l2, 8,fm14,fm10,fm06,fm02,fp02   ,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l1,l2,12,fm13,fm09,fm05,fm01,fp03,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l2,l3, 0,fm12,fm08,fm04,fp00,fp04,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l2,l3, 4,fm11,fm07,fm03,fp01,fp05,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l2,l3, 8,fm10,fm06,fm02,fp02,fp06,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l2,l3,12,fm09,fm05,fm01,fp03,fp07,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l3,l4, 0,fm08,fm04,fp00,fp04,fp08,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l3,l4, 4,fm07,fm03,fp01,fp05,fp09,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l3,l4, 8,fm06,fm02,fp02,fp06,fp10,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l3,l4,12,fm05,fm01,fp03,fp07,fp11,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l4,l5, 0,fm04,fp00,fp04,fp08,fp12,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l4,l5, 4,fm03,fp01,fp05,fp09,fp13,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l4,l5, 8,fm02,fp02,fp06,fp10,fp14,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L4_NEXT(l4,l5,12,fm01,fp03,fp07,fp11     ,r0,r1,r2,r3   ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l5,lz, 0,fp00,fp04,fp08,fp12     ,r0,r1,r2,r3   ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l5,lz, 4,fp01,fp05,fp09,fp13     ,r0,r1,r2,r3   ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l5,lz, 8,fp02,fp06,fp10,fp14     ,r0,r1,r2,r3   ,r0,r1,r2,r3); \
	CONV_L3_NEXT(l5,lz,12,fp03,fp07,fp11          ,r0,r1,r2      ,r0,r1,r2);

#define CONV_L5C(l1,l2,l3,l4,l5,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14) \
	CONV_L3_NEXT(lz,l1, 4          ,fm11,fm07,fm03      ,r2,r3,r4      ,r2,r3,r4); \
	CONV_L4_NEXT(lz,l1, 8     ,fm14,fm10,fm06,fm02   ,r1,r2,r3,r4   ,r1,r2,r3,r4); \
	CONV_L4_NEXT(lz,l1,12     ,fm13,fm09,fm05,fm01   ,r1,r2,r3,r4   ,r1,r2,r3,r4); \
	CONV_L4_NEXT(l1,l2, 0     ,fm12,fm08,fm04,fp00   ,r1,r2,r3,r4   ,r1,r2,r3,r4); \
	CONV_L4_NEXT(l1,l2, 4     ,fm11,fm07,fm03,fp01   ,r1,r2,r3,r4   ,r1,r2,r3,r4); \
	CONV_L5_NEXT(l1,l2, 8,fm14,fm10,fm06,fm02,fp02,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l1,l2,12,fm13,fm09,fm05,fm01,fp03,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l2,l3, 0,fm12,fm08,fm04,fp00,fp04,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l2,l3, 4,fm11,fm07,fm03,fp01,fp05,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l2,l3, 8,fm10,fm06,fm02,fp02,fp06,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l2,l3,12,fm09,fm05,fm01,fp03,fp07,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l3,l4, 0,fm08,fm04,fp00,fp04,fp08,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l3,l4, 4,fm07,fm03,fp01,fp05,fp09,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l3,l4, 8,fm06,fm02,fp02,fp06,fp10,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l3,l4,12,fm05,fm01,fp03,fp07,fp11,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l4,l5, 0,fm04,fp00,fp04,fp08,fp12,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l4,l5, 4,fm03,fp01,fp05,fp09,fp13,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l4,l5, 8,fm02,fp02,fp06,fp10,fp14,r0,r1,r2,r3,r4,r0,r1,r2,r3,r4); \
	CONV_L4_NEXT(l4,l5,12,fm01,fp03,fp07,fp11     ,r0,r1,r2,r3   ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l5,lz, 0,fp00,fp04,fp08,fp12     ,r0,r1,r2,r3   ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l5,lz, 4,fp01,fp05,fp09,fp13     ,r0,r1,r2,r3   ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l5,lz, 8,fp02,fp06,fp10,fp14     ,r0,r1,r2,r3   ,r0,r1,r2,r3); \
	CONV_L3_NEXT(l5,lz,12,fp03,fp07,fp11          ,r0,r1,r2      ,r0,r1,r2);

#define CONV_L6(l1,l2,l3,l4,l5,l6,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14) \
	CONV_L3_INIT(lz,l1, 4               ,fm11,fm07,fm03                           ,r3,r4,r5); \
	CONV_L4_INIT(lz,l1, 8          ,fm14,fm10,fm06,fm02         ,r3,r4,r5      ,r2,r3,r4,r5); \
	CONV_L4_NEXT(lz,l1,12          ,fm13,fm09,fm05,fm01      ,r2,r3,r4,r5      ,r2,r3,r4,r5); \
	CONV_L4_NEXT(l1,l2, 0          ,fm12,fm08,fm04,fp00      ,r2,r3,r4,r5      ,r2,r3,r4,r5); \
	CONV_L4_NEXT(l1,l2, 4          ,fm11,fm07,fm03,fp01      ,r2,r3,r4,r5      ,r2,r3,r4,r5); \
	CONV_L5_INIT(l1,l2, 8     ,fm14,fm10,fm06,fm02,fp02      ,r2,r3,r4,r5   ,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l1,l2,12     ,fm13,fm09,fm05,fm01,fp03   ,r1,r2,r3,r4,r5   ,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l2,l3, 0     ,fm12,fm08,fm04,fp00,fp04   ,r1,r2,r3,r4,r5   ,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l2,l3, 4     ,fm11,fm07,fm03,fp01,fp05   ,r1,r2,r3,r4,r5   ,r1,r2,r3,r4,r5); \
	CONV_L6_INIT(l2,l3, 8,fm14,fm10,fm06,fm02,fp02,fp06   ,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l2,l3,12,fm13,fm09,fm05,fm01,fp03,fp07,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l3,l4, 0,fm12,fm08,fm04,fp00,fp04,fp08,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l3,l4, 4,fm11,fm07,fm03,fp01,fp05,fp09,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l3,l4, 8,fm10,fm06,fm02,fp02,fp06,fp10,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l3,l4,12,fm09,fm05,fm01,fp03,fp07,fp11,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l4,l5, 0,fm08,fm04,fp00,fp04,fp08,fp12,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l4,l5, 4,fm07,fm03,fp01,fp05,fp09,fp13,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l4,l5, 8,fm06,fm02,fp02,fp06,fp10,fp14,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l4,l5,12,fm05,fm01,fp03,fp07,fp11     ,r0,r1,r2,r3,r4   ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l5,l6, 0,fm04,fp00,fp04,fp08,fp12     ,r0,r1,r2,r3,r4   ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l5,l6, 4,fm03,fp01,fp05,fp09,fp13     ,r0,r1,r2,r3,r4   ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l5,l6, 8,fm02,fp02,fp06,fp10,fp14     ,r0,r1,r2,r3,r4   ,r0,r1,r2,r3,r4); \
	CONV_L4_NEXT(l5,l6,12,fm01,fp03,fp07,fp11          ,r0,r1,r2,r3      ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l6,lz, 0,fp00,fp04,fp08,fp12          ,r0,r1,r2,r3      ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l6,lz, 4,fp01,fp05,fp09,fp13          ,r0,r1,r2,r3      ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l6,lz, 8,fp02,fp06,fp10,fp14          ,r0,r1,r2,r3      ,r0,r1,r2,r3); \
	CONV_L3_NEXT(l6,lz,12,fp03,fp07,fp11               ,r0,r1,r2         ,r0,r1,r2);

#define CONV_L6C(l1,l2,l3,l4,l5,l6,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14) \
	CONV_L3_NEXT(lz,l1, 4               ,fm11,fm07,fm03         ,r3,r4,r5         ,r3,r4,r5); \
	CONV_L4_NEXT(lz,l1, 8          ,fm14,fm10,fm06,fm02      ,r2,r3,r4,r5      ,r2,r3,r4,r5); \
	CONV_L4_NEXT(lz,l1,12          ,fm13,fm09,fm05,fm01      ,r2,r3,r4,r5      ,r2,r3,r4,r5); \
	CONV_L4_NEXT(l1,l2, 0          ,fm12,fm08,fm04,fp00      ,r2,r3,r4,r5      ,r2,r3,r4,r5); \
	CONV_L4_NEXT(l1,l2, 4          ,fm11,fm07,fm03,fp01      ,r2,r3,r4,r5      ,r2,r3,r4,r5); \
	CONV_L5_NEXT(l1,l2, 8     ,fm14,fm10,fm06,fm02,fp02   ,r1,r2,r3,r4,r5   ,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l1,l2,12     ,fm13,fm09,fm05,fm01,fp03   ,r1,r2,r3,r4,r5   ,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l2,l3, 0     ,fm12,fm08,fm04,fp00,fp04   ,r1,r2,r3,r4,r5   ,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l2,l3, 4     ,fm11,fm07,fm03,fp01,fp05   ,r1,r2,r3,r4,r5   ,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l2,l3, 8,fm14,fm10,fm06,fm02,fp02,fp06,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l2,l3,12,fm13,fm09,fm05,fm01,fp03,fp07,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l3,l4, 0,fm12,fm08,fm04,fp00,fp04,fp08,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l3,l4, 4,fm11,fm07,fm03,fp01,fp05,fp09,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l3,l4, 8,fm10,fm06,fm02,fp02,fp06,fp10,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l3,l4,12,fm09,fm05,fm01,fp03,fp07,fp11,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l4,l5, 0,fm08,fm04,fp00,fp04,fp08,fp12,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l4,l5, 4,fm07,fm03,fp01,fp05,fp09,fp13,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l4,l5, 8,fm06,fm02,fp02,fp06,fp10,fp14,r0,r1,r2,r3,r4,r5,r0,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l4,l5,12,fm05,fm01,fp03,fp07,fp11     ,r0,r1,r2,r3,r4   ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l5,l6, 0,fm04,fp00,fp04,fp08,fp12     ,r0,r1,r2,r3,r4   ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l5,l6, 4,fm03,fp01,fp05,fp09,fp13     ,r0,r1,r2,r3,r4   ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l5,l6, 8,fm02,fp02,fp06,fp10,fp14     ,r0,r1,r2,r3,r4   ,r0,r1,r2,r3,r4); \
	CONV_L4_NEXT(l5,l6,12,fm01,fp03,fp07,fp11          ,r0,r1,r2,r3      ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l6,lz, 0,fp00,fp04,fp08,fp12          ,r0,r1,r2,r3      ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l6,lz, 4,fp01,fp05,fp09,fp13          ,r0,r1,r2,r3      ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l6,lz, 8,fp02,fp06,fp10,fp14          ,r0,r1,r2,r3      ,r0,r1,r2,r3); \
	CONV_L3_NEXT(l6,lz,12,fp03,fp07,fp11               ,r0,r1,r2         ,r0,r1,r2);

#define CONV_L7(l1,l2,l3,l4,l5,l6,l7,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14) \
	CONV_L3_INIT(lz,l1, 4                    ,fm11,fm07,fm03                                  ,r4,r5,r6); \
	CONV_L4_INIT(lz,l1, 8               ,fm14,fm10,fm06,fm02             ,r4,r5,r6         ,r3,r4,r5,r6); \
	CONV_L4_NEXT(lz,l1,12               ,fm13,fm09,fm05,fm01          ,r3,r4,r5,r6         ,r3,r4,r5,r6); \
	CONV_L4_NEXT(l1,l2, 0               ,fm12,fm08,fm04,fp00          ,r3,r4,r5,r6         ,r3,r4,r5,r6); \
	CONV_L4_NEXT(l1,l2, 4               ,fm11,fm07,fm03,fp01          ,r3,r4,r5,r6         ,r3,r4,r5,r6); \
	CONV_L5_INIT(l1,l2, 8          ,fm14,fm10,fm06,fm02,fp02          ,r3,r4,r5,r6      ,r2,r3,r4,r5,r6); \
	CONV_L5_NEXT(l1,l2,12          ,fm13,fm09,fm05,fm01,fp03       ,r2,r3,r4,r5,r6      ,r2,r3,r4,r5,r6); \
	CONV_L5_NEXT(l2,l3, 0          ,fm12,fm08,fm04,fp00,fp04       ,r2,r3,r4,r5,r6      ,r2,r3,r4,r5,r6); \
	CONV_L5_NEXT(l2,l3, 4          ,fm11,fm07,fm03,fp01,fp05       ,r2,r3,r4,r5,r6      ,r2,r3,r4,r5,r6); \
	CONV_L6_INIT(l2,l3, 8     ,fm14,fm10,fm06,fm02,fp02,fp06       ,r2,r3,r4,r5,r6   ,r1,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l2,l3,12     ,fm13,fm09,fm05,fm01,fp03,fp07    ,r1,r2,r3,r4,r5,r6   ,r1,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l3,l4, 0     ,fm12,fm08,fm04,fp00,fp04,fp08    ,r1,r2,r3,r4,r5,r6   ,r1,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l3,l4, 4     ,fm11,fm07,fm03,fp01,fp05,fp09    ,r1,r2,r3,r4,r5,r6   ,r1,r2,r3,r4,r5,r6); \
	CONV_L7_INIT(l3,l4, 8,fm14,fm10,fm06,fm02,fp02,fp06,fp10   ,r1,r2,r3,r4,r5,r6,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l3,l4,12,fm13,fm09,fm05,fm01,fp03,fp07,fp11,r0,r1,r2,r3,r4,r5,r6,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l4,l5, 0,fm12,fm08,fm04,fp00,fp04,fp08,fp12,r0,r1,r2,r3,r4,r5,r6,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l4,l5, 4,fm11,fm07,fm03,fp01,fp05,fp09,fp13,r0,r1,r2,r3,r4,r5,r6,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l4,l5, 8,fm10,fm06,fm02,fp02,fp06,fp10,fp14,r0,r1,r2,r3,r4,r5,r6,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l4,l5,12,fm09,fm05,fm01,fp03,fp07,fp11     ,r0,r1,r2,r3,r4,r5   ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l5,l6, 0,fm08,fm04,fp00,fp04,fp08,fp12     ,r0,r1,r2,r3,r4,r5   ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l5,l6, 4,fm07,fm03,fp01,fp05,fp09,fp13     ,r0,r1,r2,r3,r4,r5   ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l5,l6, 8,fm06,fm02,fp02,fp06,fp10,fp14     ,r0,r1,r2,r3,r4,r5   ,r0,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l5,l6,12,fm05,fm01,fp03,fp07,fp11          ,r0,r1,r2,r3,r4      ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l6,l7, 0,fm04,fp00,fp04,fp08,fp12          ,r0,r1,r2,r3,r4      ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l6,l7, 4,fm03,fp01,fp05,fp09,fp13          ,r0,r1,r2,r3,r4      ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l6,l7, 8,fm02,fp02,fp06,fp10,fp14          ,r0,r1,r2,r3,r4      ,r0,r1,r2,r3,r4); \
	CONV_L4_NEXT(l6,l7,12,fm01,fp03,fp07,fp11               ,r0,r1,r2,r3         ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l7,lz, 0,fp00,fp04,fp08,fp12               ,r0,r1,r2,r3         ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l7,lz, 4,fp01,fp05,fp09,fp13               ,r0,r1,r2,r3         ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l7,lz, 8,fp02,fp06,fp10,fp14               ,r0,r1,r2,r3         ,r0,r1,r2,r3); \
	CONV_L3_NEXT(l7,lz,12,fp03,fp07,fp11                    ,r0,r1,r2            ,r0,r1,r2);

#define CONV_L7C(l1,l2,l3,l4,l5,l6,l7,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14) \
	CONV_L3_NEXT(lz,l1, 4                    ,fm11,fm07,fm03             ,r4,r5,r6            ,r4,r5,r6); \
	CONV_L4_NEXT(lz,l1, 8               ,fm14,fm10,fm06,fm02          ,r3,r4,r5,r6         ,r3,r4,r5,r6); \
	CONV_L4_NEXT(lz,l1,12               ,fm13,fm09,fm05,fm01          ,r3,r4,r5,r6         ,r3,r4,r5,r6); \
	CONV_L4_NEXT(l1,l2, 0               ,fm12,fm08,fm04,fp00          ,r3,r4,r5,r6         ,r3,r4,r5,r6); \
	CONV_L4_NEXT(l1,l2, 4               ,fm11,fm07,fm03,fp01          ,r3,r4,r5,r6         ,r3,r4,r5,r6); \
	CONV_L5_NEXT(l1,l2, 8          ,fm14,fm10,fm06,fm02,fp02       ,r2,r3,r4,r5,r6      ,r2,r3,r4,r5,r6); \
	CONV_L5_NEXT(l1,l2,12          ,fm13,fm09,fm05,fm01,fp03       ,r2,r3,r4,r5,r6      ,r2,r3,r4,r5,r6); \
	CONV_L5_NEXT(l2,l3, 0          ,fm12,fm08,fm04,fp00,fp04       ,r2,r3,r4,r5,r6      ,r2,r3,r4,r5,r6); \
	CONV_L5_NEXT(l2,l3, 4          ,fm11,fm07,fm03,fp01,fp05       ,r2,r3,r4,r5,r6      ,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l2,l3, 8     ,fm14,fm10,fm06,fm02,fp02,fp06    ,r1,r2,r3,r4,r5,r6   ,r1,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l2,l3,12     ,fm13,fm09,fm05,fm01,fp03,fp07    ,r1,r2,r3,r4,r5,r6   ,r1,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l3,l4, 0     ,fm12,fm08,fm04,fp00,fp04,fp08    ,r1,r2,r3,r4,r5,r6   ,r1,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l3,l4, 4     ,fm11,fm07,fm03,fp01,fp05,fp09    ,r1,r2,r3,r4,r5,r6   ,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l3,l4, 8,fm14,fm10,fm06,fm02,fp02,fp06,fp10,r0,r1,r2,r3,r4,r5,r6,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l3,l4,12,fm13,fm09,fm05,fm01,fp03,fp07,fp11,r0,r1,r2,r3,r4,r5,r6,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l4,l5, 0,fm12,fm08,fm04,fp00,fp04,fp08,fp12,r0,r1,r2,r3,r4,r5,r6,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l4,l5, 4,fm11,fm07,fm03,fp01,fp05,fp09,fp13,r0,r1,r2,r3,r4,r5,r6,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l4,l5, 8,fm10,fm06,fm02,fp02,fp06,fp10,fp14,r0,r1,r2,r3,r4,r5,r6,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l4,l5,12,fm09,fm05,fm01,fp03,fp07,fp11     ,r0,r1,r2,r3,r4,r5   ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l5,l6, 0,fm08,fm04,fp00,fp04,fp08,fp12     ,r0,r1,r2,r3,r4,r5   ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l5,l6, 4,fm07,fm03,fp01,fp05,fp09,fp13     ,r0,r1,r2,r3,r4,r5   ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l5,l6, 8,fm06,fm02,fp02,fp06,fp10,fp14     ,r0,r1,r2,r3,r4,r5   ,r0,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l5,l6,12,fm05,fm01,fp03,fp07,fp11          ,r0,r1,r2,r3,r4      ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l6,l7, 0,fm04,fp00,fp04,fp08,fp12          ,r0,r1,r2,r3,r4      ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l6,l7, 4,fm03,fp01,fp05,fp09,fp13          ,r0,r1,r2,r3,r4      ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l6,l7, 8,fm02,fp02,fp06,fp10,fp14          ,r0,r1,r2,r3,r4      ,r0,r1,r2,r3,r4); \
	CONV_L4_NEXT(l6,l7,12,fm01,fp03,fp07,fp11               ,r0,r1,r2,r3         ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l7,lz, 0,fp00,fp04,fp08,fp12               ,r0,r1,r2,r3         ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l7,lz, 4,fp01,fp05,fp09,fp13               ,r0,r1,r2,r3         ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l7,lz, 8,fp02,fp06,fp10,fp14               ,r0,r1,r2,r3         ,r0,r1,r2,r3); \
	CONV_L3_NEXT(l7,lz,12,fp03,fp07,fp11                    ,r0,r1,r2            ,r0,r1,r2);

#define CONV_L8_PART1(l1,l2,l3,l4,l5,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,r1,r2,r3,r4,r5,r6,r7) \
	CONV_L3_INIT(lz,l1, 4                    ,fm11,fm07,fm03                                 ,r5,r6,r7); \
	CONV_L4_INIT(lz,l1, 8               ,fm14,fm10,fm06,fm02            ,r5,r6,r7         ,r4,r5,r6,r7); \
	CONV_L4_NEXT(lz,l1,12               ,fm13,fm09,fm05,fm01         ,r4,r5,r6,r7         ,r4,r5,r6,r7); \
	CONV_L4_NEXT(l1,l2, 0               ,fm12,fm08,fm04,fp00         ,r4,r5,r6,r7         ,r4,r5,r6,r7); \
	CONV_L4_NEXT(l1,l2, 4               ,fm11,fm07,fm03,fp01         ,r4,r5,r6,r7         ,r4,r5,r6,r7); \
	CONV_L5_INIT(l1,l2, 8          ,fm14,fm10,fm06,fm02,fp02         ,r4,r5,r6,r7      ,r3,r4,r5,r6,r7); \
	CONV_L5_NEXT(l1,l2,12          ,fm13,fm09,fm05,fm01,fp03      ,r3,r4,r5,r6,r7      ,r3,r4,r5,r6,r7); \
	CONV_L5_NEXT(l2,l3, 0          ,fm12,fm08,fm04,fp00,fp04      ,r3,r4,r5,r6,r7      ,r3,r4,r5,r6,r7); \
	CONV_L5_NEXT(l2,l3, 4          ,fm11,fm07,fm03,fp01,fp05      ,r3,r4,r5,r6,r7      ,r3,r4,r5,r6,r7); \
	CONV_L6_INIT(l2,l3, 8     ,fm14,fm10,fm06,fm02,fp02,fp06      ,r3,r4,r5,r6,r7   ,r2,r3,r4,r5,r6,r7); \
	CONV_L6_NEXT(l2,l3,12     ,fm13,fm09,fm05,fm01,fp03,fp07   ,r2,r3,r4,r5,r6,r7   ,r2,r3,r4,r5,r6,r7); \
	CONV_L6_NEXT(l3,l4, 0     ,fm12,fm08,fm04,fp00,fp04,fp08   ,r2,r3,r4,r5,r6,r7   ,r2,r3,r4,r5,r6,r7); \
	CONV_L6_NEXT(l3,l4, 4     ,fm11,fm07,fm03,fp01,fp05,fp09   ,r2,r3,r4,r5,r6,r7   ,r2,r3,r4,r5,r6,r7); \
	CONV_L7_INIT(l3,l4, 8,fm14,fm10,fm06,fm02,fp02,fp06,fp10   ,r2,r3,r4,r5,r6,r7,r1,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l3,l4,12,fm13,fm09,fm05,fm01,fp03,fp07,fp11,r1,r2,r3,r4,r5,r6,r7,r1,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l4,l5, 0,fm12,fm08,fm04,fp00,fp04,fp08,fp12,r1,r2,r3,r4,r5,r6,r7,r1,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l4,l5, 4,fm11,fm07,fm03,fp01,fp05,fp09,fp13,r1,r2,r3,r4,r5,r6,r7,r1,r2,r3,r4,r5,r6,r7);

#define CONV_L8_PART1C(l1,l2,l3,l4,l5,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,r1,r2,r3,r4,r5,r6,r7) \
	CONV_L3_NEXT(lz,l1, 4                    ,fm11,fm07,fm03            ,r5,r6,r7            ,r5,r6,r7); \
	CONV_L4_NEXT(lz,l1, 8               ,fm14,fm10,fm06,fm02         ,r4,r5,r6,r7         ,r4,r5,r6,r7); \
	CONV_L4_NEXT(lz,l1,12               ,fm13,fm09,fm05,fm01         ,r4,r5,r6,r7         ,r4,r5,r6,r7); \
	CONV_L4_NEXT(l1,l2, 0               ,fm12,fm08,fm04,fp00         ,r4,r5,r6,r7         ,r4,r5,r6,r7); \
	CONV_L4_NEXT(l1,l2, 4               ,fm11,fm07,fm03,fp01         ,r4,r5,r6,r7         ,r4,r5,r6,r7); \
	CONV_L5_NEXT(l1,l2, 8          ,fm14,fm10,fm06,fm02,fp02      ,r3,r4,r5,r6,r7      ,r3,r4,r5,r6,r7); \
	CONV_L5_NEXT(l1,l2,12          ,fm13,fm09,fm05,fm01,fp03      ,r3,r4,r5,r6,r7      ,r3,r4,r5,r6,r7); \
	CONV_L5_NEXT(l2,l3, 0          ,fm12,fm08,fm04,fp00,fp04      ,r3,r4,r5,r6,r7      ,r3,r4,r5,r6,r7); \
	CONV_L5_NEXT(l2,l3, 4          ,fm11,fm07,fm03,fp01,fp05      ,r3,r4,r5,r6,r7      ,r3,r4,r5,r6,r7); \
	CONV_L6_NEXT(l2,l3, 8     ,fm14,fm10,fm06,fm02,fp02,fp06   ,r2,r3,r4,r5,r6,r7   ,r2,r3,r4,r5,r6,r7); \
	CONV_L6_NEXT(l2,l3,12     ,fm13,fm09,fm05,fm01,fp03,fp07   ,r2,r3,r4,r5,r6,r7   ,r2,r3,r4,r5,r6,r7); \
	CONV_L6_NEXT(l3,l4, 0     ,fm12,fm08,fm04,fp00,fp04,fp08   ,r2,r3,r4,r5,r6,r7   ,r2,r3,r4,r5,r6,r7); \
	CONV_L6_NEXT(l3,l4, 4     ,fm11,fm07,fm03,fp01,fp05,fp09   ,r2,r3,r4,r5,r6,r7   ,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l3,l4, 8,fm14,fm10,fm06,fm02,fp02,fp06,fp10,r1,r2,r3,r4,r5,r6,r7,r1,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l3,l4,12,fm13,fm09,fm05,fm01,fp03,fp07,fp11,r1,r2,r3,r4,r5,r6,r7,r1,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l4,l5, 0,fm12,fm08,fm04,fp00,fp04,fp08,fp12,r1,r2,r3,r4,r5,r6,r7,r1,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l4,l5, 4,fm11,fm07,fm03,fp01,fp05,fp09,fp13,r1,r2,r3,r4,r5,r6,r7,r1,r2,r3,r4,r5,r6,r7);

#define CONV_L8_PART2(l4,l5,l6,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14,r0,r1,r2,r3,r4,r5,r6,r7) \
	CONV_L8_INIT(l4,l5, 8,fm14,fm10,fm06,fm02,fp02,fp06,fp10,fp14   ,r1,r2,r3,r4,r5,r6,r7,r0,r1,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l4,l5,12,fm13,fm09,fm05,fm01,fp03,fp07,fp11     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l5,l6, 0,fm12,fm08,fm04,fp00,fp04,fp08,fp12     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l5,l6, 4,fm11,fm07,fm03,fp01,fp05,fp09,fp13     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6);

#define CONV_L8_PART2C(l4,l5,l6,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14,r0,r1,r2,r3,r4,r5,r6,r7) \
	CONV_L8_NEXT(l4,l5, 8,fm14,fm10,fm06,fm02,fp02,fp06,fp10,fp14,r0,r1,r2,r3,r4,r5,r6,r7,r0,r1,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l4,l5,12,fm13,fm09,fm05,fm01,fp03,fp07,fp11     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l5,l6, 0,fm12,fm08,fm04,fp00,fp04,fp08,fp12     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l5,l6, 4,fm11,fm07,fm03,fp01,fp05,fp09,fp13     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6);

#define CONV_L8_PART3(l4,l5,l6,l7,l8,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14,r0,r1,r2,r3,r4,r5,r6,r7) \
	CONV_L8_INIT(l4,l5, 8,fm14,fm10,fm06,fm02,fp02,fp06,fp10,fp14   ,r1,r2,r3,r4,r5,r6,r7,r0,r1,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l4,l5,12,fm13,fm09,fm05,fm01,fp03,fp07,fp11     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l5,l6, 0,fm12,fm08,fm04,fp00,fp04,fp08,fp12     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l5,l6, 4,fm11,fm07,fm03,fp01,fp05,fp09,fp13     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l5,l6, 8,fm10,fm06,fm02,fp02,fp06,fp10,fp14     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l5,l6,12,fm09,fm05,fm01,fp03,fp07,fp11          ,r0,r1,r2,r3,r4,r5      ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l6,l7, 0,fm08,fm04,fp00,fp04,fp08,fp12          ,r0,r1,r2,r3,r4,r5      ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l6,l7, 4,fm07,fm03,fp01,fp05,fp09,fp13          ,r0,r1,r2,r3,r4,r5      ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l6,l7, 8,fm06,fm02,fp02,fp06,fp10,fp14          ,r0,r1,r2,r3,r4,r5      ,r0,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l6,l7,12,fm05,fm01,fp03,fp07,fp11               ,r0,r1,r2,r3,r4         ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l7,l8, 0,fm04,fp00,fp04,fp08,fp12               ,r0,r1,r2,r3,r4         ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l7,l8, 4,fm03,fp01,fp05,fp09,fp13               ,r0,r1,r2,r3,r4         ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l7,l8, 8,fm02,fp02,fp06,fp10,fp14               ,r0,r1,r2,r3,r4         ,r0,r1,r2,r3,r4); \
	CONV_L4_NEXT(l7,l8,12,fm01,fp03,fp07,fp11                    ,r0,r1,r2,r3            ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l8,lz, 0,fp00,fp04,fp08,fp12                    ,r0,r1,r2,r3            ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l8,lz, 4,fp01,fp05,fp09,fp13                    ,r0,r1,r2,r3            ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l8,lz, 8,fp02,fp06,fp10,fp14                    ,r0,r1,r2,r3            ,r0,r1,r2,r3); \
	CONV_L3_NEXT(l8,lz,12,fp03,fp07,fp11                         ,r0,r1,r2               ,r0,r1,r2);

#define CONV_L8_PART3C(l4,l5,l6,l7,l8,fm14,fm13,fm12,fm11,fm10,fm09,fm08,fm07,fm06,fm05,fm04,fm03,fm02,fm01, \
		fp00,fp01,fp02,fp03,fp04,fp05,fp06,fp07,fp08,fp09,fp10,fp11,fp12,fp13,fp14,r0,r1,r2,r3,r4,r5,r6,r7) \
	CONV_L8_NEXT(l4,l5, 8,fm14,fm10,fm06,fm02,fp02,fp06,fp10,fp14,r0,r1,r2,r3,r4,r5,r6,r7,r0,r1,r2,r3,r4,r5,r6,r7); \
	CONV_L7_NEXT(l4,l5,12,fm13,fm09,fm05,fm01,fp03,fp07,fp11     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l5,l6, 0,fm12,fm08,fm04,fp00,fp04,fp08,fp12     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l5,l6, 4,fm11,fm07,fm03,fp01,fp05,fp09,fp13     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L7_NEXT(l5,l6, 8,fm10,fm06,fm02,fp02,fp06,fp10,fp14     ,r0,r1,r2,r3,r4,r5,r6   ,r0,r1,r2,r3,r4,r5,r6); \
	CONV_L6_NEXT(l5,l6,12,fm09,fm05,fm01,fp03,fp07,fp11          ,r0,r1,r2,r3,r4,r5      ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l6,l7, 0,fm08,fm04,fp00,fp04,fp08,fp12          ,r0,r1,r2,r3,r4,r5      ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l6,l7, 4,fm07,fm03,fp01,fp05,fp09,fp13          ,r0,r1,r2,r3,r4,r5      ,r0,r1,r2,r3,r4,r5); \
	CONV_L6_NEXT(l6,l7, 8,fm06,fm02,fp02,fp06,fp10,fp14          ,r0,r1,r2,r3,r4,r5      ,r0,r1,r2,r3,r4,r5); \
	CONV_L5_NEXT(l6,l7,12,fm05,fm01,fp03,fp07,fp11               ,r0,r1,r2,r3,r4         ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l7,l8, 0,fm04,fp00,fp04,fp08,fp12               ,r0,r1,r2,r3,r4         ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l7,l8, 4,fm03,fp01,fp05,fp09,fp13               ,r0,r1,r2,r3,r4         ,r0,r1,r2,r3,r4); \
	CONV_L5_NEXT(l7,l8, 8,fm02,fp02,fp06,fp10,fp14               ,r0,r1,r2,r3,r4         ,r0,r1,r2,r3,r4); \
	CONV_L4_NEXT(l7,l8,12,fm01,fp03,fp07,fp11                    ,r0,r1,r2,r3            ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l8,lz, 0,fp00,fp04,fp08,fp12                    ,r0,r1,r2,r3            ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l8,lz, 4,fp01,fp05,fp09,fp13                    ,r0,r1,r2,r3            ,r0,r1,r2,r3); \
	CONV_L4_NEXT(l8,lz, 8,fp02,fp06,fp10,fp14                    ,r0,r1,r2,r3            ,r0,r1,r2,r3); \
	CONV_L3_NEXT(l8,lz,12,fp03,fp07,fp11                         ,r0,r1,r2               ,r0,r1,r2);

#define FILTER_A \
	am14=_mm_set1_ps(fa[ 0]); \
	am13=_mm_set1_ps(fa[ 1]); \
	am12=_mm_set1_ps(fa[ 2]); \
	am11=_mm_set1_ps(fa[ 3]); \
	am10=_mm_set1_ps(fa[ 4]); \
	am09=_mm_set1_ps(fa[ 5]); \
	am08=_mm_set1_ps(fa[ 6]); \
	am07=_mm_set1_ps(fa[ 7]); \
	am06=_mm_set1_ps(fa[ 8]); \
	am05=_mm_set1_ps(fa[ 9]); \
	am04=_mm_set1_ps(fa[10]); \
	am03=_mm_set1_ps(fa[11]); \
	am02=_mm_set1_ps(fa[12]); \
	am01=_mm_set1_ps(fa[13]); \
	ap00=_mm_set1_ps(fa[14]); \
	ap01=_mm_set1_ps(fa[15]); \
	ap02=_mm_set1_ps(fa[16]); \
	ap03=_mm_set1_ps(fa[17]); \
	ap04=_mm_set1_ps(fa[18]); \
	ap05=_mm_set1_ps(fa[19]); \
	ap06=_mm_set1_ps(fa[20]); \
	ap07=_mm_set1_ps(fa[21]); \
	ap08=_mm_set1_ps(fa[22]); \
	ap09=_mm_set1_ps(fa[23]); \
	ap10=_mm_set1_ps(fa[24]); \
	ap11=_mm_set1_ps(fa[25]); \
	ap12=_mm_set1_ps(fa[26]); \
	ap13=_mm_set1_ps(fa[27]); \
	ap14=_mm_set1_ps(fa[28]);

#define FILTER_B \
	bm14=_mm_set1_ps(fb[ 0]); \
	bm13=_mm_set1_ps(fb[ 1]); \
	bm12=_mm_set1_ps(fb[ 2]); \
	bm11=_mm_set1_ps(fb[ 3]); \
	bm10=_mm_set1_ps(fb[ 4]); \
	bm09=_mm_set1_ps(fb[ 5]); \
	bm08=_mm_set1_ps(fb[ 6]); \
	bm07=_mm_set1_ps(fb[ 7]); \
	bm06=_mm_set1_ps(fb[ 8]); \
	bm05=_mm_set1_ps(fb[ 9]); \
	bm04=_mm_set1_ps(fb[10]); \
	bm03=_mm_set1_ps(fb[11]); \
	bm02=_mm_set1_ps(fb[12]); \
	bm01=_mm_set1_ps(fb[13]); \
	bp00=_mm_set1_ps(fb[14]); \
	bp01=_mm_set1_ps(fb[15]); \
	bp02=_mm_set1_ps(fb[16]); \
	bp03=_mm_set1_ps(fb[17]); \
	bp04=_mm_set1_ps(fb[18]); \
	bp05=_mm_set1_ps(fb[19]); \
	bp06=_mm_set1_ps(fb[20]); \
	bp07=_mm_set1_ps(fb[21]); \
	bp08=_mm_set1_ps(fb[22]); \
	bp09=_mm_set1_ps(fb[23]); \
	bp10=_mm_set1_ps(fb[24]); \
	bp11=_mm_set1_ps(fb[25]); \
	bp12=_mm_set1_ps(fb[26]); \
	bp13=_mm_set1_ps(fb[27]); \
	bp14=_mm_set1_ps(fb[28]);

#define FILTER_C \
	cm14=_mm_set1_ps(fc[ 0]); \
	cm13=_mm_set1_ps(fc[ 1]); \
	cm12=_mm_set1_ps(fc[ 2]); \
	cm11=_mm_set1_ps(fc[ 3]); \
	cm10=_mm_set1_ps(fc[ 4]); \
	cm09=_mm_set1_ps(fc[ 5]); \
	cm08=_mm_set1_ps(fc[ 6]); \
	cm07=_mm_set1_ps(fc[ 7]); \
	cm06=_mm_set1_ps(fc[ 8]); \
	cm05=_mm_set1_ps(fc[ 9]); \
	cm04=_mm_set1_ps(fc[10]); \
	cm03=_mm_set1_ps(fc[11]); \
	cm02=_mm_set1_ps(fc[12]); \
	cm01=_mm_set1_ps(fc[13]); \
	cp00=_mm_set1_ps(fc[14]); \
	cp01=_mm_set1_ps(fc[15]); \
	cp02=_mm_set1_ps(fc[16]); \
	cp03=_mm_set1_ps(fc[17]); \
	cp04=_mm_set1_ps(fc[18]); \
	cp05=_mm_set1_ps(fc[19]); \
	cp06=_mm_set1_ps(fc[20]); \
	cp07=_mm_set1_ps(fc[21]); \
	cp08=_mm_set1_ps(fc[22]); \
	cp09=_mm_set1_ps(fc[23]); \
	cp10=_mm_set1_ps(fc[24]); \
	cp11=_mm_set1_ps(fc[25]); \
	cp12=_mm_set1_ps(fc[26]); \
	cp13=_mm_set1_ps(fc[27]); \
	cp14=_mm_set1_ps(fc[28]);

#define FILTER_E \
	em14=_mm_set1_ps(fe[ 0]); \
	em13=_mm_set1_ps(fe[ 1]); \
	em12=_mm_set1_ps(fe[ 2]); \
	em11=_mm_set1_ps(fe[ 3]); \
	em10=_mm_set1_ps(fe[ 4]); \
	em09=_mm_set1_ps(fe[ 5]); \
	em08=_mm_set1_ps(fe[ 6]); \
	em07=_mm_set1_ps(fe[ 7]); \
	em06=_mm_set1_ps(fe[ 8]); \
	em05=_mm_set1_ps(fe[ 9]); \
	em04=_mm_set1_ps(fe[10]); \
	em03=_mm_set1_ps(fe[11]); \
	em02=_mm_set1_ps(fe[12]); \
	em01=_mm_set1_ps(fe[13]); \
	ep00=_mm_set1_ps(fe[14]); \
	ep01=_mm_set1_ps(fe[15]); \
	ep02=_mm_set1_ps(fe[16]); \
	ep03=_mm_set1_ps(fe[17]); \
	ep04=_mm_set1_ps(fe[18]); \
	ep05=_mm_set1_ps(fe[19]); \
	ep06=_mm_set1_ps(fe[20]); \
	ep07=_mm_set1_ps(fe[21]); \
	ep08=_mm_set1_ps(fe[22]); \
	ep09=_mm_set1_ps(fe[23]); \
	ep10=_mm_set1_ps(fe[24]); \
	ep11=_mm_set1_ps(fe[25]); \
	ep12=_mm_set1_ps(fe[26]); \
	ep13=_mm_set1_ps(fe[27]); \
	ep14=_mm_set1_ps(fe[28]);

//****************************************************************************************
void convolut_ib_sse_cxyz_(int *n1,int *n2,int *n3,int *ibyz_c,int *ibxz_c,int *ibxy_c,
	float *cprecr,float *x_c,float *y_c,float *fa) {
	int mm,n1p,n2p,n3p,nseg1,nseg2,nseg3,n1t,n2t,n3t,istart;
	__m128  am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01,ap00,ap01;
	__m128  ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14;
	__m128  r0,r1,r2,r3,r4,r5,r6,r7,x,y,cprecrsse;
	__m128i lz,l1,l2,l3,l4,l5,l6,l7,l8;
	int i1,i2,i3,loffset;
	float wai[1024] __attribute__ ((aligned (16)));
	float wao[1024] __attribute__ ((aligned (16)));

	cprecrsse=_mm_set1_ps(cprecr[0]);
	FILTER_A;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*n1+1);
	nseg2=(*n2+1);
	nseg3=nseg1*nseg2;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #pragma omp parallel for schedule(static,1)\
           default(shared) \
           private(i3,i2,i1,istart,n1t,loffset,mm,n1p,wai,wao)\
           private(lz,l1,l2,l3,l4,l5,l6,l7,l8,r0,r1,r2,r3,r4,r5,r6,r7,x,y)
	for(i3=0;i3<*n3+1;i3++) {
	for(i2=0;i2<*n2+1;i2++) {
		istart=ibyz_c[0+2*(i2+nseg2*i3)];
		n1t=ibyz_c[1+2*(i2+nseg2*i3)]-istart+1;
		if(n1t>0) {
		loffset=istart+nseg1*(i2+i3*nseg2);
		//----------------------------------------------------------------------
		for(i1=0;i1<n1t-4;i1+=4) {
			wai[i1+0]=x_c[loffset+i1+0];
			wai[i1+1]=x_c[loffset+i1+1];
			wai[i1+2]=x_c[loffset+i1+2];
			wai[i1+3]=x_c[loffset+i1+3];
		}
		mm=n1t%4;
		if(mm==0) {
			n1p=n1t;
			wai[n1p-4]=x_c[loffset+n1p-4];
			wai[n1p-3]=x_c[loffset+n1p-3];
			wai[n1p-2]=x_c[loffset+n1p-2];
			wai[n1p-1]=x_c[loffset+n1p-1];
		}
		else if(mm==1) {
			n1p=n1t+3;
			wai[n1p-4]=x_c[loffset+n1p-4];
			wai[n1p-3]=0.0;
			wai[n1p-2]=0.0;
			wai[n1p-1]=0.0;
		}
		else if(mm==2) {
			n1p=n1t+2;
			wai[n1p-4]=x_c[loffset+n1p-4];
			wai[n1p-3]=x_c[loffset+n1p-3];
			wai[n1p-2]=0.0;
			wai[n1p-1]=0.0;
		}
		else if(mm==3) {
			n1p=n1t+1;
			wai[n1p-4]=x_c[loffset+n1p-4];
			wai[n1p-3]=x_c[loffset+n1p-3];
			wai[n1p-2]=x_c[loffset+n1p-2];
			wai[n1p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n1t<5) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			CONV_L1(l1,am03,am02,am01,ap00,ap01,ap02,ap03);
			y=_mm_mul_ps(_mm_castsi128_ps(l1),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao[n1p- 4],r0);
		}
		else if(n1t<9) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			CONV_L2(l1,l2,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			y=_mm_mul_ps(_mm_castsi128_ps(l1),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l2),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao[n1p- 8],r1);
			_mm_store_ps(&wao[n1p- 4],r0);
		}
		else if(n1t<13) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			CONV_L3(l1,l2,l3,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			y=_mm_mul_ps(_mm_castsi128_ps(l1),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l2),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l3),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao[n1p-12],r2);
			_mm_store_ps(&wao[n1p- 8],r1);
			_mm_store_ps(&wao[n1p- 4],r0);
		}
		else if(n1t<17) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			CONV_L4(l1,l2,l3,l4,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			y=_mm_mul_ps(_mm_castsi128_ps(l1),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l2),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l3),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l4),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao[n1p-16],r3);
			_mm_store_ps(&wao[n1p-12],r2);
			_mm_store_ps(&wao[n1p- 8],r1);
			_mm_store_ps(&wao[n1p- 4],r0);
		}
		else if(n1t<21) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			CONV_L5(l1,l2,l3,l4,l5,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			y=_mm_mul_ps(_mm_castsi128_ps(l1),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l2),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l3),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l4),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l5),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao[n1p-20],r4);
			_mm_store_ps(&wao[n1p-16],r3);
			_mm_store_ps(&wao[n1p-12],r2);
			_mm_store_ps(&wao[n1p- 8],r1);
			_mm_store_ps(&wao[n1p- 4],r0);
		}
		else if(n1t<25) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			l6=_mm_castps_si128(_mm_load_ps(wai+20));
			CONV_L6(l1,l2,l3,l4,l5,l6,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			y=_mm_mul_ps(_mm_castsi128_ps(l1),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l2),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l3),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l4),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l5),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l6),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao[n1p-24],r5);
			_mm_store_ps(&wao[n1p-20],r4);
			_mm_store_ps(&wao[n1p-16],r3);
			_mm_store_ps(&wao[n1p-12],r2);
			_mm_store_ps(&wao[n1p- 8],r1);
			_mm_store_ps(&wao[n1p- 4],r0);
		}
		else if(n1t<29) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			l6=_mm_castps_si128(_mm_load_ps(wai+20));
			l7=_mm_castps_si128(_mm_load_ps(wai+24));
			CONV_L7(l1,l2,l3,l4,l5,l6,l7,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			y=_mm_mul_ps(_mm_castsi128_ps(l1),cprecrsse);COUNT;r6=_mm_add_ps(r6,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l2),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l3),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l4),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l5),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l6),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l7),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao[n1p-28],r6);
			_mm_store_ps(&wao[n1p-24],r5);
			_mm_store_ps(&wao[n1p-20],r4);
			_mm_store_ps(&wao[n1p-16],r3);
			_mm_store_ps(&wao[n1p-12],r2);
			_mm_store_ps(&wao[n1p- 8],r1);
			_mm_store_ps(&wao[n1p- 4],r0);
		}
		else {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			CONV_L8_PART1(l1,l2,l3,l4,l5,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r1,r2,r3,r4,r5,r6,r7);
			for(i1=28;i1<n1p-4;i1+=4) {
				l6=_mm_castps_si128(_mm_load_ps(wai+i1-8));
				CONV_L8_PART2(l4,l5,l6,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r0,r1,r2,r3,r4,r5,r6,r7);
				y=_mm_mul_ps(_mm_castsi128_ps(l1),cprecrsse);COUNT;r7=_mm_add_ps(r7,y);COUNT;
				_mm_store_ps(&wao[i1-28],r7);
				l1=l2;l2=l3;l3=l4; //this is need because of cprecrsse 
				l4=l5;l5=l6;
				UP_L8(r0,r1,r2,r3,r4,r5,r6,r7);
			}
			l6=_mm_castps_si128(_mm_load_ps(wai+n1p-12));
			l7=_mm_castps_si128(_mm_load_ps(wai+n1p- 8));
			l8=_mm_castps_si128(_mm_load_ps(wai+n1p- 4));
			CONV_L8_PART3(l4,l5,l6,l7,l8,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r0,r1,r2,r3,r4,r5,r6,r7);
			y=_mm_mul_ps(_mm_castsi128_ps(l1),cprecrsse);COUNT;r7=_mm_add_ps(r7,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l2),cprecrsse);COUNT;r6=_mm_add_ps(r6,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l3),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l4),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l5),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l6),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l7),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l8),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao[n1p-32],r7);
			_mm_store_ps(&wao[n1p-28],r6);
			_mm_store_ps(&wao[n1p-24],r5);
			_mm_store_ps(&wao[n1p-20],r4);
			_mm_store_ps(&wao[n1p-16],r3);
			_mm_store_ps(&wao[n1p-12],r2);
			_mm_store_ps(&wao[n1p- 8],r1);
			_mm_store_ps(&wao[n1p- 4],r0);
		}
		//----------------------------------------------------------------------
		for(i1=0;i1<(n1p-4);i1+=4) {
			y_c[loffset+i1+0]=wao[i1+0]; y_c[loffset+i1+1]=wao[i1+1];
			y_c[loffset+i1+2]=wao[i1+2]; y_c[loffset+i1+3]=wao[i1+3];
		}
		for(;i1<n1t;i1++) y_c[loffset+i1]=wao[i1];
		} //end of if for n1t>0
		//----------------------------------------------------------------------
	} //end of loop over i2
	} //end of loop over i3
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #pragma omp parallel for schedule(static,1) \
           default(shared) \
           private(i3,i2,i1,istart,n2t,loffset,mm,n1p,wai,wao)\
           private(lz,l1,l2,l3,l4,l5,l6,l7,l8,r0,r1,r2,r3,r4,r5,r6,r7,x,y)
	for(i3=0;i3<*n3+1;i3++) {
	for(i1=0;i1<*n1+1;i1++) {
		istart=ibxz_c[0+2*(i1+nseg1*i3)];
		n2t=ibxz_c[1+2*(i1+nseg1*i3)]-istart+1;
		if(n2t>0) {
                loffset=i1+nseg1*(istart+i3*nseg2);
		//----------------------------------------------------------------------
		for(i2=0;i2<n2t-4;i2+=4) {
			wai[i2+0]=x_c[loffset+(i2+0)*nseg1];
			wai[i2+1]=x_c[loffset+(i2+1)*nseg1];
			wai[i2+2]=x_c[loffset+(i2+2)*nseg1];
			wai[i2+3]=x_c[loffset+(i2+3)*nseg1];
		}
		mm=n2t%4;
		if(mm==0) {
			n2p=n2t;
			wai[n2p-4]=x_c[loffset+(n2p-4)*nseg1];
			wai[n2p-3]=x_c[loffset+(n2p-3)*nseg1];
			wai[n2p-2]=x_c[loffset+(n2p-2)*nseg1];
			wai[n2p-1]=x_c[loffset+(n2p-1)*nseg1];
		}
		else if(mm==1) {
			n2p=n2t+3;
			wai[n2p-4]=x_c[loffset+(n2p-4)*nseg1];
			wai[n2p-3]=0.0;
			wai[n2p-2]=0.0;
			wai[n2p-1]=0.0;
		}
		else if(mm==2) {
			n2p=n2t+2;
			wai[n2p-4]=x_c[loffset+(n2p-4)*nseg1];
			wai[n2p-3]=x_c[loffset+(n2p-3)*nseg1];
			wai[n2p-2]=0.0;
			wai[n2p-1]=0.0;
		}
		else if(mm==3) {
			n2p=n2t+1;
			wai[n2p-4]=x_c[loffset+(n2p-4)*nseg1];
			wai[n2p-3]=x_c[loffset+(n2p-3)*nseg1];
			wai[n2p-2]=x_c[loffset+(n2p-2)*nseg1];
			wai[n2p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n2t<5) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			CONV_L1(l1,am03,am02,am01,ap00,ap01,ap02,ap03);
			_mm_store_ps(&wao[n2p- 4],r0);
		}
		else if(n2t<9) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			CONV_L2(l1,l2,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			_mm_store_ps(&wao[n2p- 8],r1);
			_mm_store_ps(&wao[n2p- 4],r0);
		}
		else if(n2t<13) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			CONV_L3(l1,l2,l3,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			_mm_store_ps(&wao[n2p-12],r2);
			_mm_store_ps(&wao[n2p- 8],r1);
			_mm_store_ps(&wao[n2p- 4],r0);
		}
		else if(n2t<17) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			CONV_L4(l1,l2,l3,l4,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			_mm_store_ps(&wao[n2p-16],r3);
			_mm_store_ps(&wao[n2p-12],r2);
			_mm_store_ps(&wao[n2p- 8],r1);
			_mm_store_ps(&wao[n2p- 4],r0);
		}
		else if(n2t<21) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			CONV_L5(l1,l2,l3,l4,l5,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			_mm_store_ps(&wao[n2p-20],r4);
			_mm_store_ps(&wao[n2p-16],r3);
			_mm_store_ps(&wao[n2p-12],r2);
			_mm_store_ps(&wao[n2p- 8],r1);
			_mm_store_ps(&wao[n2p- 4],r0);
		}
		else if(n2t<25) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			l6=_mm_castps_si128(_mm_load_ps(wai+20));
			CONV_L6(l1,l2,l3,l4,l5,l6,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			_mm_store_ps(&wao[n2p-24],r5);
			_mm_store_ps(&wao[n2p-20],r4);
			_mm_store_ps(&wao[n2p-16],r3);
			_mm_store_ps(&wao[n2p-12],r2);
			_mm_store_ps(&wao[n2p- 8],r1);
			_mm_store_ps(&wao[n2p- 4],r0);
		}
		else if(n2t<29) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			l6=_mm_castps_si128(_mm_load_ps(wai+20));
			l7=_mm_castps_si128(_mm_load_ps(wai+24));
			CONV_L7(l1,l2,l3,l4,l5,l6,l7,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			_mm_store_ps(&wao[n2p-28],r6);
			_mm_store_ps(&wao[n2p-24],r5);
			_mm_store_ps(&wao[n2p-20],r4);
			_mm_store_ps(&wao[n2p-16],r3);
			_mm_store_ps(&wao[n2p-12],r2);
			_mm_store_ps(&wao[n2p- 8],r1);
			_mm_store_ps(&wao[n2p- 4],r0);
		}
		else {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			CONV_L8_PART1(l1,l2,l3,l4,l5,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r1,r2,r3,r4,r5,r6,r7);
			for(i2=28;i2<n2p-4;i2+=4) {
				l6=_mm_castps_si128(_mm_load_ps(wai+i2-8));
				CONV_L8_PART2(l4,l5,l6,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r0,r1,r2,r3,r4,r5,r6,r7);
				_mm_store_ps(&wao[i2-28],r7);
				l4=l5;l5=l6;
				UP_L8(r0,r1,r2,r3,r4,r5,r6,r7);
			}
			l6=_mm_castps_si128(_mm_load_ps(wai+n2p-12));
			l7=_mm_castps_si128(_mm_load_ps(wai+n2p- 8));
			l8=_mm_castps_si128(_mm_load_ps(wai+n2p- 4));
			CONV_L8_PART3(l4,l5,l6,l7,l8,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r0,r1,r2,r3,r4,r5,r6,r7);
			_mm_store_ps(&wao[n2p-32],r7);
			_mm_store_ps(&wao[n2p-28],r6);
			_mm_store_ps(&wao[n2p-24],r5);
			_mm_store_ps(&wao[n2p-20],r4);
			_mm_store_ps(&wao[n2p-16],r3);
			_mm_store_ps(&wao[n2p-12],r2);
			_mm_store_ps(&wao[n2p- 8],r1);
			_mm_store_ps(&wao[n2p- 4],r0);
		}
		//----------------------------------------------------------------------
		for(i2=0;i2<(n2p-4);i2+=4) {
			y_c[loffset+(i2+0)*nseg1]=y_c[loffset+(i2+0)*nseg1]+wao[i2+0];
			y_c[loffset+(i2+1)*nseg1]=y_c[loffset+(i2+1)*nseg1]+wao[i2+1];
			y_c[loffset+(i2+2)*nseg1]=y_c[loffset+(i2+2)*nseg1]+wao[i2+2];
			y_c[loffset+(i2+3)*nseg1]=y_c[loffset+(i2+3)*nseg1]+wao[i2+3];
		}
		for(;i2<n2t;i2++) y_c[loffset+i2*nseg1]=y_c[loffset+i2*nseg1]+wao[i2];
		} //end of if for n2t>0
		//----------------------------------------------------------------------
	} //end of loop over i1
	} //end of loop over i3
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n3t,loffset,mm,n1p,wai,wao,n3p)\
           private(lz,l1,l2,l3,l4,l5,l6,l7,l8,r0,r1,r2,r3,r4,r5,r6,r7,x,y)
	for(i2=0;i2<*n2+1;i2++) {
	for(i1=0;i1<*n1+1;i1++) {
		istart=ibxy_c[0+2*(i1+nseg1*i2)];
		n3t=ibxy_c[1+2*(i1+nseg1*i2)]-istart+1;
		if(n3t>0) {
		loffset=i1+nseg1*(i2+istart*nseg2);
		//----------------------------------------------------------------------
		for(i3=0;i3<n3t-4;i3+=4) {
			wai[i3+0]=x_c[loffset+(i3+0)*nseg3];
			wai[i3+1]=x_c[loffset+(i3+1)*nseg3];
			wai[i3+2]=x_c[loffset+(i3+2)*nseg3];
			wai[i3+3]=x_c[loffset+(i3+3)*nseg3];
		}
		mm=n3t%4;
		if(mm==0) {
			n3p=n3t;
			wai[n3p-4]=x_c[loffset+(n3p-4)*nseg3];
			wai[n3p-3]=x_c[loffset+(n3p-3)*nseg3];
			wai[n3p-2]=x_c[loffset+(n3p-2)*nseg3];
			wai[n3p-1]=x_c[loffset+(n3p-1)*nseg3];
		}
		else if(mm==1) {
			n3p=n3t+3;
			wai[n3p-4]=x_c[loffset+(n3p-4)*nseg3];
			wai[n3p-3]=0.0;
			wai[n3p-2]=0.0;
			wai[n3p-1]=0.0;
		}
		else if(mm==2) {
			n3p=n3t+2;
			wai[n3p-4]=x_c[loffset+(n3p-4)*nseg3];
			wai[n3p-3]=x_c[loffset+(n3p-3)*nseg3];
			wai[n3p-2]=0.0;
			wai[n3p-1]=0.0;
		}
		else if(mm==3) {
			n3p=n3t+1;
			wai[n3p-4]=x_c[loffset+(n3p-4)*nseg3];
			wai[n3p-3]=x_c[loffset+(n3p-3)*nseg3];
			wai[n3p-2]=x_c[loffset+(n3p-2)*nseg3];
			wai[n3p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n3t<5) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			CONV_L1(l1,am03,am02,am01,ap00,ap01,ap02,ap03);
			_mm_store_ps(&wao[n3p- 4],r0);
		}
		else if(n3t<9) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			CONV_L2(l1,l2,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			_mm_store_ps(&wao[n3p- 8],r1);
			_mm_store_ps(&wao[n3p- 4],r0);
		}
		else if(n3t<13) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			CONV_L3(l1,l2,l3,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			_mm_store_ps(&wao[n3p-12],r2);
			_mm_store_ps(&wao[n3p- 8],r1);
			_mm_store_ps(&wao[n3p- 4],r0);
		}
		else if(n3t<17) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			CONV_L4(l1,l2,l3,l4,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			_mm_store_ps(&wao[n3p-16],r3);
			_mm_store_ps(&wao[n3p-12],r2);
			_mm_store_ps(&wao[n3p- 8],r1);
			_mm_store_ps(&wao[n3p- 4],r0);
		}
		else if(n3t<21) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			CONV_L5(l1,l2,l3,l4,l5,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			_mm_store_ps(&wao[n3p-20],r4);
			_mm_store_ps(&wao[n3p-16],r3);
			_mm_store_ps(&wao[n3p-12],r2);
			_mm_store_ps(&wao[n3p- 8],r1);
			_mm_store_ps(&wao[n3p- 4],r0);
		}
		else if(n3t<25) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			l6=_mm_castps_si128(_mm_load_ps(wai+20));
			CONV_L6(l1,l2,l3,l4,l5,l6,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			_mm_store_ps(&wao[n3p-24],r5);
			_mm_store_ps(&wao[n3p-20],r4);
			_mm_store_ps(&wao[n3p-16],r3);
			_mm_store_ps(&wao[n3p-12],r2);
			_mm_store_ps(&wao[n3p- 8],r1);
			_mm_store_ps(&wao[n3p- 4],r0);
		}
		else if(n3t<29) {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			l6=_mm_castps_si128(_mm_load_ps(wai+20));
			l7=_mm_castps_si128(_mm_load_ps(wai+24));
			CONV_L7(l1,l2,l3,l4,l5,l6,l7,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			_mm_store_ps(&wao[n3p-28],r6);
			_mm_store_ps(&wao[n3p-24],r5);
			_mm_store_ps(&wao[n3p-20],r4);
			_mm_store_ps(&wao[n3p-16],r3);
			_mm_store_ps(&wao[n3p-12],r2);
			_mm_store_ps(&wao[n3p- 8],r1);
			_mm_store_ps(&wao[n3p- 4],r0);
		}
		else {
			l1=_mm_castps_si128(_mm_load_ps(wai));
			l2=_mm_castps_si128(_mm_load_ps(wai+4));
			l3=_mm_castps_si128(_mm_load_ps(wai+8));
			l4=_mm_castps_si128(_mm_load_ps(wai+12));
			l5=_mm_castps_si128(_mm_load_ps(wai+16));
			CONV_L8_PART1(l1,l2,l3,l4,l5,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r1,r2,r3,r4,r5,r6,r7);
			for(i3=28;i3<n3p-4;i3+=4) {
				l6=_mm_castps_si128(_mm_load_ps(wai+i3-8));
				CONV_L8_PART2(l4,l5,l6,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r0,r1,r2,r3,r4,r5,r6,r7);
				_mm_store_ps(&wao[i3-28],r7);
				l4=l5;l5=l6;
				UP_L8(r0,r1,r2,r3,r4,r5,r6,r7);
			}
			l6=_mm_castps_si128(_mm_load_ps(wai+n3p-12));
			l7=_mm_castps_si128(_mm_load_ps(wai+n3p- 8));
			l8=_mm_castps_si128(_mm_load_ps(wai+n3p- 4));
			CONV_L8_PART3(l4,l5,l6,l7,l8,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r0,r1,r2,r3,r4,r5,r6,r7);
			_mm_store_ps(&wao[n3p-32],r7);
			_mm_store_ps(&wao[n3p-28],r6);
			_mm_store_ps(&wao[n3p-24],r5);
			_mm_store_ps(&wao[n3p-20],r4);
			_mm_store_ps(&wao[n3p-16],r3);
			_mm_store_ps(&wao[n3p-12],r2);
			_mm_store_ps(&wao[n3p- 8],r1);
			_mm_store_ps(&wao[n3p- 4],r0);
		}
		//----------------------------------------------------------------------
		for(i3=0;i3<(n3p-4);i3+=4) {
			y_c[loffset+(i3+0)*nseg3]=y_c[loffset+(i3+0)*nseg3]+wao[i3+0];
			y_c[loffset+(i3+1)*nseg3]=y_c[loffset+(i3+1)*nseg3]+wao[i3+1];
			y_c[loffset+(i3+2)*nseg3]=y_c[loffset+(i3+2)*nseg3]+wao[i3+2];
			y_c[loffset+(i3+3)*nseg3]=y_c[loffset+(i3+3)*nseg3]+wao[i3+3];
		}
		for(;i3<n3t;i3++) y_c[loffset+i3*nseg3]=y_c[loffset+i3*nseg3]+wao[i3];
		} //end of if for n3t>0
	} //end of loop over i1
	} //end of loop over i2
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
} // end of function convolut_ib_sse_cxyz_
//****************************************************************************************
void convolut_ib_sse_fx_wss_(int *n1,int *n2,int *n3,int *ibyz_f,float *cprecr,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n1p,nseg1,nseg2,n1t,istart;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y,cprecrsse;
	__m128i lz;
	__m128i l11,l12,l13,l14,l15,l16,l17,l18;
	__m128 r10,r11,r12,r13,r14,r15,r16,r17;
	int i1,i2,i3,loffset1;
	float wai1[1024] __attribute__ ((aligned (16)));
	float wao1[1024] __attribute__ ((aligned (16)));
	cprecrsse=_mm_set1_ps(cprecr[0]);
	FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along X )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n1t,loffset1,mm,n1p,wai1,wao1)\
           private(l11,l12,l13,l14,l15,l16,l17,l18,r10,r11,r12,r13,r14,r15,r16,r17)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i3=*nfl3;i3<=*nfu3;i3++) {
	for(i2=*nfl2;i2<=*nfu2;i2++) {
		istart=ibyz_f[0+2*(i2+(*n2+1)*i3)];
		n1t=ibyz_f[1+2*(i2+(*n2+1)*i3)]-istart+1;
		if(n1t>0) {
		loffset1=7*(istart-*nfl1+nseg1*(i2-*nfl2+(i3-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i1=0;i1<n1t-4;i1+=4) {
			wai1[i1+0]=x_f[loffset1+(i1+0)*7];
			wai1[i1+1]=x_f[loffset1+(i1+1)*7];
			wai1[i1+2]=x_f[loffset1+(i1+2)*7];
			wai1[i1+3]=x_f[loffset1+(i1+3)*7];
		}
		mm=n1t%4;
		if(mm==0) {
			n1p=n1t;
			wai1[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai1[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai1[n1p-2]=x_f[loffset1+(n1p-2)*7];
			wai1[n1p-1]=x_f[loffset1+(n1p-1)*7];
		}
		else if(mm==1) {
			n1p=n1t+3;
			wai1[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai1[n1p-3]=0.0;
			wai1[n1p-2]=0.0;
			wai1[n1p-1]=0.0;
		}
		else if(mm==2) {
			n1p=n1t+2;
			wai1[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai1[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai1[n1p-2]=0.0;
			wai1[n1p-1]=0.0;
		}
		else if(mm==3) {
			n1p=n1t+1;
			wai1[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai1[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai1[n1p-2]=x_f[loffset1+(n1p-2)*7];
			wai1[n1p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n1t<5) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			CONV_L1(l11,em03,em02,em01,ep00,ep01,ep02,ep03);
			y=_mm_mul_ps(_mm_castsi128_ps(l11),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao1[n1p- 4],r0);
		}
		else if(n1t<9) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			CONV_L2(l11,l12,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			y=_mm_mul_ps(_mm_castsi128_ps(l11),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l12),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao1[n1p- 8],r1);
			_mm_store_ps(&wao1[n1p- 4],r0);
		}
		else if(n1t<13) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			CONV_L3(l11,l12,l13,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			y=_mm_mul_ps(_mm_castsi128_ps(l11),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l12),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l13),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao1[n1p-12],r2);
			_mm_store_ps(&wao1[n1p- 8],r1);
			_mm_store_ps(&wao1[n1p- 4],r0);
		}
		else if(n1t<17) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));

			CONV_L4(l11,l12,l13,l14,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l11),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l12),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l13),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l14),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao1[n1p-16],r3);
			_mm_store_ps(&wao1[n1p-12],r2);
			_mm_store_ps(&wao1[n1p- 8],r1);
			_mm_store_ps(&wao1[n1p- 4],r0);
		}
		else if(n1t<21) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));

			CONV_L5(l11,l12,l13,l14,l15,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l11),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l12),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l13),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l14),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l15),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao1[n1p-20],r4);
			_mm_store_ps(&wao1[n1p-16],r3);
			_mm_store_ps(&wao1[n1p-12],r2);
			_mm_store_ps(&wao1[n1p- 8],r1);
			_mm_store_ps(&wao1[n1p- 4],r0);
		}
		else if(n1t<25) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));
			l16=_mm_castps_si128(_mm_load_ps(wai1+20));

			CONV_L6(l11,l12,l13,l14,l15,l16,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l11),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l12),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l13),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l14),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l15),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l16),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao1[n1p-24],r5);
			_mm_store_ps(&wao1[n1p-20],r4);
			_mm_store_ps(&wao1[n1p-16],r3);
			_mm_store_ps(&wao1[n1p-12],r2);
			_mm_store_ps(&wao1[n1p- 8],r1);
			_mm_store_ps(&wao1[n1p- 4],r0);
		}
		else if(n1t<29) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));
			l16=_mm_castps_si128(_mm_load_ps(wai1+20));
			l17=_mm_castps_si128(_mm_load_ps(wai1+24));

			CONV_L7(l11,l12,l13,l14,l15,l16,l17,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l11),cprecrsse);COUNT;r6=_mm_add_ps(r6,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l12),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l13),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l14),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l15),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l16),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l17),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao1[n1p-28],r6);
			_mm_store_ps(&wao1[n1p-24],r5);
			_mm_store_ps(&wao1[n1p-20],r4);
			_mm_store_ps(&wao1[n1p-16],r3);
			_mm_store_ps(&wao1[n1p-12],r2);
			_mm_store_ps(&wao1[n1p- 8],r1);
			_mm_store_ps(&wao1[n1p- 4],r0);
		}
		else {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));

			CONV_L8_PART1(l11,l12,l13,l14,l15,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r11,r12,r13,r14,r15,r16,r17);

			for(i1=28;i1<n1p-4;i1+=4) {
				l16=_mm_castps_si128(_mm_load_ps(wai1+i1-8));

				CONV_L8_PART2(l14,l15,l16,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r10,r11,r12,r13,r14,r15,r16,r17);
				y=_mm_mul_ps(_mm_castsi128_ps(l11),cprecrsse);COUNT;r17=_mm_add_ps(r17,y);COUNT;
				_mm_store_ps(&wao1[i1-28],r17);

				l11=l12;l12=l13;l13=l14; //this is need because of cprecrsse 
				l14=l15;l15=l16;
				UP_L8(r10,r11,r12,r13,r14,r15,r16,r17);

			}
			l16=_mm_castps_si128(_mm_load_ps(wai1+n1p-12));
			l17=_mm_castps_si128(_mm_load_ps(wai1+n1p- 8));
			l18=_mm_castps_si128(_mm_load_ps(wai1+n1p- 4));

			CONV_L8_PART3(l14,l15,l16,l17,l18,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r10,r11,r12,r13,r14,r15,r16,r17);
			y=_mm_mul_ps(_mm_castsi128_ps(l11),cprecrsse);COUNT;r17=_mm_add_ps(r17,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l12),cprecrsse);COUNT;r16=_mm_add_ps(r16,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l13),cprecrsse);COUNT;r15=_mm_add_ps(r15,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l14),cprecrsse);COUNT;r14=_mm_add_ps(r14,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l15),cprecrsse);COUNT;r13=_mm_add_ps(r13,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l16),cprecrsse);COUNT;r12=_mm_add_ps(r12,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l17),cprecrsse);COUNT;r11=_mm_add_ps(r11,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l18),cprecrsse);COUNT;r10=_mm_add_ps(r10,y);COUNT;
			_mm_store_ps(&wao1[n1p-32],r17);
			_mm_store_ps(&wao1[n1p-28],r16);
			_mm_store_ps(&wao1[n1p-24],r15);
			_mm_store_ps(&wao1[n1p-20],r14);
			_mm_store_ps(&wao1[n1p-16],r13);
			_mm_store_ps(&wao1[n1p-12],r12);
			_mm_store_ps(&wao1[n1p- 8],r11);
			_mm_store_ps(&wao1[n1p- 4],r10);

		}
		//----------------------------------------------------------------------
		for(i1=0;i1<(n1p-4);i1+=4) {
			y_f[loffset1+(i1+0)*7]=wao1[i1+0];
			y_f[loffset1+(i1+1)*7]=wao1[i1+1];
			y_f[loffset1+(i1+2)*7]=wao1[i1+2];
			y_f[loffset1+(i1+3)*7]=wao1[i1+3];
		}
		for(;i1<n1t;i1++) {
			y_f[loffset1+i1*7]=wao1[i1];
		}
		} //end of if for n1t>0
		//----------------------------------------------------------------------
	} //end of loop over i2
	} //end of loop over i3
} // end of function convolut_ib_sse_fx_wss_
//****************************************************************************************
void convolut_ib_sse_fx_sws_wws_(int *n1,int *n2,int *n3,int *ibyz_f,float *cprecr,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n1p,nseg1,nseg2,n1t,istart;
	__m128 am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01,ap00,ap01;
	__m128 ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14;
	__m128 bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01;
	__m128 bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14;
	__m128 cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01;
	__m128 cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y,cprecrsse;
	__m128i lz;
	__m128i l21,l22,l23,l24,l25,l26,l27,l28;
	__m128i l31,l32,l33,l34,l35,l36,l37,l38;
	__m128 r20,r21,r22,r23,r24,r25,r26,r27;
	__m128 r30,r31,r32,r33,r34,r35,r36,r37;
	int i1,i2,i3,loffset1,loffset2;
	float wai2[1024] __attribute__ ((aligned (16)));
	float wai3[1024] __attribute__ ((aligned (16)));
	float wao2[1024] __attribute__ ((aligned (16)));
	float wao3[1024] __attribute__ ((aligned (16)));
	cprecrsse=_mm_set1_ps(cprecr[0]);
	FILTER_A; FILTER_B; FILTER_C; FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along X )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n1t,loffset1,loffset2,mm,n1p,wai2,wai3,wao2,wao3)\
           private(l21,l22,l23,l24,l25,l26,l27,l28,l31,l32,l33,l34,l35,l36,l37,l38)\
           private(r20,r21,r22,r23,r24,r25,r26,r27,r30,r31,r32,r33,r34,r35,r36,r37)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i3=*nfl3;i3<=*nfu3;i3++) {
	for(i2=*nfl2;i2<=*nfu2;i2++) {
		istart=ibyz_f[0+2*(i2+(*n2+1)*i3)];
		n1t=ibyz_f[1+2*(i2+(*n2+1)*i3)]-istart+1;
		if(n1t>0) {
		loffset1=1+7*(istart-*nfl1+nseg1*(i2-*nfl2+(i3-*nfl3)*nseg2));
		loffset2=2+7*(istart-*nfl1+nseg1*(i2-*nfl2+(i3-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i1=0;i1<n1t-4;i1+=4) {
			wai2[i1+0]=x_f[loffset1+(i1+0)*7];
			wai3[i1+0]=x_f[loffset2+(i1+0)*7];
			wai2[i1+1]=x_f[loffset1+(i1+1)*7];
			wai3[i1+1]=x_f[loffset2+(i1+1)*7];
			wai2[i1+2]=x_f[loffset1+(i1+2)*7];
			wai3[i1+2]=x_f[loffset2+(i1+2)*7];
			wai2[i1+3]=x_f[loffset1+(i1+3)*7];
			wai3[i1+3]=x_f[loffset2+(i1+3)*7];
		}
		mm=n1t%4;
		if(mm==0) {
			n1p=n1t;
			wai2[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai3[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai2[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai3[n1p-3]=x_f[loffset2+(n1p-3)*7];
			wai2[n1p-2]=x_f[loffset1+(n1p-2)*7];
			wai3[n1p-2]=x_f[loffset2+(n1p-2)*7];
			wai2[n1p-1]=x_f[loffset1+(n1p-1)*7];
			wai3[n1p-1]=x_f[loffset2+(n1p-1)*7];
		}
		else if(mm==1) {
			n1p=n1t+3;
			wai2[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai3[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai2[n1p-3]=0.0;
			wai3[n1p-3]=0.0;
			wai2[n1p-2]=0.0;
			wai3[n1p-2]=0.0;
			wai2[n1p-1]=0.0;
			wai3[n1p-1]=0.0;
		}
		else if(mm==2) {
			n1p=n1t+2;
			wai2[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai3[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai2[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai3[n1p-3]=x_f[loffset2+(n1p-3)*7];
			wai2[n1p-2]=0.0;
			wai3[n1p-2]=0.0;
			wai2[n1p-1]=0.0;
			wai3[n1p-1]=0.0;
		}
		else if(mm==3) {
			n1p=n1t+1;
			wai2[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai3[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai2[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai3[n1p-3]=x_f[loffset2+(n1p-3)*7];
			wai2[n1p-2]=x_f[loffset1+(n1p-2)*7];
			wai3[n1p-2]=x_f[loffset2+(n1p-2)*7];
			wai2[n1p-1]=0.0;
			wai3[n1p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n1t<5) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			CONV_L1(l21,am03,am02,am01,ap00,ap01,ap02,ap03);
			CONV_L1C(l31,bm03,bm02,bm01,bp00,bp01,bp02,bp03);
			y=_mm_mul_ps(_mm_castsi128_ps(l21),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao2[n1p- 4],r0);
			CONV_L1(l21,cm03,cm02,cm01,cp00,cp01,cp02,cp03);
			CONV_L1C(l31,em03,em02,em01,ep00,ep01,ep02,ep03);
			y=_mm_mul_ps(_mm_castsi128_ps(l31),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao3[n1p- 4],r0);
		}
		else if(n1t<9) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			CONV_L2(l21,l22,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			CONV_L2C(l31,l32,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07);
			y=_mm_mul_ps(_mm_castsi128_ps(l21),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l22),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao2[n1p- 8],r1);
			_mm_store_ps(&wao2[n1p- 4],r0);
			CONV_L2(l21,l22,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07);
			CONV_L2C(l31,l32,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			y=_mm_mul_ps(_mm_castsi128_ps(l31),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l32),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao3[n1p- 8],r1);
			_mm_store_ps(&wao3[n1p- 4],r0);
		}
		else if(n1t<13) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			CONV_L3(l21,l22,l23,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			CONV_L3C(l31,l32,l33,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11);
			y=_mm_mul_ps(_mm_castsi128_ps(l21),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l22),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l23),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao2[n1p-12],r2);
			_mm_store_ps(&wao2[n1p- 8],r1);
			_mm_store_ps(&wao2[n1p- 4],r0);
			CONV_L3(l21,l22,l23,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11);
			CONV_L3C(l31,l32,l33,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			y=_mm_mul_ps(_mm_castsi128_ps(l31),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l32),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l33),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao3[n1p-12],r2);
			_mm_store_ps(&wao3[n1p- 8],r1);
			_mm_store_ps(&wao3[n1p- 4],r0);
		}
		else if(n1t<17) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));

			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));

			CONV_L4(l21,l22,l23,l24,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L4C(l31,l32,l33,l34,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l21),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l22),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l23),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l24),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao2[n1p-16],r3);
			_mm_store_ps(&wao2[n1p-12],r2);
			_mm_store_ps(&wao2[n1p- 8],r1);
			_mm_store_ps(&wao2[n1p- 4],r0);
			CONV_L4(l21,l22,l23,l24,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L4C(l31,l32,l33,l34,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l31),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l32),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l33),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l34),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao3[n1p-16],r3);
			_mm_store_ps(&wao3[n1p-12],r2);
			_mm_store_ps(&wao3[n1p- 8],r1);
			_mm_store_ps(&wao3[n1p- 4],r0);
		}
		else if(n1t<21) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));

			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));

			CONV_L5(l21,l22,l23,l24,l25,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L5C(l31,l32,l33,l34,l35,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l21),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l22),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l23),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l24),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l25),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao2[n1p-20],r4);
			_mm_store_ps(&wao2[n1p-16],r3);
			_mm_store_ps(&wao2[n1p-12],r2);
			_mm_store_ps(&wao2[n1p- 8],r1);
			_mm_store_ps(&wao2[n1p- 4],r0);
			CONV_L5(l21,l22,l23,l24,l25,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L5C(l31,l32,l33,l34,l35,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l31),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l32),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l33),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l34),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l35),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao3[n1p-20],r4);
			_mm_store_ps(&wao3[n1p-16],r3);
			_mm_store_ps(&wao3[n1p-12],r2);
			_mm_store_ps(&wao3[n1p- 8],r1);
			_mm_store_ps(&wao3[n1p- 4],r0);
		}
		else if(n1t<25) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));
			l26=_mm_castps_si128(_mm_load_ps(wai2+20));

			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));
			l36=_mm_castps_si128(_mm_load_ps(wai3+20));

			CONV_L6(l21,l22,l23,l24,l25,l26,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L6C(l31,l32,l33,l34,l35,l36,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l21),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l22),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l23),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l24),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l25),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l26),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao2[n1p-24],r5);
			_mm_store_ps(&wao2[n1p-20],r4);
			_mm_store_ps(&wao2[n1p-16],r3);
			_mm_store_ps(&wao2[n1p-12],r2);
			_mm_store_ps(&wao2[n1p- 8],r1);
			_mm_store_ps(&wao2[n1p- 4],r0);
			CONV_L6(l21,l22,l23,l24,l25,l26,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L6C(l31,l32,l33,l34,l35,l36,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l31),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l32),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l33),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l34),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l35),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l36),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao3[n1p-24],r5);
			_mm_store_ps(&wao3[n1p-20],r4);
			_mm_store_ps(&wao3[n1p-16],r3);
			_mm_store_ps(&wao3[n1p-12],r2);
			_mm_store_ps(&wao3[n1p- 8],r1);
			_mm_store_ps(&wao3[n1p- 4],r0);
		}
		else if(n1t<29) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));
			l26=_mm_castps_si128(_mm_load_ps(wai2+20));
			l27=_mm_castps_si128(_mm_load_ps(wai2+24));

			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));
			l36=_mm_castps_si128(_mm_load_ps(wai3+20));
			l37=_mm_castps_si128(_mm_load_ps(wai3+24));

			CONV_L7(l21,l22,l23,l24,l25,l26,l27,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L7C(l31,l32,l33,l34,l35,l36,l37,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l21),cprecrsse);COUNT;r6=_mm_add_ps(r6,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l22),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l23),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l24),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l25),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l26),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l27),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao2[n1p-28],r6);
			_mm_store_ps(&wao2[n1p-24],r5);
			_mm_store_ps(&wao2[n1p-20],r4);
			_mm_store_ps(&wao2[n1p-16],r3);
			_mm_store_ps(&wao2[n1p-12],r2);
			_mm_store_ps(&wao2[n1p- 8],r1);
			_mm_store_ps(&wao2[n1p- 4],r0);
			CONV_L7(l21,l22,l23,l24,l25,l26,l27,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L7C(l31,l32,l33,l34,l35,l36,l37,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l31),cprecrsse);COUNT;r6=_mm_add_ps(r6,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l32),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l33),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l34),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l35),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l36),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l37),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao3[n1p-28],r6);
			_mm_store_ps(&wao3[n1p-24],r5);
			_mm_store_ps(&wao3[n1p-20],r4);
			_mm_store_ps(&wao3[n1p-16],r3);
			_mm_store_ps(&wao3[n1p-12],r2);
			_mm_store_ps(&wao3[n1p- 8],r1);
			_mm_store_ps(&wao3[n1p- 4],r0);
		}
		else {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));

			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));

			CONV_L8_PART1(l21,l22,l23,l24,l25,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r21,r22,r23,r24,r25,r26,r27);
			CONV_L8_PART1C(l31,l32,l33,l34,l35,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,r21,r22,r23,r24,r25,r26,r27);

			CONV_L8_PART1(l21,l22,l23,l24,l25,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,r31,r32,r33,r34,r35,r36,r37);
			CONV_L8_PART1C(l31,l32,l33,l34,l35,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r31,r32,r33,r34,r35,r36,r37);

			for(i1=28;i1<n1p-4;i1+=4) {
				l26=_mm_castps_si128(_mm_load_ps(wai2+i1-8));
				l36=_mm_castps_si128(_mm_load_ps(wai3+i1-8));

				CONV_L8_PART2(l24,l25,l26,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r20,r21,r22,r23,r24,r25,r26,r27);
				CONV_L8_PART2C(l34,l35,l36,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
					bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r20,r21,r22,r23,r24,r25,r26,r27);
				y=_mm_mul_ps(_mm_castsi128_ps(l21),cprecrsse);COUNT;r27=_mm_add_ps(r27,y);COUNT;
				_mm_store_ps(&wao2[i1-28],r27);
				CONV_L8_PART2(l24,l25,l26,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
					cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r30,r31,r32,r33,r34,r35,r36,r37);
				CONV_L8_PART2C(l34,l35,l36,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r30,r31,r32,r33,r34,r35,r36,r37);
				y=_mm_mul_ps(_mm_castsi128_ps(l31),cprecrsse);COUNT;r37=_mm_add_ps(r37,y);COUNT;
				_mm_store_ps(&wao3[i1-28],r37);

				l21=l22;l22=l23;l23=l24; //this is need because of cprecrsse 
				l24=l25;l25=l26;
				UP_L8(r20,r21,r22,r23,r24,r25,r26,r27);

				l31=l32;l32=l33;l33=l34; //this is need because of cprecrsse 
				l34=l35;l35=l36;
				UP_L8(r30,r31,r32,r33,r34,r35,r36,r37);

			}
			l26=_mm_castps_si128(_mm_load_ps(wai2+n1p-12));
			l27=_mm_castps_si128(_mm_load_ps(wai2+n1p- 8));
			l28=_mm_castps_si128(_mm_load_ps(wai2+n1p- 4));

			l36=_mm_castps_si128(_mm_load_ps(wai3+n1p-12));
			l37=_mm_castps_si128(_mm_load_ps(wai3+n1p- 8));
			l38=_mm_castps_si128(_mm_load_ps(wai3+n1p- 4));

			CONV_L8_PART3(l24,l25,l26,l27,l28,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r20,r21,r22,r23,r24,r25,r26,r27);
			CONV_L8_PART3C(l34,l35,l36,l37,l38,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r20,r21,r22,r23,r24,r25,r26,r27);
			y=_mm_mul_ps(_mm_castsi128_ps(l21),cprecrsse);COUNT;r27=_mm_add_ps(r27,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l22),cprecrsse);COUNT;r26=_mm_add_ps(r26,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l23),cprecrsse);COUNT;r25=_mm_add_ps(r25,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l24),cprecrsse);COUNT;r24=_mm_add_ps(r24,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l25),cprecrsse);COUNT;r23=_mm_add_ps(r23,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l26),cprecrsse);COUNT;r22=_mm_add_ps(r22,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l27),cprecrsse);COUNT;r21=_mm_add_ps(r21,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l28),cprecrsse);COUNT;r20=_mm_add_ps(r20,y);COUNT;
			_mm_store_ps(&wao2[n1p-32],r27);
			_mm_store_ps(&wao2[n1p-28],r26);
			_mm_store_ps(&wao2[n1p-24],r25);
			_mm_store_ps(&wao2[n1p-20],r24);
			_mm_store_ps(&wao2[n1p-16],r23);
			_mm_store_ps(&wao2[n1p-12],r22);
			_mm_store_ps(&wao2[n1p- 8],r21);
			_mm_store_ps(&wao2[n1p- 4],r20);
			CONV_L8_PART3(l24,l25,l26,l27,l28,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r30,r31,r32,r33,r34,r35,r36,r37);
			CONV_L8_PART3C(l34,l35,l36,l37,l38,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r30,r31,r32,r33,r34,r35,r36,r37);
			y=_mm_mul_ps(_mm_castsi128_ps(l31),cprecrsse);COUNT;r37=_mm_add_ps(r37,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l32),cprecrsse);COUNT;r36=_mm_add_ps(r36,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l33),cprecrsse);COUNT;r35=_mm_add_ps(r35,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l34),cprecrsse);COUNT;r34=_mm_add_ps(r34,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l35),cprecrsse);COUNT;r33=_mm_add_ps(r33,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l36),cprecrsse);COUNT;r32=_mm_add_ps(r32,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l37),cprecrsse);COUNT;r31=_mm_add_ps(r31,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l38),cprecrsse);COUNT;r30=_mm_add_ps(r30,y);COUNT;
			_mm_store_ps(&wao3[n1p-32],r37);
			_mm_store_ps(&wao3[n1p-28],r36);
			_mm_store_ps(&wao3[n1p-24],r35);
			_mm_store_ps(&wao3[n1p-20],r34);
			_mm_store_ps(&wao3[n1p-16],r33);
			_mm_store_ps(&wao3[n1p-12],r32);
			_mm_store_ps(&wao3[n1p- 8],r31);
			_mm_store_ps(&wao3[n1p- 4],r30);
		}
		//----------------------------------------------------------------------
		for(i1=0;i1<(n1p-4);i1+=4) {
			y_f[loffset1+(i1+0)*7]=wao2[i1+0];
			y_f[loffset2+(i1+0)*7]=wao3[i1+0];
			y_f[loffset1+(i1+1)*7]=wao2[i1+1];
			y_f[loffset2+(i1+1)*7]=wao3[i1+1];
			y_f[loffset1+(i1+2)*7]=wao2[i1+2];
			y_f[loffset2+(i1+2)*7]=wao3[i1+2];
			y_f[loffset1+(i1+3)*7]=wao2[i1+3];
			y_f[loffset2+(i1+3)*7]=wao3[i1+3];
		}
		for(;i1<n1t;i1++) {
			y_f[loffset1+i1*7]=wao2[i1];
			y_f[loffset2+i1*7]=wao3[i1];
		}
		} //end of if for n1t>0
		//----------------------------------------------------------------------
	} //end of loop over i2
	} //end of loop over i3
} // end of function convolut_ib_sse_fx_sws_wws_
//****************************************************************************************
void convolut_ib_sse_fx_ssw_wsw_(int *n1,int *n2,int *n3,int *ibyz_f,float *cprecr,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n1p,nseg1,nseg2,n1t,istart;
	__m128 am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01,ap00,ap01;
	__m128 ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14;
	__m128 bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01;
	__m128 bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14;
	__m128 cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01;
	__m128 cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y,cprecrsse;
	__m128i lz;
	__m128i l41,l42,l43,l44,l45,l46,l47,l48;
	__m128i l51,l52,l53,l54,l55,l56,l57,l58;
	__m128 r40,r41,r42,r43,r44,r45,r46,r47;
	__m128 r50,r51,r52,r53,r54,r55,r56,r57;
	int i1,i2,i3,loffset1,loffset2;
	float wai4[1024] __attribute__ ((aligned (16)));
	float wai5[1024] __attribute__ ((aligned (16)));
	float wao4[1024] __attribute__ ((aligned (16)));
	float wao5[1024] __attribute__ ((aligned (16)));
	cprecrsse=_mm_set1_ps(cprecr[0]);
	FILTER_A; FILTER_B; FILTER_C; FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along X )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n1t,loffset1,loffset2,mm,n1p,wai4,wai5,wao4,wao5)\
           private(l41,l42,l43,l44,l45,l46,l47,l48,l51,l52,l53,l54,l55,l56,l57,l58)\
           private(r40,r41,r42,r43,r44,r45,r46,r47,r50,r51,r52,r53,r54,r55,r56,r57)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i3=*nfl3;i3<=*nfu3;i3++) {
	for(i2=*nfl2;i2<=*nfu2;i2++) {
		istart=ibyz_f[0+2*(i2+(*n2+1)*i3)];
		n1t=ibyz_f[1+2*(i2+(*n2+1)*i3)]-istart+1;
		if(n1t>0) {
		loffset1=3+7*(istart-*nfl1+nseg1*(i2-*nfl2+(i3-*nfl3)*nseg2));
		loffset2=4+7*(istart-*nfl1+nseg1*(i2-*nfl2+(i3-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i1=0;i1<n1t-4;i1+=4) {
			wai4[i1+0]=x_f[loffset1+(i1+0)*7];
			wai5[i1+0]=x_f[loffset2+(i1+0)*7];
			wai4[i1+1]=x_f[loffset1+(i1+1)*7];
			wai5[i1+1]=x_f[loffset2+(i1+1)*7];
			wai4[i1+2]=x_f[loffset1+(i1+2)*7];
			wai5[i1+2]=x_f[loffset2+(i1+2)*7];
			wai4[i1+3]=x_f[loffset1+(i1+3)*7];
			wai5[i1+3]=x_f[loffset2+(i1+3)*7];
		}
		mm=n1t%4;
		if(mm==0) {
			n1p=n1t;
			wai4[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai5[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai4[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai5[n1p-3]=x_f[loffset2+(n1p-3)*7];
			wai4[n1p-2]=x_f[loffset1+(n1p-2)*7];
			wai5[n1p-2]=x_f[loffset2+(n1p-2)*7];
			wai4[n1p-1]=x_f[loffset1+(n1p-1)*7];
			wai5[n1p-1]=x_f[loffset2+(n1p-1)*7];
		}
		else if(mm==1) {
			n1p=n1t+3;
			wai4[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai5[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai4[n1p-3]=0.0;
			wai5[n1p-3]=0.0;
			wai4[n1p-2]=0.0;
			wai5[n1p-2]=0.0;
			wai4[n1p-1]=0.0;
			wai5[n1p-1]=0.0;
		}
		else if(mm==2) {
			n1p=n1t+2;
			wai4[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai5[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai4[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai5[n1p-3]=x_f[loffset2+(n1p-3)*7];
			wai4[n1p-2]=0.0;
			wai5[n1p-2]=0.0;
			wai4[n1p-1]=0.0;
			wai5[n1p-1]=0.0;
		}
		else if(mm==3) {
			n1p=n1t+1;
			wai4[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai5[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai4[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai5[n1p-3]=x_f[loffset2+(n1p-3)*7];
			wai4[n1p-2]=x_f[loffset1+(n1p-2)*7];
			wai5[n1p-2]=x_f[loffset2+(n1p-2)*7];
			wai4[n1p-1]=0.0;
			wai5[n1p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n1t<5) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			CONV_L1(l41,am03,am02,am01,ap00,ap01,ap02,ap03);
			CONV_L1C(l51,bm03,bm02,bm01,bp00,bp01,bp02,bp03);
			y=_mm_mul_ps(_mm_castsi128_ps(l41),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao4[n1p- 4],r0);
			CONV_L1(l41,cm03,cm02,cm01,cp00,cp01,cp02,cp03);
			CONV_L1C(l51,em03,em02,em01,ep00,ep01,ep02,ep03);
			y=_mm_mul_ps(_mm_castsi128_ps(l51),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao5[n1p- 4],r0);
		}
		else if(n1t<9) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			CONV_L2(l41,l42,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			CONV_L2C(l51,l52,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07);
			y=_mm_mul_ps(_mm_castsi128_ps(l41),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l42),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao4[n1p- 8],r1);
			_mm_store_ps(&wao4[n1p- 4],r0);
			CONV_L2(l41,l42,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07);
			CONV_L2C(l51,l52,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			y=_mm_mul_ps(_mm_castsi128_ps(l51),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l52),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao5[n1p- 8],r1);
			_mm_store_ps(&wao5[n1p- 4],r0);
		}
		else if(n1t<13) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			CONV_L3(l41,l42,l43,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			CONV_L3C(l51,l52,l53,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11);
			y=_mm_mul_ps(_mm_castsi128_ps(l41),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l42),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l43),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao4[n1p-12],r2);
			_mm_store_ps(&wao4[n1p- 8],r1);
			_mm_store_ps(&wao4[n1p- 4],r0);
			CONV_L3(l41,l42,l43,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11);
			CONV_L3C(l51,l52,l53,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			y=_mm_mul_ps(_mm_castsi128_ps(l51),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l52),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l53),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao5[n1p-12],r2);
			_mm_store_ps(&wao5[n1p- 8],r1);
			_mm_store_ps(&wao5[n1p- 4],r0);
		}
		else if(n1t<17) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));

			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));

			CONV_L4(l41,l42,l43,l44,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L4C(l51,l52,l53,l54,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l41),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l42),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l43),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l44),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao4[n1p-16],r3);
			_mm_store_ps(&wao4[n1p-12],r2);
			_mm_store_ps(&wao4[n1p- 8],r1);
			_mm_store_ps(&wao4[n1p- 4],r0);
			CONV_L4(l41,l42,l43,l44,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L4C(l51,l52,l53,l54,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l51),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l52),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l53),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l54),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao5[n1p-16],r3);
			_mm_store_ps(&wao5[n1p-12],r2);
			_mm_store_ps(&wao5[n1p- 8],r1);
			_mm_store_ps(&wao5[n1p- 4],r0);
		}
		else if(n1t<21) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));

			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));

			CONV_L5(l41,l42,l43,l44,l45,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L5C(l51,l52,l53,l54,l55,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l41),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l42),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l43),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l44),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l45),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao4[n1p-20],r4);
			_mm_store_ps(&wao4[n1p-16],r3);
			_mm_store_ps(&wao4[n1p-12],r2);
			_mm_store_ps(&wao4[n1p- 8],r1);
			_mm_store_ps(&wao4[n1p- 4],r0);
			CONV_L5(l41,l42,l43,l44,l45,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L5C(l51,l52,l53,l54,l55,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l51),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l52),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l53),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l54),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l55),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao5[n1p-20],r4);
			_mm_store_ps(&wao5[n1p-16],r3);
			_mm_store_ps(&wao5[n1p-12],r2);
			_mm_store_ps(&wao5[n1p- 8],r1);
			_mm_store_ps(&wao5[n1p- 4],r0);
		}
		else if(n1t<25) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			l46=_mm_castps_si128(_mm_load_ps(wai4+20));

			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			l56=_mm_castps_si128(_mm_load_ps(wai5+20));

			CONV_L6(l41,l42,l43,l44,l45,l46,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L6C(l51,l52,l53,l54,l55,l56,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l41),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l42),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l43),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l44),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l45),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l46),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao4[n1p-24],r5);
			_mm_store_ps(&wao4[n1p-20],r4);
			_mm_store_ps(&wao4[n1p-16],r3);
			_mm_store_ps(&wao4[n1p-12],r2);
			_mm_store_ps(&wao4[n1p- 8],r1);
			_mm_store_ps(&wao4[n1p- 4],r0);
			CONV_L6(l41,l42,l43,l44,l45,l46,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L6C(l51,l52,l53,l54,l55,l56,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l51),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l52),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l53),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l54),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l55),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l56),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao5[n1p-24],r5);
			_mm_store_ps(&wao5[n1p-20],r4);
			_mm_store_ps(&wao5[n1p-16],r3);
			_mm_store_ps(&wao5[n1p-12],r2);
			_mm_store_ps(&wao5[n1p- 8],r1);
			_mm_store_ps(&wao5[n1p- 4],r0);
		}
		else if(n1t<29) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			l46=_mm_castps_si128(_mm_load_ps(wai4+20));
			l47=_mm_castps_si128(_mm_load_ps(wai4+24));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			l56=_mm_castps_si128(_mm_load_ps(wai5+20));
			l57=_mm_castps_si128(_mm_load_ps(wai5+24));
			CONV_L7(l41,l42,l43,l44,l45,l46,l47,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L7C(l51,l52,l53,l54,l55,l56,l57,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l41),cprecrsse);COUNT;r6=_mm_add_ps(r6,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l42),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l43),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l44),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l45),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l46),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l47),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao4[n1p-28],r6);
			_mm_store_ps(&wao4[n1p-24],r5);
			_mm_store_ps(&wao4[n1p-20],r4);
			_mm_store_ps(&wao4[n1p-16],r3);
			_mm_store_ps(&wao4[n1p-12],r2);
			_mm_store_ps(&wao4[n1p- 8],r1);
			_mm_store_ps(&wao4[n1p- 4],r0);
			CONV_L7(l41,l42,l43,l44,l45,l46,l47,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L7C(l51,l52,l53,l54,l55,l56,l57,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l51),cprecrsse);COUNT;r6=_mm_add_ps(r6,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l52),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l53),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l54),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l55),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l56),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l57),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao5[n1p-28],r6);
			_mm_store_ps(&wao5[n1p-24],r5);
			_mm_store_ps(&wao5[n1p-20],r4);
			_mm_store_ps(&wao5[n1p-16],r3);
			_mm_store_ps(&wao5[n1p-12],r2);
			_mm_store_ps(&wao5[n1p- 8],r1);
			_mm_store_ps(&wao5[n1p- 4],r0);
		}
		else {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			CONV_L8_PART1(l41,l42,l43,l44,l45,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r41,r42,r43,r44,r45,r46,r47);
			CONV_L8_PART1C(l51,l52,l53,l54,l55,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,r41,r42,r43,r44,r45,r46,r47);
			CONV_L8_PART1(l41,l42,l43,l44,l45,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,r51,r52,r53,r54,r55,r56,r57);
			CONV_L8_PART1C(l51,l52,l53,l54,l55,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r51,r52,r53,r54,r55,r56,r57);
			for(i1=28;i1<n1p-4;i1+=4) {
				l46=_mm_castps_si128(_mm_load_ps(wai4+i1-8));
				l56=_mm_castps_si128(_mm_load_ps(wai5+i1-8));
				CONV_L8_PART2(l44,l45,l46,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r40,r41,r42,r43,r44,r45,r46,r47);
				CONV_L8_PART2C(l54,l55,l56,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
					bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r40,r41,r42,r43,r44,r45,r46,r47);
				y=_mm_mul_ps(_mm_castsi128_ps(l41),cprecrsse);COUNT;r47=_mm_add_ps(r47,y);COUNT;
				_mm_store_ps(&wao4[i1-28],r47);
				CONV_L8_PART2(l44,l45,l46,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
					cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r50,r51,r52,r53,r54,r55,r56,r57);
				CONV_L8_PART2C(l54,l55,l56,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r50,r51,r52,r53,r54,r55,r56,r57);
				y=_mm_mul_ps(_mm_castsi128_ps(l51),cprecrsse);COUNT;r57=_mm_add_ps(r57,y);COUNT;
				_mm_store_ps(&wao5[i1-28],r57);
				l41=l42;l42=l43;l43=l44; //this is need because of cprecrsse 
				l44=l45;l45=l46;
				UP_L8(r40,r41,r42,r43,r44,r45,r46,r47);
				l51=l52;l52=l53;l53=l54; //this is need because of cprecrsse 
				l54=l55;l55=l56;
				UP_L8(r50,r51,r52,r53,r54,r55,r56,r57);
			}
			l46=_mm_castps_si128(_mm_load_ps(wai4+n1p-12));
			l47=_mm_castps_si128(_mm_load_ps(wai4+n1p- 8));
			l48=_mm_castps_si128(_mm_load_ps(wai4+n1p- 4));
			l56=_mm_castps_si128(_mm_load_ps(wai5+n1p-12));
			l57=_mm_castps_si128(_mm_load_ps(wai5+n1p- 8));
			l58=_mm_castps_si128(_mm_load_ps(wai5+n1p- 4));
			CONV_L8_PART3(l44,l45,l46,l47,l48,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r40,r41,r42,r43,r44,r45,r46,r47);
			CONV_L8_PART3C(l54,l55,l56,l57,l58,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r40,r41,r42,r43,r44,r45,r46,r47);
			y=_mm_mul_ps(_mm_castsi128_ps(l41),cprecrsse);COUNT;r47=_mm_add_ps(r47,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l42),cprecrsse);COUNT;r46=_mm_add_ps(r46,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l43),cprecrsse);COUNT;r45=_mm_add_ps(r45,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l44),cprecrsse);COUNT;r44=_mm_add_ps(r44,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l45),cprecrsse);COUNT;r43=_mm_add_ps(r43,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l46),cprecrsse);COUNT;r42=_mm_add_ps(r42,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l47),cprecrsse);COUNT;r41=_mm_add_ps(r41,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l48),cprecrsse);COUNT;r40=_mm_add_ps(r40,y);COUNT;
			_mm_store_ps(&wao4[n1p-32],r47);
			_mm_store_ps(&wao4[n1p-28],r46);
			_mm_store_ps(&wao4[n1p-24],r45);
			_mm_store_ps(&wao4[n1p-20],r44);
			_mm_store_ps(&wao4[n1p-16],r43);
			_mm_store_ps(&wao4[n1p-12],r42);
			_mm_store_ps(&wao4[n1p- 8],r41);
			_mm_store_ps(&wao4[n1p- 4],r40);
			CONV_L8_PART3(l44,l45,l46,l47,l48,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r50,r51,r52,r53,r54,r55,r56,r57);
			CONV_L8_PART3C(l54,l55,l56,l57,l58,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r50,r51,r52,r53,r54,r55,r56,r57);
			y=_mm_mul_ps(_mm_castsi128_ps(l51),cprecrsse);COUNT;r57=_mm_add_ps(r57,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l52),cprecrsse);COUNT;r56=_mm_add_ps(r56,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l53),cprecrsse);COUNT;r55=_mm_add_ps(r55,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l54),cprecrsse);COUNT;r54=_mm_add_ps(r54,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l55),cprecrsse);COUNT;r53=_mm_add_ps(r53,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l56),cprecrsse);COUNT;r52=_mm_add_ps(r52,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l57),cprecrsse);COUNT;r51=_mm_add_ps(r51,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l58),cprecrsse);COUNT;r50=_mm_add_ps(r50,y);COUNT;
			_mm_store_ps(&wao5[n1p-32],r57);
			_mm_store_ps(&wao5[n1p-28],r56);
			_mm_store_ps(&wao5[n1p-24],r55);
			_mm_store_ps(&wao5[n1p-20],r54);
			_mm_store_ps(&wao5[n1p-16],r53);
			_mm_store_ps(&wao5[n1p-12],r52);
			_mm_store_ps(&wao5[n1p- 8],r51);
			_mm_store_ps(&wao5[n1p- 4],r50);
		}
		//----------------------------------------------------------------------
		for(i1=0;i1<(n1p-4);i1+=4) {
			y_f[loffset1+(i1+0)*7]=wao4[i1+0];
			y_f[loffset2+(i1+0)*7]=wao5[i1+0];
			y_f[loffset1+(i1+1)*7]=wao4[i1+1];
			y_f[loffset2+(i1+1)*7]=wao5[i1+1];
			y_f[loffset1+(i1+2)*7]=wao4[i1+2];
			y_f[loffset2+(i1+2)*7]=wao5[i1+2];
			y_f[loffset1+(i1+3)*7]=wao4[i1+3];
			y_f[loffset2+(i1+3)*7]=wao5[i1+3];
		}
		for(;i1<n1t;i1++) {
			y_f[loffset1+i1*7]=wao4[i1];
			y_f[loffset2+i1*7]=wao5[i1];
		}
		} //end of if for n1t>0
		//----------------------------------------------------------------------
	} //end of loop over i2
	} //end of loop over i3
} // end of function convolut_ib_sse_fx_ssw_wsw_
//****************************************************************************************
void convolut_ib_sse_fx_sww_www_(int *n1,int *n2,int *n3,int *ibyz_f,float *cprecr,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n1p,nseg1,nseg2,n1t,istart;
	__m128 am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01,ap00,ap01;
	__m128 ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14;
	__m128 bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01;
	__m128 bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14;
	__m128 cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01;
	__m128 cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y,cprecrsse;
	__m128i lz;
	__m128i l61,l62,l63,l64,l65,l66,l67,l68;
	__m128i l71,l72,l73,l74,l75,l76,l77,l78;
	__m128 r60,r61,r62,r63,r64,r65,r66,r67;
	__m128 r70,r71,r72,r73,r74,r75,r76,r77;
	int i1,i2,i3,loffset1,loffset2;
	float wai6[1024] __attribute__ ((aligned (16)));
	float wai7[1024] __attribute__ ((aligned (16)));
	float wao6[1024] __attribute__ ((aligned (16)));
	float wao7[1024] __attribute__ ((aligned (16)));
	cprecrsse=_mm_set1_ps(cprecr[0]);
	FILTER_A; FILTER_B; FILTER_C; FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along X )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n1t,loffset1,loffset2,mm,n1p,wai6,wai7,wao6,wao7)\
           private(l61,l62,l63,l64,l65,l66,l67,l68,l71,l72,l73,l74,l75,l76,l77,l78)\
           private(r60,r61,r62,r63,r64,r65,r66,r67,r70,r71,r72,r73,r74,r75,r76,r77)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i3=*nfl3;i3<=*nfu3;i3++) {
	for(i2=*nfl2;i2<=*nfu2;i2++) {
		istart=ibyz_f[0+2*(i2+(*n2+1)*i3)];
		n1t=ibyz_f[1+2*(i2+(*n2+1)*i3)]-istart+1;
		if(n1t>0) {
		loffset1=5+7*(istart-*nfl1+nseg1*(i2-*nfl2+(i3-*nfl3)*nseg2));
		loffset2=6+7*(istart-*nfl1+nseg1*(i2-*nfl2+(i3-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i1=0;i1<n1t-4;i1+=4) {
			wai6[i1+0]=x_f[loffset1+(i1+0)*7];
			wai7[i1+0]=x_f[loffset2+(i1+0)*7];
			wai6[i1+1]=x_f[loffset1+(i1+1)*7];
			wai7[i1+1]=x_f[loffset2+(i1+1)*7];
			wai6[i1+2]=x_f[loffset1+(i1+2)*7];
			wai7[i1+2]=x_f[loffset2+(i1+2)*7];
			wai6[i1+3]=x_f[loffset1+(i1+3)*7];
			wai7[i1+3]=x_f[loffset2+(i1+3)*7];
		}
		mm=n1t%4;
		if(mm==0) {
			n1p=n1t;
			wai6[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai7[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai6[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai7[n1p-3]=x_f[loffset2+(n1p-3)*7];
			wai6[n1p-2]=x_f[loffset1+(n1p-2)*7];
			wai7[n1p-2]=x_f[loffset2+(n1p-2)*7];
			wai6[n1p-1]=x_f[loffset1+(n1p-1)*7];
			wai7[n1p-1]=x_f[loffset2+(n1p-1)*7];
		}
		else if(mm==1) {
			n1p=n1t+3;
			wai6[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai7[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai6[n1p-3]=0.0;
			wai7[n1p-3]=0.0;
			wai6[n1p-2]=0.0;
			wai7[n1p-2]=0.0;
			wai6[n1p-1]=0.0;
			wai7[n1p-1]=0.0;
		}
		else if(mm==2) {
			n1p=n1t+2;
			wai6[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai7[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai6[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai7[n1p-3]=x_f[loffset2+(n1p-3)*7];
			wai6[n1p-2]=0.0;
			wai7[n1p-2]=0.0;
			wai6[n1p-1]=0.0;
			wai7[n1p-1]=0.0;
		}
		else if(mm==3) {
			n1p=n1t+1;
			wai6[n1p-4]=x_f[loffset1+(n1p-4)*7];
			wai7[n1p-4]=x_f[loffset2+(n1p-4)*7];
			wai6[n1p-3]=x_f[loffset1+(n1p-3)*7];
			wai7[n1p-3]=x_f[loffset2+(n1p-3)*7];
			wai6[n1p-2]=x_f[loffset1+(n1p-2)*7];
			wai7[n1p-2]=x_f[loffset2+(n1p-2)*7];
			wai6[n1p-1]=0.0;
			wai7[n1p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n1t<5) {
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			CONV_L1(l61,am03,am02,am01,ap00,ap01,ap02,ap03);
			CONV_L1C(l71,bm03,bm02,bm01,bp00,bp01,bp02,bp03);
			y=_mm_mul_ps(_mm_castsi128_ps(l61),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao6[n1p- 4],r0);
			CONV_L1(l61,cm03,cm02,cm01,cp00,cp01,cp02,cp03);
			CONV_L1C(l71,em03,em02,em01,ep00,ep01,ep02,ep03);
			y=_mm_mul_ps(_mm_castsi128_ps(l71),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao7[n1p- 4],r0);
		}
		else if(n1t<9) {
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			CONV_L2(l61,l62,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			CONV_L2C(l71,l72,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07);
			y=_mm_mul_ps(_mm_castsi128_ps(l61),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l62),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao6[n1p- 8],r1);
			_mm_store_ps(&wao6[n1p- 4],r0);
			CONV_L2(l61,l62,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07);
			CONV_L2C(l71,l72,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			y=_mm_mul_ps(_mm_castsi128_ps(l71),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l72),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao7[n1p- 8],r1);
			_mm_store_ps(&wao7[n1p- 4],r0);
		}
		else if(n1t<13) {
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			CONV_L3(l61,l62,l63,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			CONV_L3C(l71,l72,l73,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11);
			y=_mm_mul_ps(_mm_castsi128_ps(l61),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l62),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l63),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao6[n1p-12],r2);
			_mm_store_ps(&wao6[n1p- 8],r1);
			_mm_store_ps(&wao6[n1p- 4],r0);
			CONV_L3(l61,l62,l63,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11);
			CONV_L3C(l71,l72,l73,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			y=_mm_mul_ps(_mm_castsi128_ps(l71),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l72),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l73),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao7[n1p-12],r2);
			_mm_store_ps(&wao7[n1p- 8],r1);
			_mm_store_ps(&wao7[n1p- 4],r0);
		}
		else if(n1t<17) {
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			CONV_L4(l61,l62,l63,l64,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L4C(l71,l72,l73,l74,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l61),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l62),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l63),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l64),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao6[n1p-16],r3);
			_mm_store_ps(&wao6[n1p-12],r2);
			_mm_store_ps(&wao6[n1p- 8],r1);
			_mm_store_ps(&wao6[n1p- 4],r0);
			CONV_L4(l61,l62,l63,l64,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L4C(l71,l72,l73,l74,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l71),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l72),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l73),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l74),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao7[n1p-16],r3);
			_mm_store_ps(&wao7[n1p-12],r2);
			_mm_store_ps(&wao7[n1p- 8],r1);
			_mm_store_ps(&wao7[n1p- 4],r0);
		}
		else if(n1t<21) {
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			CONV_L5(l61,l62,l63,l64,l65,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L5C(l71,l72,l73,l74,l75,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l61),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l62),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l63),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l64),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l65),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao6[n1p-20],r4);
			_mm_store_ps(&wao6[n1p-16],r3);
			_mm_store_ps(&wao6[n1p-12],r2);
			_mm_store_ps(&wao6[n1p- 8],r1);
			_mm_store_ps(&wao6[n1p- 4],r0);
			CONV_L5(l61,l62,l63,l64,l65,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L5C(l71,l72,l73,l74,l75,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l71),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l72),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l73),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l74),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l75),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao7[n1p-20],r4);
			_mm_store_ps(&wao7[n1p-16],r3);
			_mm_store_ps(&wao7[n1p-12],r2);
			_mm_store_ps(&wao7[n1p- 8],r1);
			_mm_store_ps(&wao7[n1p- 4],r0);
		}
		else if(n1t<25) {
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			l66=_mm_castps_si128(_mm_load_ps(wai6+20));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			l76=_mm_castps_si128(_mm_load_ps(wai7+20));
			CONV_L6(l61,l62,l63,l64,l65,l66,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L6C(l71,l72,l73,l74,l75,l76,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l61),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l62),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l63),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l64),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l65),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l66),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao6[n1p-24],r5);
			_mm_store_ps(&wao6[n1p-20],r4);
			_mm_store_ps(&wao6[n1p-16],r3);
			_mm_store_ps(&wao6[n1p-12],r2);
			_mm_store_ps(&wao6[n1p- 8],r1);
			_mm_store_ps(&wao6[n1p- 4],r0);
			CONV_L6(l61,l62,l63,l64,l65,l66,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L6C(l71,l72,l73,l74,l75,l76,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l71),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l72),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l73),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l74),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l75),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l76),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao7[n1p-24],r5);
			_mm_store_ps(&wao7[n1p-20],r4);
			_mm_store_ps(&wao7[n1p-16],r3);
			_mm_store_ps(&wao7[n1p-12],r2);
			_mm_store_ps(&wao7[n1p- 8],r1);
			_mm_store_ps(&wao7[n1p- 4],r0);
		}
		else if(n1t<29) {
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			l66=_mm_castps_si128(_mm_load_ps(wai6+20));
			l67=_mm_castps_si128(_mm_load_ps(wai6+24));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			l76=_mm_castps_si128(_mm_load_ps(wai7+20));
			l77=_mm_castps_si128(_mm_load_ps(wai7+24));
			CONV_L7(l61,l62,l63,l64,l65,l66,l67,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L7C(l71,l72,l73,l74,l75,l76,l77,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			y=_mm_mul_ps(_mm_castsi128_ps(l61),cprecrsse);COUNT;r6=_mm_add_ps(r6,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l62),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l63),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l64),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l65),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l66),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l67),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao6[n1p-28],r6);
			_mm_store_ps(&wao6[n1p-24],r5);
			_mm_store_ps(&wao6[n1p-20],r4);
			_mm_store_ps(&wao6[n1p-16],r3);
			_mm_store_ps(&wao6[n1p-12],r2);
			_mm_store_ps(&wao6[n1p- 8],r1);
			_mm_store_ps(&wao6[n1p- 4],r0);
			CONV_L7(l61,l62,l63,l64,l65,l66,l67,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L7C(l71,l72,l73,l74,l75,l76,l77,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			y=_mm_mul_ps(_mm_castsi128_ps(l71),cprecrsse);COUNT;r6=_mm_add_ps(r6,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l72),cprecrsse);COUNT;r5=_mm_add_ps(r5,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l73),cprecrsse);COUNT;r4=_mm_add_ps(r4,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l74),cprecrsse);COUNT;r3=_mm_add_ps(r3,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l75),cprecrsse);COUNT;r2=_mm_add_ps(r2,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l76),cprecrsse);COUNT;r1=_mm_add_ps(r1,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l77),cprecrsse);COUNT;r0=_mm_add_ps(r0,y);COUNT;
			_mm_store_ps(&wao7[n1p-28],r6);
			_mm_store_ps(&wao7[n1p-24],r5);
			_mm_store_ps(&wao7[n1p-20],r4);
			_mm_store_ps(&wao7[n1p-16],r3);
			_mm_store_ps(&wao7[n1p-12],r2);
			_mm_store_ps(&wao7[n1p- 8],r1);
			_mm_store_ps(&wao7[n1p- 4],r0);
		}
		else {
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			CONV_L8_PART1(l61,l62,l63,l64,l65,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r61,r62,r63,r64,r65,r66,r67);
			CONV_L8_PART1C(l71,l72,l73,l74,l75,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,r61,r62,r63,r64,r65,r66,r67);

			CONV_L8_PART1(l61,l62,l63,l64,l65,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,r71,r72,r73,r74,r75,r76,r77);
			CONV_L8_PART1C(l71,l72,l73,l74,l75,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r71,r72,r73,r74,r75,r76,r77);
			for(i1=28;i1<n1p-4;i1+=4) {
				l66=_mm_castps_si128(_mm_load_ps(wai6+i1-8));
				l76=_mm_castps_si128(_mm_load_ps(wai7+i1-8));
				CONV_L8_PART2(l64,l65,l66,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r60,r61,r62,r63,r64,r65,r66,r67);
				CONV_L8_PART2C(l74,l75,l76,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
					bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r60,r61,r62,r63,r64,r65,r66,r67);
				y=_mm_mul_ps(_mm_castsi128_ps(l61),cprecrsse);COUNT;r67=_mm_add_ps(r67,y);COUNT;
				_mm_store_ps(&wao6[i1-28],r67);
				CONV_L8_PART2(l64,l65,l66,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
					cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r70,r71,r72,r73,r74,r75,r76,r77);
				CONV_L8_PART2C(l74,l75,l76,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r70,r71,r72,r73,r74,r75,r76,r77);
				y=_mm_mul_ps(_mm_castsi128_ps(l71),cprecrsse);COUNT;r77=_mm_add_ps(r77,y);COUNT;
				_mm_store_ps(&wao7[i1-28],r77);
				l61=l62;l62=l63;l63=l64; //this is need because of cprecrsse 
				l64=l65;l65=l66;
				UP_L8(r60,r61,r62,r63,r64,r65,r66,r67);
				l71=l72;l72=l73;l73=l74; //this is need because of cprecrsse 
				l74=l75;l75=l76;
				UP_L8(r70,r71,r72,r73,r74,r75,r76,r77);
			}
			l66=_mm_castps_si128(_mm_load_ps(wai6+n1p-12));
			l67=_mm_castps_si128(_mm_load_ps(wai6+n1p- 8));
			l68=_mm_castps_si128(_mm_load_ps(wai6+n1p- 4));
			l76=_mm_castps_si128(_mm_load_ps(wai7+n1p-12));
			l77=_mm_castps_si128(_mm_load_ps(wai7+n1p- 8));
			l78=_mm_castps_si128(_mm_load_ps(wai7+n1p- 4));
			CONV_L8_PART3(l64,l65,l66,l67,l68,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r60,r61,r62,r63,r64,r65,r66,r67);
			CONV_L8_PART3C(l74,l75,l76,l77,l78,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r60,r61,r62,r63,r64,r65,r66,r67);
			y=_mm_mul_ps(_mm_castsi128_ps(l61),cprecrsse);COUNT;r67=_mm_add_ps(r67,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l62),cprecrsse);COUNT;r66=_mm_add_ps(r66,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l63),cprecrsse);COUNT;r65=_mm_add_ps(r65,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l64),cprecrsse);COUNT;r64=_mm_add_ps(r64,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l65),cprecrsse);COUNT;r63=_mm_add_ps(r63,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l66),cprecrsse);COUNT;r62=_mm_add_ps(r62,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l67),cprecrsse);COUNT;r61=_mm_add_ps(r61,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l68),cprecrsse);COUNT;r60=_mm_add_ps(r60,y);COUNT;
			_mm_store_ps(&wao6[n1p-32],r67);
			_mm_store_ps(&wao6[n1p-28],r66);
			_mm_store_ps(&wao6[n1p-24],r65);
			_mm_store_ps(&wao6[n1p-20],r64);
			_mm_store_ps(&wao6[n1p-16],r63);
			_mm_store_ps(&wao6[n1p-12],r62);
			_mm_store_ps(&wao6[n1p- 8],r61);
			_mm_store_ps(&wao6[n1p- 4],r60);
			CONV_L8_PART3(l64,l65,l66,l67,l68,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r70,r71,r72,r73,r74,r75,r76,r77);
			CONV_L8_PART3C(l74,l75,l76,l77,l78,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r70,r71,r72,r73,r74,r75,r76,r77);
			y=_mm_mul_ps(_mm_castsi128_ps(l71),cprecrsse);COUNT;r77=_mm_add_ps(r77,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l72),cprecrsse);COUNT;r76=_mm_add_ps(r76,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l73),cprecrsse);COUNT;r75=_mm_add_ps(r75,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l74),cprecrsse);COUNT;r74=_mm_add_ps(r74,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l75),cprecrsse);COUNT;r73=_mm_add_ps(r73,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l76),cprecrsse);COUNT;r72=_mm_add_ps(r72,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l77),cprecrsse);COUNT;r71=_mm_add_ps(r71,y);COUNT;
			y=_mm_mul_ps(_mm_castsi128_ps(l78),cprecrsse);COUNT;r70=_mm_add_ps(r70,y);COUNT;
			_mm_store_ps(&wao7[n1p-32],r77);
			_mm_store_ps(&wao7[n1p-28],r76);
			_mm_store_ps(&wao7[n1p-24],r75);
			_mm_store_ps(&wao7[n1p-20],r74);
			_mm_store_ps(&wao7[n1p-16],r73);
			_mm_store_ps(&wao7[n1p-12],r72);
			_mm_store_ps(&wao7[n1p- 8],r71);
			_mm_store_ps(&wao7[n1p- 4],r70);
		}
		//----------------------------------------------------------------------
		for(i1=0;i1<(n1p-4);i1+=4) {
			y_f[loffset1+(i1+0)*7]=wao6[i1+0];
			y_f[loffset2+(i1+0)*7]=wao7[i1+0];
			y_f[loffset1+(i1+1)*7]=wao6[i1+1];
			y_f[loffset2+(i1+1)*7]=wao7[i1+1];
			y_f[loffset1+(i1+2)*7]=wao6[i1+2];
			y_f[loffset2+(i1+2)*7]=wao7[i1+2];
			y_f[loffset1+(i1+3)*7]=wao6[i1+3];
			y_f[loffset2+(i1+3)*7]=wao7[i1+3];
		}
		for(;i1<n1t;i1++) {
			y_f[loffset1+i1*7]=wao6[i1];
			y_f[loffset2+i1*7]=wao7[i1];
		}
		} //end of if for n1t>0
		//----------------------------------------------------------------------
	} //end of loop over i2
	} //end of loop over i3
} // end of function convolut_ib_sse_fx_sww_www_
//****************************************************************************************
void convolut_ib_sse_fy_sws_(int *n1,int *n2,int *n3,int *ibxz_f,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n2p,nseg1,nseg2,nseg3,n2t,istart;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y;
	__m128i lz;
	__m128i l21,l22,l23,l24,l25,l26,l27,l28;
	__m128 r20,r21,r22,r23,r24,r25,r26,r27;
	int i1,i2,i3,loffset1;
	float wai2[1024] __attribute__ ((aligned (16)));
	float wao2[1024] __attribute__ ((aligned (16)));
	//----------------------------------
	FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	nseg3=nseg1*7;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along Y )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n2t,loffset1,mm,n2p,wai2,wao2)\
           private(l21,l22,l23,l24,l25,l26,l27,l28,r20,r21,r22,r23,r24,r25,r26,r27)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i3=*nfl3;i3<=*nfu3;i3++) {
	for(i1=*nfl1;i1<=*nfu1;i1++) {
		istart=ibxz_f[0+2*(i1+(*n1+1)*i3)];
		n2t=ibxz_f[1+2*(i1+(*n1+1)*i3)]-istart+1;
		if(n2t>0) {
		loffset1=1+7*(i1-*nfl1+nseg1*(istart-*nfl2+(i3-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i2=0;i2<n2t-4;i2+=4) {
			wai2[i2+0]=x_f[loffset1+(i2+0)*nseg3];
			wai2[i2+1]=x_f[loffset1+(i2+1)*nseg3];
			wai2[i2+2]=x_f[loffset1+(i2+2)*nseg3];
			wai2[i2+3]=x_f[loffset1+(i2+3)*nseg3];
		}
		mm=n2t%4;
		if(mm==0) {
			n2p=n2t;
			wai2[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai2[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai2[n2p-2]=x_f[loffset1+(n2p-2)*nseg3];
			wai2[n2p-1]=x_f[loffset1+(n2p-1)*nseg3];
		}
		else if(mm==1) {
			n2p=n2t+3;
			wai2[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai2[n2p-3]=0.0;
			wai2[n2p-2]=0.0;
			wai2[n2p-1]=0.0;
		}
		else if(mm==2) {
			n2p=n2t+2;
			wai2[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai2[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai2[n2p-2]=0.0;
			wai2[n2p-1]=0.0;
		}
		else if(mm==3) {
			n2p=n2t+1;
			wai2[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai2[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai2[n2p-2]=x_f[loffset1+(n2p-2)*nseg3];
			wai2[n2p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n2t<5) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			CONV_L1(l21,em03,em02,em01,ep00,ep01,ep02,ep03);
			_mm_store_ps(&wao2[n2p- 4],r0);
		}
		else if(n2t<9) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			CONV_L2(l21,l22,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			_mm_store_ps(&wao2[n2p- 8],r1);
			_mm_store_ps(&wao2[n2p- 4],r0);
		}
		else if(n2t<13) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			CONV_L3(l21,l22,l23,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			_mm_store_ps(&wao2[n2p-12],r2);
			_mm_store_ps(&wao2[n2p- 8],r1);
			_mm_store_ps(&wao2[n2p- 4],r0);
		}
		else if(n2t<17) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			CONV_L4(l21,l22,l23,l24,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao2[n2p-16],r3);
			_mm_store_ps(&wao2[n2p-12],r2);
			_mm_store_ps(&wao2[n2p- 8],r1);
			_mm_store_ps(&wao2[n2p- 4],r0);
		}
		else if(n2t<21) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));
			CONV_L5(l21,l22,l23,l24,l25,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao2[n2p-20],r4);
			_mm_store_ps(&wao2[n2p-16],r3);
			_mm_store_ps(&wao2[n2p-12],r2);
			_mm_store_ps(&wao2[n2p- 8],r1);
			_mm_store_ps(&wao2[n2p- 4],r0);
		}
		else if(n2t<25) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));
			l26=_mm_castps_si128(_mm_load_ps(wai2+20));
			CONV_L6(l21,l22,l23,l24,l25,l26,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao2[n2p-24],r5);
			_mm_store_ps(&wao2[n2p-20],r4);
			_mm_store_ps(&wao2[n2p-16],r3);
			_mm_store_ps(&wao2[n2p-12],r2);
			_mm_store_ps(&wao2[n2p- 8],r1);
			_mm_store_ps(&wao2[n2p- 4],r0);
		}
		else if(n2t<29) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));
			l26=_mm_castps_si128(_mm_load_ps(wai2+20));
			l27=_mm_castps_si128(_mm_load_ps(wai2+24));
			CONV_L7(l21,l22,l23,l24,l25,l26,l27,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao2[n2p-28],r6);
			_mm_store_ps(&wao2[n2p-24],r5);
			_mm_store_ps(&wao2[n2p-20],r4);
			_mm_store_ps(&wao2[n2p-16],r3);
			_mm_store_ps(&wao2[n2p-12],r2);
			_mm_store_ps(&wao2[n2p- 8],r1);
			_mm_store_ps(&wao2[n2p- 4],r0);
		}
		else {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));
			CONV_L8_PART1(l21,l22,l23,l24,l25,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r21,r22,r23,r24,r25,r26,r27);
			for(i2=28;i2<n2p-4;i2+=4) {
				l26=_mm_castps_si128(_mm_load_ps(wai2+i2-8));
				CONV_L8_PART2(l24,l25,l26,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r20,r21,r22,r23,r24,r25,r26,r27);
				_mm_store_ps(&wao2[i2-28],r27);
				l24=l25;l25=l26;
				UP_L8(r20,r21,r22,r23,r24,r25,r26,r27);
			}
			l26=_mm_castps_si128(_mm_load_ps(wai2+n2p-12));
			l27=_mm_castps_si128(_mm_load_ps(wai2+n2p- 8));
			l28=_mm_castps_si128(_mm_load_ps(wai2+n2p- 4));
			CONV_L8_PART3(l24,l25,l26,l27,l28,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r20,r21,r22,r23,r24,r25,r26,r27);
			_mm_store_ps(&wao2[n2p-32],r27);
			_mm_store_ps(&wao2[n2p-28],r26);
			_mm_store_ps(&wao2[n2p-24],r25);
			_mm_store_ps(&wao2[n2p-20],r24);
			_mm_store_ps(&wao2[n2p-16],r23);
			_mm_store_ps(&wao2[n2p-12],r22);
			_mm_store_ps(&wao2[n2p- 8],r21);
			_mm_store_ps(&wao2[n2p- 4],r20);
		}
		//----------------------------------------------------------------------
		for(i2=0;i2<(n2p-4);i2+=4) {
			y_f[loffset1+(i2+0)*nseg3]=y_f[loffset1+(i2+0)*nseg3]+wao2[i2+0];
			y_f[loffset1+(i2+1)*nseg3]=y_f[loffset1+(i2+1)*nseg3]+wao2[i2+1];
			y_f[loffset1+(i2+2)*nseg3]=y_f[loffset1+(i2+2)*nseg3]+wao2[i2+2];
			y_f[loffset1+(i2+3)*nseg3]=y_f[loffset1+(i2+3)*nseg3]+wao2[i2+3];
		}
		for(;i2<n2t;i2++) {
			y_f[loffset1+i2*nseg3]=y_f[loffset1+i2*nseg3]+wao2[i2];
		}
		} //end of if for n2t>0
		//----------------------------------------------------------------------
	} //end of loop over i1
	} //end of loop over i3
} // end of function convolut_ib_sse_fy_sws_
//****************************************************************************************
void convolut_ib_sse_fy_wss_wws_(int *n1,int *n2,int *n3,int *ibxz_f,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n2p,nseg1,nseg2,nseg3,n2t,istart;
	__m128 am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01,ap00,ap01;
	__m128 ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14;
	__m128 bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01;
	__m128 bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14;
	__m128 cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01;
	__m128 cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y;
	__m128i lz;
	__m128i l11,l12,l13,l14,l15,l16,l17,l18;
	__m128i l31,l32,l33,l34,l35,l36,l37,l38;
	__m128 r10,r11,r12,r13,r14,r15,r16,r17;
	__m128 r30,r31,r32,r33,r34,r35,r36,r37;
	int i1,i2,i3,loffset1,loffset2;
	float wai1[1024] __attribute__ ((aligned (16)));
	float wai3[1024] __attribute__ ((aligned (16)));
	float wao1[1024] __attribute__ ((aligned (16)));
	float wao3[1024] __attribute__ ((aligned (16)));
	//----------------------------------
	FILTER_A; FILTER_B; FILTER_C; FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	nseg3=nseg1*7;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along Y )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n2t,loffset1,loffset2,mm,n2p,wai1,wai3,wao1,wao3)\
           private(l11,l12,l13,l14,l15,l16,l17,l18,l31,l32,l33,l34,l35,l36,l37,l38)\
           private(r10,r11,r12,r13,r14,r15,r16,r17,r30,r31,r32,r33,r34,r35,r36,r37)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i3=*nfl3;i3<=*nfu3;i3++) {
	for(i1=*nfl1;i1<=*nfu1;i1++) {
		istart=ibxz_f[0+2*(i1+(*n1+1)*i3)];
		n2t=ibxz_f[1+2*(i1+(*n1+1)*i3)]-istart+1;
		if(n2t>0) {
		loffset1=0+7*(i1-*nfl1+nseg1*(istart-*nfl2+(i3-*nfl3)*nseg2));
		loffset2=2+7*(i1-*nfl1+nseg1*(istart-*nfl2+(i3-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i2=0;i2<n2t-4;i2+=4) {
			wai1[i2+0]=x_f[loffset1+(i2+0)*nseg3];
			wai3[i2+0]=x_f[loffset2+(i2+0)*nseg3];
			wai1[i2+1]=x_f[loffset1+(i2+1)*nseg3];
			wai3[i2+1]=x_f[loffset2+(i2+1)*nseg3];
			wai1[i2+2]=x_f[loffset1+(i2+2)*nseg3];
			wai3[i2+2]=x_f[loffset2+(i2+2)*nseg3];
			wai1[i2+3]=x_f[loffset1+(i2+3)*nseg3];
			wai3[i2+3]=x_f[loffset2+(i2+3)*nseg3];
		}
		mm=n2t%4;
		if(mm==0) {
			n2p=n2t;
			wai1[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai3[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai1[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai3[n2p-3]=x_f[loffset2+(n2p-3)*nseg3];
			wai1[n2p-2]=x_f[loffset1+(n2p-2)*nseg3];
			wai3[n2p-2]=x_f[loffset2+(n2p-2)*nseg3];
			wai1[n2p-1]=x_f[loffset1+(n2p-1)*nseg3];
			wai3[n2p-1]=x_f[loffset2+(n2p-1)*nseg3];
		}
		else if(mm==1) {
			n2p=n2t+3;
			wai1[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai3[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai1[n2p-3]=0.0;
			wai3[n2p-3]=0.0;
			wai1[n2p-2]=0.0;
			wai3[n2p-2]=0.0;
			wai1[n2p-1]=0.0;
			wai3[n2p-1]=0.0;
		}
		else if(mm==2) {
			n2p=n2t+2;
			wai1[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai3[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai1[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai3[n2p-3]=x_f[loffset2+(n2p-3)*nseg3];
			wai1[n2p-2]=0.0;
			wai3[n2p-2]=0.0;
			wai1[n2p-1]=0.0;
			wai3[n2p-1]=0.0;
		}
		else if(mm==3) {
			n2p=n2t+1;
			wai1[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai3[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai1[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai3[n2p-3]=x_f[loffset2+(n2p-3)*nseg3];
			wai1[n2p-2]=x_f[loffset1+(n2p-2)*nseg3];
			wai3[n2p-2]=x_f[loffset2+(n2p-2)*nseg3];
			wai1[n2p-1]=0.0;
			wai3[n2p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n2t<5) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			CONV_L1(l11,am03,am02,am01,ap00,ap01,ap02,ap03);
			CONV_L1C(l31,bm03,bm02,bm01,bp00,bp01,bp02,bp03);
			_mm_store_ps(&wao1[n2p- 4],r0);
			CONV_L1(l11,cm03,cm02,cm01,cp00,cp01,cp02,cp03);
			CONV_L1C(l31,em03,em02,em01,ep00,ep01,ep02,ep03);
			_mm_store_ps(&wao3[n2p- 4],r0);
		}
		else if(n2t<9) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			CONV_L2(l11,l12,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			CONV_L2C(l31,l32,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07);
			_mm_store_ps(&wao1[n2p- 8],r1);
			_mm_store_ps(&wao1[n2p- 4],r0);
			CONV_L2(l11,l12,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07);
			CONV_L2C(l31,l32,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			_mm_store_ps(&wao3[n2p- 8],r1);
			_mm_store_ps(&wao3[n2p- 4],r0);
		}
		else if(n2t<13) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			CONV_L3(l11,l12,l13,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			CONV_L3C(l31,l32,l33,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11);
			_mm_store_ps(&wao1[n2p-12],r2);
			_mm_store_ps(&wao1[n2p- 8],r1);
			_mm_store_ps(&wao1[n2p- 4],r0);
			CONV_L3(l11,l12,l13,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11);
			CONV_L3C(l31,l32,l33,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			_mm_store_ps(&wao3[n2p-12],r2);
			_mm_store_ps(&wao3[n2p- 8],r1);
			_mm_store_ps(&wao3[n2p- 4],r0);
		}
		else if(n2t<17) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			CONV_L4(l11,l12,l13,l14,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L4C(l31,l32,l33,l34,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao1[n2p-16],r3);
			_mm_store_ps(&wao1[n2p-12],r2);
			_mm_store_ps(&wao1[n2p- 8],r1);
			_mm_store_ps(&wao1[n2p- 4],r0);
			CONV_L4(l11,l12,l13,l14,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L4C(l31,l32,l33,l34,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao3[n2p-16],r3);
			_mm_store_ps(&wao3[n2p-12],r2);
			_mm_store_ps(&wao3[n2p- 8],r1);
			_mm_store_ps(&wao3[n2p- 4],r0);
		}
		else if(n2t<21) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));
			CONV_L5(l11,l12,l13,l14,l15,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L5C(l31,l32,l33,l34,l35,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao1[n2p-20],r4);
			_mm_store_ps(&wao1[n2p-16],r3);
			_mm_store_ps(&wao1[n2p-12],r2);
			_mm_store_ps(&wao1[n2p- 8],r1);
			_mm_store_ps(&wao1[n2p- 4],r0);
			CONV_L5(l11,l12,l13,l14,l15,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L5C(l31,l32,l33,l34,l35,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao3[n2p-20],r4);
			_mm_store_ps(&wao3[n2p-16],r3);
			_mm_store_ps(&wao3[n2p-12],r2);
			_mm_store_ps(&wao3[n2p- 8],r1);
			_mm_store_ps(&wao3[n2p- 4],r0);
		}
		else if(n2t<25) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));
			l16=_mm_castps_si128(_mm_load_ps(wai1+20));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));
			l36=_mm_castps_si128(_mm_load_ps(wai3+20));
			CONV_L6(l11,l12,l13,l14,l15,l16,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L6C(l31,l32,l33,l34,l35,l36,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao1[n2p-24],r5);
			_mm_store_ps(&wao1[n2p-20],r4);
			_mm_store_ps(&wao1[n2p-16],r3);
			_mm_store_ps(&wao1[n2p-12],r2);
			_mm_store_ps(&wao1[n2p- 8],r1);
			_mm_store_ps(&wao1[n2p- 4],r0);
			CONV_L6(l11,l12,l13,l14,l15,l16,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L6C(l31,l32,l33,l34,l35,l36,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao3[n2p-24],r5);
			_mm_store_ps(&wao3[n2p-20],r4);
			_mm_store_ps(&wao3[n2p-16],r3);
			_mm_store_ps(&wao3[n2p-12],r2);
			_mm_store_ps(&wao3[n2p- 8],r1);
			_mm_store_ps(&wao3[n2p- 4],r0);
		}
		else if(n2t<29) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));
			l16=_mm_castps_si128(_mm_load_ps(wai1+20));
			l17=_mm_castps_si128(_mm_load_ps(wai1+24));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));
			l36=_mm_castps_si128(_mm_load_ps(wai3+20));
			l37=_mm_castps_si128(_mm_load_ps(wai3+24));
			CONV_L7(l11,l12,l13,l14,l15,l16,l17,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L7C(l31,l32,l33,l34,l35,l36,l37,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao1[n2p-28],r6);
			_mm_store_ps(&wao1[n2p-24],r5);
			_mm_store_ps(&wao1[n2p-20],r4);
			_mm_store_ps(&wao1[n2p-16],r3);
			_mm_store_ps(&wao1[n2p-12],r2);
			_mm_store_ps(&wao1[n2p- 8],r1);
			_mm_store_ps(&wao1[n2p- 4],r0);
			CONV_L7(l11,l12,l13,l14,l15,l16,l17,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L7C(l31,l32,l33,l34,l35,l36,l37,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao3[n2p-28],r6);
			_mm_store_ps(&wao3[n2p-24],r5);
			_mm_store_ps(&wao3[n2p-20],r4);
			_mm_store_ps(&wao3[n2p-16],r3);
			_mm_store_ps(&wao3[n2p-12],r2);
			_mm_store_ps(&wao3[n2p- 8],r1);
			_mm_store_ps(&wao3[n2p- 4],r0);
		}
		else {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));
			CONV_L8_PART1(l11,l12,l13,l14,l15,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r11,r12,r13,r14,r15,r16,r17);
			CONV_L8_PART1C(l31,l32,l33,l34,l35,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,r11,r12,r13,r14,r15,r16,r17);
			CONV_L8_PART1(l11,l12,l13,l14,l15,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,r31,r32,r33,r34,r35,r36,r37);
			CONV_L8_PART1C(l31,l32,l33,l34,l35,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r31,r32,r33,r34,r35,r36,r37);
			for(i2=28;i2<n2p-4;i2+=4) {
				l16=_mm_castps_si128(_mm_load_ps(wai1+i2-8));
				l36=_mm_castps_si128(_mm_load_ps(wai3+i2-8));
				CONV_L8_PART2(l14,l15,l16,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r10,r11,r12,r13,r14,r15,r16,r17);
				CONV_L8_PART2C(l34,l35,l36,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
					bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r10,r11,r12,r13,r14,r15,r16,r17);
				_mm_store_ps(&wao1[i2-28],r17);
				CONV_L8_PART2(l14,l15,l16,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
					cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r30,r31,r32,r33,r34,r35,r36,r37);
				CONV_L8_PART2C(l34,l35,l36,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r30,r31,r32,r33,r34,r35,r36,r37);
				_mm_store_ps(&wao3[i2-28],r37);
				l14=l15;l15=l16;
				UP_L8(r10,r11,r12,r13,r14,r15,r16,r17);
				l34=l35;l35=l36;
				UP_L8(r30,r31,r32,r33,r34,r35,r36,r37);
			}
			l16=_mm_castps_si128(_mm_load_ps(wai1+n2p-12));
			l17=_mm_castps_si128(_mm_load_ps(wai1+n2p- 8));
			l18=_mm_castps_si128(_mm_load_ps(wai1+n2p- 4));
			l36=_mm_castps_si128(_mm_load_ps(wai3+n2p-12));
			l37=_mm_castps_si128(_mm_load_ps(wai3+n2p- 8));
			l38=_mm_castps_si128(_mm_load_ps(wai3+n2p- 4));
			CONV_L8_PART3(l14,l15,l16,l17,l18,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r10,r11,r12,r13,r14,r15,r16,r17);
			CONV_L8_PART3C(l34,l35,l36,l37,l38,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r10,r11,r12,r13,r14,r15,r16,r17);
			_mm_store_ps(&wao1[n2p-32],r17);
			_mm_store_ps(&wao1[n2p-28],r16);
			_mm_store_ps(&wao1[n2p-24],r15);
			_mm_store_ps(&wao1[n2p-20],r14);
			_mm_store_ps(&wao1[n2p-16],r13);
			_mm_store_ps(&wao1[n2p-12],r12);
			_mm_store_ps(&wao1[n2p- 8],r11);
			_mm_store_ps(&wao1[n2p- 4],r10);
			CONV_L8_PART3(l14,l15,l16,l17,l18,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r30,r31,r32,r33,r34,r35,r36,r37);
			CONV_L8_PART3C(l34,l35,l36,l37,l38,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r30,r31,r32,r33,r34,r35,r36,r37);
			_mm_store_ps(&wao3[n2p-32],r37);
			_mm_store_ps(&wao3[n2p-28],r36);
			_mm_store_ps(&wao3[n2p-24],r35);
			_mm_store_ps(&wao3[n2p-20],r34);
			_mm_store_ps(&wao3[n2p-16],r33);
			_mm_store_ps(&wao3[n2p-12],r32);
			_mm_store_ps(&wao3[n2p- 8],r31);
			_mm_store_ps(&wao3[n2p- 4],r30);
		}
		//----------------------------------------------------------------------
		for(i2=0;i2<(n2p-4);i2+=4) {
			y_f[loffset1+(i2+0)*nseg3]=y_f[loffset1+(i2+0)*nseg3]+wao1[i2+0];
			y_f[loffset2+(i2+0)*nseg3]=y_f[loffset2+(i2+0)*nseg3]+wao3[i2+0];
			y_f[loffset1+(i2+1)*nseg3]=y_f[loffset1+(i2+1)*nseg3]+wao1[i2+1];
			y_f[loffset2+(i2+1)*nseg3]=y_f[loffset2+(i2+1)*nseg3]+wao3[i2+1];
			y_f[loffset1+(i2+2)*nseg3]=y_f[loffset1+(i2+2)*nseg3]+wao1[i2+2];
			y_f[loffset2+(i2+2)*nseg3]=y_f[loffset2+(i2+2)*nseg3]+wao3[i2+2];
			y_f[loffset1+(i2+3)*nseg3]=y_f[loffset1+(i2+3)*nseg3]+wao1[i2+3];
			y_f[loffset2+(i2+3)*nseg3]=y_f[loffset2+(i2+3)*nseg3]+wao3[i2+3];
		}
		for(;i2<n2t;i2++) {
			y_f[loffset1+i2*nseg3]=y_f[loffset1+i2*nseg3]+wao1[i2];
			y_f[loffset2+i2*nseg3]=y_f[loffset2+i2*nseg3]+wao3[i2];
		}
		} //end of if for n2t>0
		//----------------------------------------------------------------------
	} //end of loop over i1
	} //end of loop over i3
} // end of function convolut_ib_sse_fy_wss_wws_
//****************************************************************************************
void convolut_ib_sse_fy_ssw_sww_(int *n1,int *n2,int *n3,int *ibxz_f,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n2p,nseg1,nseg2,nseg3,n2t,istart;
	__m128 am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01,ap00,ap01;
	__m128 ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14;
	__m128 bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01;
	__m128 bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14;
	__m128 cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01;
	__m128 cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y;
	__m128i lz;
	__m128i l41,l42,l43,l44,l45,l46,l47,l48;
	__m128i l61,l62,l63,l64,l65,l66,l67,l68;
	__m128 r40,r41,r42,r43,r44,r45,r46,r47;
	__m128 r60,r61,r62,r63,r64,r65,r66,r67;
	int i1,i2,i3,loffset1,loffset2;
	float wai4[1024] __attribute__ ((aligned (16)));
	float wai6[1024] __attribute__ ((aligned (16)));
	float wao4[1024] __attribute__ ((aligned (16)));
	float wao6[1024] __attribute__ ((aligned (16)));
	//----------------------------------
	FILTER_A; FILTER_B; FILTER_C; FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	nseg3=nseg1*7;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along Y )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n2t,loffset1,loffset2,mm,n2p,wai4,wai6,wao4,wao6)\
           private(l41,l42,l43,l44,l45,l46,l47,l48,l61,l62,l63,l64,l65,l66,l67,l68)\
           private(r40,r41,r42,r43,r44,r45,r46,r47,r60,r61,r62,r63,r64,r65,r66,r67)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i3=*nfl3;i3<=*nfu3;i3++) {
	for(i1=*nfl1;i1<=*nfu1;i1++) {
		istart=ibxz_f[0+2*(i1+(*n1+1)*i3)];
		n2t=ibxz_f[1+2*(i1+(*n1+1)*i3)]-istart+1;
		if(n2t>0) {
		loffset1=3+7*(i1-*nfl1+nseg1*(istart-*nfl2+(i3-*nfl3)*nseg2));
		loffset2=5+7*(i1-*nfl1+nseg1*(istart-*nfl2+(i3-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i2=0;i2<n2t-4;i2+=4) {
			wai4[i2+0]=x_f[loffset1+(i2+0)*nseg3];
			wai6[i2+0]=x_f[loffset2+(i2+0)*nseg3];
			wai4[i2+1]=x_f[loffset1+(i2+1)*nseg3];
			wai6[i2+1]=x_f[loffset2+(i2+1)*nseg3];
			wai4[i2+2]=x_f[loffset1+(i2+2)*nseg3];
			wai6[i2+2]=x_f[loffset2+(i2+2)*nseg3];
			wai4[i2+3]=x_f[loffset1+(i2+3)*nseg3];
			wai6[i2+3]=x_f[loffset2+(i2+3)*nseg3];
		}
		mm=n2t%4;
		if(mm==0) {
			n2p=n2t;
			wai4[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai6[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai4[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai6[n2p-3]=x_f[loffset2+(n2p-3)*nseg3];
			wai4[n2p-2]=x_f[loffset1+(n2p-2)*nseg3];
			wai6[n2p-2]=x_f[loffset2+(n2p-2)*nseg3];
			wai4[n2p-1]=x_f[loffset1+(n2p-1)*nseg3];
			wai6[n2p-1]=x_f[loffset2+(n2p-1)*nseg3];
		}
		else if(mm==1) {
			n2p=n2t+3;
			wai4[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai6[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai4[n2p-3]=0.0;
			wai6[n2p-3]=0.0;
			wai4[n2p-2]=0.0;
			wai6[n2p-2]=0.0;
			wai4[n2p-1]=0.0;
			wai6[n2p-1]=0.0;
		}
		else if(mm==2) {
			n2p=n2t+2;
			wai4[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai6[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai4[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai6[n2p-3]=x_f[loffset2+(n2p-3)*nseg3];
			wai4[n2p-2]=0.0;
			wai6[n2p-2]=0.0;
			wai4[n2p-1]=0.0;
			wai6[n2p-1]=0.0;
		}
		else if(mm==3) {
			n2p=n2t+1;
			wai4[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai6[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai4[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai6[n2p-3]=x_f[loffset2+(n2p-3)*nseg3];
			wai4[n2p-2]=x_f[loffset1+(n2p-2)*nseg3];
			wai6[n2p-2]=x_f[loffset2+(n2p-2)*nseg3];
			wai4[n2p-1]=0.0;
			wai6[n2p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n2t<5) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			CONV_L1(l41,am03,am02,am01,ap00,ap01,ap02,ap03);
			CONV_L1C(l61,bm03,bm02,bm01,bp00,bp01,bp02,bp03);
			_mm_store_ps(&wao4[n2p- 4],r0);
			CONV_L1(l41,cm03,cm02,cm01,cp00,cp01,cp02,cp03);
			CONV_L1C(l61,em03,em02,em01,ep00,ep01,ep02,ep03);
			_mm_store_ps(&wao6[n2p- 4],r0);
		}
		else if(n2t<9) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			CONV_L2(l41,l42,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			CONV_L2C(l61,l62,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07);
			_mm_store_ps(&wao4[n2p- 8],r1);
			_mm_store_ps(&wao4[n2p- 4],r0);
			CONV_L2(l41,l42,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07);
			CONV_L2C(l61,l62,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			_mm_store_ps(&wao6[n2p- 8],r1);
			_mm_store_ps(&wao6[n2p- 4],r0);
		}
		else if(n2t<13) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			CONV_L3(l41,l42,l43,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			CONV_L3C(l61,l62,l63,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11);
			_mm_store_ps(&wao4[n2p-12],r2);
			_mm_store_ps(&wao4[n2p- 8],r1);
			_mm_store_ps(&wao4[n2p- 4],r0);
			CONV_L3(l41,l42,l43,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11);
			CONV_L3C(l61,l62,l63,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			_mm_store_ps(&wao6[n2p-12],r2);
			_mm_store_ps(&wao6[n2p- 8],r1);
			_mm_store_ps(&wao6[n2p- 4],r0);
		}
		else if(n2t<17) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			CONV_L4(l41,l42,l43,l44,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L4C(l61,l62,l63,l64,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao4[n2p-16],r3);
			_mm_store_ps(&wao4[n2p-12],r2);
			_mm_store_ps(&wao4[n2p- 8],r1);
			_mm_store_ps(&wao4[n2p- 4],r0);
			CONV_L4(l41,l42,l43,l44,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L4C(l61,l62,l63,l64,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao6[n2p-16],r3);
			_mm_store_ps(&wao6[n2p-12],r2);
			_mm_store_ps(&wao6[n2p- 8],r1);
			_mm_store_ps(&wao6[n2p- 4],r0);
		}
		else if(n2t<21) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			CONV_L5(l41,l42,l43,l44,l45,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L5C(l61,l62,l63,l64,l65,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao4[n2p-20],r4);
			_mm_store_ps(&wao4[n2p-16],r3);
			_mm_store_ps(&wao4[n2p-12],r2);
			_mm_store_ps(&wao4[n2p- 8],r1);
			_mm_store_ps(&wao4[n2p- 4],r0);
			CONV_L5(l41,l42,l43,l44,l45,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L5C(l61,l62,l63,l64,l65,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao6[n2p-20],r4);
			_mm_store_ps(&wao6[n2p-16],r3);
			_mm_store_ps(&wao6[n2p-12],r2);
			_mm_store_ps(&wao6[n2p- 8],r1);
			_mm_store_ps(&wao6[n2p- 4],r0);
		}
		else if(n2t<25) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			l46=_mm_castps_si128(_mm_load_ps(wai4+20));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			l66=_mm_castps_si128(_mm_load_ps(wai6+20));
			CONV_L6(l41,l42,l43,l44,l45,l46,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L6C(l61,l62,l63,l64,l65,l66,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao4[n2p-24],r5);
			_mm_store_ps(&wao4[n2p-20],r4);
			_mm_store_ps(&wao4[n2p-16],r3);
			_mm_store_ps(&wao4[n2p-12],r2);
			_mm_store_ps(&wao4[n2p- 8],r1);
			_mm_store_ps(&wao4[n2p- 4],r0);
			CONV_L6(l41,l42,l43,l44,l45,l46,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L6C(l61,l62,l63,l64,l65,l66,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao6[n2p-24],r5);
			_mm_store_ps(&wao6[n2p-20],r4);
			_mm_store_ps(&wao6[n2p-16],r3);
			_mm_store_ps(&wao6[n2p-12],r2);
			_mm_store_ps(&wao6[n2p- 8],r1);
			_mm_store_ps(&wao6[n2p- 4],r0);
		}
		else if(n2t<29) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			l46=_mm_castps_si128(_mm_load_ps(wai4+20));
			l47=_mm_castps_si128(_mm_load_ps(wai4+24));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			l66=_mm_castps_si128(_mm_load_ps(wai6+20));
			l67=_mm_castps_si128(_mm_load_ps(wai6+24));
			CONV_L7(l41,l42,l43,l44,l45,l46,l47,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L7C(l61,l62,l63,l64,l65,l66,l67,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao4[n2p-28],r6);
			_mm_store_ps(&wao4[n2p-24],r5);
			_mm_store_ps(&wao4[n2p-20],r4);
			_mm_store_ps(&wao4[n2p-16],r3);
			_mm_store_ps(&wao4[n2p-12],r2);
			_mm_store_ps(&wao4[n2p- 8],r1);
			_mm_store_ps(&wao4[n2p- 4],r0);
			CONV_L7(l41,l42,l43,l44,l45,l46,l47,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L7C(l61,l62,l63,l64,l65,l66,l67,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao6[n2p-28],r6);
			_mm_store_ps(&wao6[n2p-24],r5);
			_mm_store_ps(&wao6[n2p-20],r4);
			_mm_store_ps(&wao6[n2p-16],r3);
			_mm_store_ps(&wao6[n2p-12],r2);
			_mm_store_ps(&wao6[n2p- 8],r1);
			_mm_store_ps(&wao6[n2p- 4],r0);
		}
		else {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			CONV_L8_PART1(l41,l42,l43,l44,l45,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r41,r42,r43,r44,r45,r46,r47);
			CONV_L8_PART1C(l61,l62,l63,l64,l65,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,r41,r42,r43,r44,r45,r46,r47);
			CONV_L8_PART1(l41,l42,l43,l44,l45,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,r61,r62,r63,r64,r65,r66,r67);
			CONV_L8_PART1C(l61,l62,l63,l64,l65,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r61,r62,r63,r64,r65,r66,r67);
			for(i2=28;i2<n2p-4;i2+=4) {
				l46=_mm_castps_si128(_mm_load_ps(wai4+i2-8));
				l66=_mm_castps_si128(_mm_load_ps(wai6+i2-8));
				CONV_L8_PART2(l44,l45,l46,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r40,r41,r42,r43,r44,r45,r46,r47);
				CONV_L8_PART2C(l64,l65,l66,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
					bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r40,r41,r42,r43,r44,r45,r46,r47);
				_mm_store_ps(&wao4[i2-28],r47);
				CONV_L8_PART2(l44,l45,l46,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
					cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r60,r61,r62,r63,r64,r65,r66,r67);
				CONV_L8_PART2C(l64,l65,l66,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r60,r61,r62,r63,r64,r65,r66,r67);
				_mm_store_ps(&wao6[i2-28],r67);
				l44=l45;l45=l46;
				UP_L8(r40,r41,r42,r43,r44,r45,r46,r47);
				l64=l65;l65=l66;
				UP_L8(r60,r61,r62,r63,r64,r65,r66,r67);
			}
			l46=_mm_castps_si128(_mm_load_ps(wai4+n2p-12));
			l47=_mm_castps_si128(_mm_load_ps(wai4+n2p- 8));
			l48=_mm_castps_si128(_mm_load_ps(wai4+n2p- 4));
			l66=_mm_castps_si128(_mm_load_ps(wai6+n2p-12));
			l67=_mm_castps_si128(_mm_load_ps(wai6+n2p- 8));
			l68=_mm_castps_si128(_mm_load_ps(wai6+n2p- 4));
			CONV_L8_PART3(l44,l45,l46,l47,l48,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r40,r41,r42,r43,r44,r45,r46,r47);
			CONV_L8_PART3C(l64,l65,l66,l67,l68,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r40,r41,r42,r43,r44,r45,r46,r47);
			_mm_store_ps(&wao4[n2p-32],r47);
			_mm_store_ps(&wao4[n2p-28],r46);
			_mm_store_ps(&wao4[n2p-24],r45);
			_mm_store_ps(&wao4[n2p-20],r44);
			_mm_store_ps(&wao4[n2p-16],r43);
			_mm_store_ps(&wao4[n2p-12],r42);
			_mm_store_ps(&wao4[n2p- 8],r41);
			_mm_store_ps(&wao4[n2p- 4],r40);
			CONV_L8_PART3(l44,l45,l46,l47,l48,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r60,r61,r62,r63,r64,r65,r66,r67);
			CONV_L8_PART3C(l64,l65,l66,l67,l68,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r60,r61,r62,r63,r64,r65,r66,r67);
			_mm_store_ps(&wao6[n2p-32],r67);
			_mm_store_ps(&wao6[n2p-28],r66);
			_mm_store_ps(&wao6[n2p-24],r65);
			_mm_store_ps(&wao6[n2p-20],r64);
			_mm_store_ps(&wao6[n2p-16],r63);
			_mm_store_ps(&wao6[n2p-12],r62);
			_mm_store_ps(&wao6[n2p- 8],r61);
			_mm_store_ps(&wao6[n2p- 4],r60);
		}
		//----------------------------------------------------------------------
		for(i2=0;i2<(n2p-4);i2+=4) {
			y_f[loffset1+(i2+0)*nseg3]=y_f[loffset1+(i2+0)*nseg3]+wao4[i2+0];
			y_f[loffset2+(i2+0)*nseg3]=y_f[loffset2+(i2+0)*nseg3]+wao6[i2+0];
			y_f[loffset1+(i2+1)*nseg3]=y_f[loffset1+(i2+1)*nseg3]+wao4[i2+1];
			y_f[loffset2+(i2+1)*nseg3]=y_f[loffset2+(i2+1)*nseg3]+wao6[i2+1];
			y_f[loffset1+(i2+2)*nseg3]=y_f[loffset1+(i2+2)*nseg3]+wao4[i2+2];
			y_f[loffset2+(i2+2)*nseg3]=y_f[loffset2+(i2+2)*nseg3]+wao6[i2+2];
			y_f[loffset1+(i2+3)*nseg3]=y_f[loffset1+(i2+3)*nseg3]+wao4[i2+3];
			y_f[loffset2+(i2+3)*nseg3]=y_f[loffset2+(i2+3)*nseg3]+wao6[i2+3];
		}
		for(;i2<n2t;i2++) {
			y_f[loffset1+i2*nseg3]=y_f[loffset1+i2*nseg3]+wao4[i2];
			y_f[loffset2+i2*nseg3]=y_f[loffset2+i2*nseg3]+wao6[i2];
		}
		} //end of if for n2t>0
		//----------------------------------------------------------------------
	} //end of loop over i1
	} //end of loop over i3
} // end of function convolut_ib_sse_fy_ssw_sww_
//****************************************************************************************
void convolut_ib_sse_fy_wsw_www_(int *n1,int *n2,int *n3,int *ibxz_f,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n2p,nseg1,nseg2,nseg3,n2t,istart;
	__m128 am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01,ap00,ap01;
	__m128 ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14;
	__m128 bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01;
	__m128 bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14;
	__m128 cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01;
	__m128 cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y;
	__m128i lz;
	__m128i l51,l52,l53,l54,l55,l56,l57,l58;
	__m128i l71,l72,l73,l74,l75,l76,l77,l78;
	__m128 r50,r51,r52,r53,r54,r55,r56,r57;
	__m128 r70,r71,r72,r73,r74,r75,r76,r77;
	int i1,i2,i3,loffset1,loffset2;
	float wai5[1024] __attribute__ ((aligned (16)));
	float wai7[1024] __attribute__ ((aligned (16)));
	float wao5[1024] __attribute__ ((aligned (16)));
	float wao7[1024] __attribute__ ((aligned (16)));
	//----------------------------------
	FILTER_A; FILTER_B; FILTER_C; FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	nseg3=nseg1*7;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along Y )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n2t,loffset1,loffset2,mm,n2p,wai5,wai7,wao5,wao7)\
           private(l51,l52,l53,l54,l55,l56,l57,l58,l71,l72,l73,l74,l75,l76,l77,l78)\
           private(r50,r51,r52,r53,r54,r55,r56,r57,r70,r71,r72,r73,r74,r75,r76,r77)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i3=*nfl3;i3<=*nfu3;i3++) {
	for(i1=*nfl1;i1<=*nfu1;i1++) {
		istart=ibxz_f[0+2*(i1+(*n1+1)*i3)];
		n2t=ibxz_f[1+2*(i1+(*n1+1)*i3)]-istart+1;
		if(n2t>0) {
		loffset1=4+7*(i1-*nfl1+nseg1*(istart-*nfl2+(i3-*nfl3)*nseg2));
		loffset2=6+7*(i1-*nfl1+nseg1*(istart-*nfl2+(i3-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i2=0;i2<n2t-4;i2+=4) {
			wai5[i2+0]=x_f[loffset1+(i2+0)*nseg3];
			wai7[i2+0]=x_f[loffset2+(i2+0)*nseg3];
			wai5[i2+1]=x_f[loffset1+(i2+1)*nseg3];
			wai7[i2+1]=x_f[loffset2+(i2+1)*nseg3];
			wai5[i2+2]=x_f[loffset1+(i2+2)*nseg3];
			wai7[i2+2]=x_f[loffset2+(i2+2)*nseg3];
			wai5[i2+3]=x_f[loffset1+(i2+3)*nseg3];
			wai7[i2+3]=x_f[loffset2+(i2+3)*nseg3];
		}
		mm=n2t%4;
		if(mm==0) {
			n2p=n2t;
			wai5[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai7[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai5[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai7[n2p-3]=x_f[loffset2+(n2p-3)*nseg3];
			wai5[n2p-2]=x_f[loffset1+(n2p-2)*nseg3];
			wai7[n2p-2]=x_f[loffset2+(n2p-2)*nseg3];
			wai5[n2p-1]=x_f[loffset1+(n2p-1)*nseg3];
			wai7[n2p-1]=x_f[loffset2+(n2p-1)*nseg3];
		}
		else if(mm==1) {
			n2p=n2t+3;
			wai5[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai7[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai5[n2p-3]=0.0;
			wai7[n2p-3]=0.0;
			wai5[n2p-2]=0.0;
			wai7[n2p-2]=0.0;
			wai5[n2p-1]=0.0;
			wai7[n2p-1]=0.0;
		}
		else if(mm==2) {
			n2p=n2t+2;
			wai5[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai7[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai5[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai7[n2p-3]=x_f[loffset2+(n2p-3)*nseg3];
			wai5[n2p-2]=0.0;
			wai7[n2p-2]=0.0;
			wai5[n2p-1]=0.0;
			wai7[n2p-1]=0.0;
		}
		else if(mm==3) {
			n2p=n2t+1;
			wai5[n2p-4]=x_f[loffset1+(n2p-4)*nseg3];
			wai7[n2p-4]=x_f[loffset2+(n2p-4)*nseg3];
			wai5[n2p-3]=x_f[loffset1+(n2p-3)*nseg3];
			wai7[n2p-3]=x_f[loffset2+(n2p-3)*nseg3];
			wai5[n2p-2]=x_f[loffset1+(n2p-2)*nseg3];
			wai7[n2p-2]=x_f[loffset2+(n2p-2)*nseg3];
			wai5[n2p-1]=0.0;
			wai7[n2p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n2t<5) {
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			CONV_L1(l51,am03,am02,am01,ap00,ap01,ap02,ap03);
			CONV_L1C(l71,bm03,bm02,bm01,bp00,bp01,bp02,bp03);
			_mm_store_ps(&wao5[n2p- 4],r0);
			CONV_L1(l51,cm03,cm02,cm01,cp00,cp01,cp02,cp03);
			CONV_L1C(l71,em03,em02,em01,ep00,ep01,ep02,ep03);
			_mm_store_ps(&wao7[n2p- 4],r0);
		}
		else if(n2t<9) {
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			CONV_L2(l51,l52,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			CONV_L2C(l71,l72,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07);
			_mm_store_ps(&wao5[n2p- 8],r1);
			_mm_store_ps(&wao5[n2p- 4],r0);
			CONV_L2(l51,l52,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07);
			CONV_L2C(l71,l72,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			_mm_store_ps(&wao7[n2p- 8],r1);
			_mm_store_ps(&wao7[n2p- 4],r0);
		}
		else if(n2t<13) {
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			CONV_L3(l51,l52,l53,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			CONV_L3C(l71,l72,l73,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11);
			_mm_store_ps(&wao5[n2p-12],r2);
			_mm_store_ps(&wao5[n2p- 8],r1);
			_mm_store_ps(&wao5[n2p- 4],r0);
			CONV_L3(l51,l52,l53,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11);
			CONV_L3C(l71,l72,l73,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			_mm_store_ps(&wao7[n2p-12],r2);
			_mm_store_ps(&wao7[n2p- 8],r1);
			_mm_store_ps(&wao7[n2p- 4],r0);
		}
		else if(n2t<17) {
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			CONV_L4(l51,l52,l53,l54,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L4C(l71,l72,l73,l74,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao5[n2p-16],r3);
			_mm_store_ps(&wao5[n2p-12],r2);
			_mm_store_ps(&wao5[n2p- 8],r1);
			_mm_store_ps(&wao5[n2p- 4],r0);
			CONV_L4(l51,l52,l53,l54,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L4C(l71,l72,l73,l74,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao7[n2p-16],r3);
			_mm_store_ps(&wao7[n2p-12],r2);
			_mm_store_ps(&wao7[n2p- 8],r1);
			_mm_store_ps(&wao7[n2p- 4],r0);
		}
		else if(n2t<21) {
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			CONV_L5(l51,l52,l53,l54,l55,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L5C(l71,l72,l73,l74,l75,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao5[n2p-20],r4);
			_mm_store_ps(&wao5[n2p-16],r3);
			_mm_store_ps(&wao5[n2p-12],r2);
			_mm_store_ps(&wao5[n2p- 8],r1);
			_mm_store_ps(&wao5[n2p- 4],r0);
			CONV_L5(l51,l52,l53,l54,l55,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L5C(l71,l72,l73,l74,l75,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao7[n2p-20],r4);
			_mm_store_ps(&wao7[n2p-16],r3);
			_mm_store_ps(&wao7[n2p-12],r2);
			_mm_store_ps(&wao7[n2p- 8],r1);
			_mm_store_ps(&wao7[n2p- 4],r0);
		}
		else if(n2t<25) {
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			l56=_mm_castps_si128(_mm_load_ps(wai5+20));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			l76=_mm_castps_si128(_mm_load_ps(wai7+20));
			CONV_L6(l51,l52,l53,l54,l55,l56,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L6C(l71,l72,l73,l74,l75,l76,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao5[n2p-24],r5);
			_mm_store_ps(&wao5[n2p-20],r4);
			_mm_store_ps(&wao5[n2p-16],r3);
			_mm_store_ps(&wao5[n2p-12],r2);
			_mm_store_ps(&wao5[n2p- 8],r1);
			_mm_store_ps(&wao5[n2p- 4],r0);
			CONV_L6(l51,l52,l53,l54,l55,l56,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L6C(l71,l72,l73,l74,l75,l76,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao7[n2p-24],r5);
			_mm_store_ps(&wao7[n2p-20],r4);
			_mm_store_ps(&wao7[n2p-16],r3);
			_mm_store_ps(&wao7[n2p-12],r2);
			_mm_store_ps(&wao7[n2p- 8],r1);
			_mm_store_ps(&wao7[n2p- 4],r0);
		}
		else if(n2t<29) {
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			l56=_mm_castps_si128(_mm_load_ps(wai5+20));
			l57=_mm_castps_si128(_mm_load_ps(wai5+24));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			l76=_mm_castps_si128(_mm_load_ps(wai7+20));
			l77=_mm_castps_si128(_mm_load_ps(wai7+24));
			CONV_L7(l51,l52,l53,l54,l55,l56,l57,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L7C(l71,l72,l73,l74,l75,l76,l77,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao5[n2p-28],r6);
			_mm_store_ps(&wao5[n2p-24],r5);
			_mm_store_ps(&wao5[n2p-20],r4);
			_mm_store_ps(&wao5[n2p-16],r3);
			_mm_store_ps(&wao5[n2p-12],r2);
			_mm_store_ps(&wao5[n2p- 8],r1);
			_mm_store_ps(&wao5[n2p- 4],r0);
			CONV_L7(l51,l52,l53,l54,l55,l56,l57,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L7C(l71,l72,l73,l74,l75,l76,l77,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao7[n2p-28],r6);
			_mm_store_ps(&wao7[n2p-24],r5);
			_mm_store_ps(&wao7[n2p-20],r4);
			_mm_store_ps(&wao7[n2p-16],r3);
			_mm_store_ps(&wao7[n2p-12],r2);
			_mm_store_ps(&wao7[n2p- 8],r1);
			_mm_store_ps(&wao7[n2p- 4],r0);
		}
		else {
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			CONV_L8_PART1(l51,l52,l53,l54,l55,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r51,r52,r53,r54,r55,r56,r57);
			CONV_L8_PART1C(l71,l72,l73,l74,l75,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,r51,r52,r53,r54,r55,r56,r57);
			CONV_L8_PART1(l51,l52,l53,l54,l55,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,r71,r72,r73,r74,r75,r76,r77);
			CONV_L8_PART1C(l71,l72,l73,l74,l75,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r71,r72,r73,r74,r75,r76,r77);
			for(i2=28;i2<n2p-4;i2+=4) {
				l56=_mm_castps_si128(_mm_load_ps(wai5+i2-8));
				l76=_mm_castps_si128(_mm_load_ps(wai7+i2-8));
				CONV_L8_PART2(l54,l55,l56,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r50,r51,r52,r53,r54,r55,r56,r57);
				CONV_L8_PART2C(l74,l75,l76,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
					bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r50,r51,r52,r53,r54,r55,r56,r57);
				_mm_store_ps(&wao5[i2-28],r57);
				CONV_L8_PART2(l54,l55,l56,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
					cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r70,r71,r72,r73,r74,r75,r76,r77);
				CONV_L8_PART2C(l74,l75,l76,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r70,r71,r72,r73,r74,r75,r76,r77);
				_mm_store_ps(&wao7[i2-28],r77);
				l54=l55;l55=l56;
				UP_L8(r50,r51,r52,r53,r54,r55,r56,r57);
				l74=l75;l75=l76;
				UP_L8(r70,r71,r72,r73,r74,r75,r76,r77);
			}
			l56=_mm_castps_si128(_mm_load_ps(wai5+n2p-12));
			l57=_mm_castps_si128(_mm_load_ps(wai5+n2p- 8));
			l58=_mm_castps_si128(_mm_load_ps(wai5+n2p- 4));
			l76=_mm_castps_si128(_mm_load_ps(wai7+n2p-12));
			l77=_mm_castps_si128(_mm_load_ps(wai7+n2p- 8));
			l78=_mm_castps_si128(_mm_load_ps(wai7+n2p- 4));
			CONV_L8_PART3(l54,l55,l56,l57,l58,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r50,r51,r52,r53,r54,r55,r56,r57);
			CONV_L8_PART3C(l74,l75,l76,l77,l78,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r50,r51,r52,r53,r54,r55,r56,r57);
			_mm_store_ps(&wao5[n2p-32],r57);
			_mm_store_ps(&wao5[n2p-28],r56);
			_mm_store_ps(&wao5[n2p-24],r55);
			_mm_store_ps(&wao5[n2p-20],r54);
			_mm_store_ps(&wao5[n2p-16],r53);
			_mm_store_ps(&wao5[n2p-12],r52);
			_mm_store_ps(&wao5[n2p- 8],r51);
			_mm_store_ps(&wao5[n2p- 4],r50);
			CONV_L8_PART3(l54,l55,l56,l57,l58,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r70,r71,r72,r73,r74,r75,r76,r77);
			CONV_L8_PART3C(l74,l75,l76,l77,l78,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r70,r71,r72,r73,r74,r75,r76,r77);
			_mm_store_ps(&wao7[n2p-32],r77);
			_mm_store_ps(&wao7[n2p-28],r76);
			_mm_store_ps(&wao7[n2p-24],r75);
			_mm_store_ps(&wao7[n2p-20],r74);
			_mm_store_ps(&wao7[n2p-16],r73);
			_mm_store_ps(&wao7[n2p-12],r72);
			_mm_store_ps(&wao7[n2p- 8],r71);
			_mm_store_ps(&wao7[n2p- 4],r70);
		}
		//----------------------------------------------------------------------
		for(i2=0;i2<(n2p-4);i2+=4) {
			y_f[loffset1+(i2+0)*nseg3]=y_f[loffset1+(i2+0)*nseg3]+wao5[i2+0];
			y_f[loffset2+(i2+0)*nseg3]=y_f[loffset2+(i2+0)*nseg3]+wao7[i2+0];
			y_f[loffset1+(i2+1)*nseg3]=y_f[loffset1+(i2+1)*nseg3]+wao5[i2+1];
			y_f[loffset2+(i2+1)*nseg3]=y_f[loffset2+(i2+1)*nseg3]+wao7[i2+1];
			y_f[loffset1+(i2+2)*nseg3]=y_f[loffset1+(i2+2)*nseg3]+wao5[i2+2];
			y_f[loffset2+(i2+2)*nseg3]=y_f[loffset2+(i2+2)*nseg3]+wao7[i2+2];
			y_f[loffset1+(i2+3)*nseg3]=y_f[loffset1+(i2+3)*nseg3]+wao5[i2+3];
			y_f[loffset2+(i2+3)*nseg3]=y_f[loffset2+(i2+3)*nseg3]+wao7[i2+3];
		}
		for(;i2<n2t;i2++) {
			y_f[loffset1+i2*nseg3]=y_f[loffset1+i2*nseg3]+wao5[i2];
			y_f[loffset2+i2*nseg3]=y_f[loffset2+i2*nseg3]+wao7[i2];
		}
		} //end of if for n2t>0
		//----------------------------------------------------------------------
	} //end of loop over i1
	} //end of loop over i3
} // end of function convolut_ib_sse_fy_wsw_www_
//****************************************************************************************
void convolut_ib_sse_fz_ssw_(int *n1,int *n2,int *n3,int *ibxy_f,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n3p,nseg1,nseg2,nseg3,n3t,istart;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y;
	__m128i lz;
	__m128i l41,l42,l43,l44,l45,l46,l47,l48;
	__m128 r40,r41,r42,r43,r44,r45,r46,r47;
	int i1,i2,i3,loffset1;
	float wai4[1024] __attribute__ ((aligned (16)));
	float wao4[1024] __attribute__ ((aligned (16)));
	//----------------------------------
	FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	nseg3=nseg1*nseg2*7;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along Z )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n3t,loffset1,mm,n3p,wai4,wao4)\
           private(l41,l42,l43,l44,l45,l46,l47,l48,r40,r41,r42,r43,r44,r45,r46,r47)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i2=*nfl2;i2<=*nfu2;i2++) {
	for(i1=*nfl1;i1<=*nfu1;i1++) {
		istart=ibxy_f[0+2*(i1+(*n1+1)*i2)];
		n3t=ibxy_f[1+2*(i1+(*n1+1)*i2)]-istart+1;
		if(n3t>0) {
		loffset1=3+7*(i1-*nfl1+nseg1*(i2-*nfl2+(istart-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i3=0;i3<n3t-4;i3+=4) {
			wai4[i3+0]=x_f[loffset1+(i3+0)*nseg3];
			wai4[i3+1]=x_f[loffset1+(i3+1)*nseg3];
			wai4[i3+2]=x_f[loffset1+(i3+2)*nseg3];
			wai4[i3+3]=x_f[loffset1+(i3+3)*nseg3];
		}
		mm=n3t%4;
		if(mm==0) {
			n3p=n3t;
			wai4[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai4[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai4[n3p-2]=x_f[loffset1+(n3p-2)*nseg3];
			wai4[n3p-1]=x_f[loffset1+(n3p-1)*nseg3];
		}
		else if(mm==1) {
			n3p=n3t+3;
			wai4[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai4[n3p-3]=0.0;
			wai4[n3p-2]=0.0;
			wai4[n3p-1]=0.0;
		}
		else if(mm==2) {
			n3p=n3t+2;
			wai4[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai4[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai4[n3p-2]=0.0;
			wai4[n3p-1]=0.0;
		}
		else if(mm==3) {
			n3p=n3t+1;
			wai4[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai4[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai4[n3p-2]=x_f[loffset1+(n3p-2)*nseg3];
			wai4[n3p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n3t<5) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			CONV_L1(l41,em03,em02,em01,ep00,ep01,ep02,ep03);
			_mm_store_ps(&wao4[n3p- 4],r0);
		}
		else if(n3t<9) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			CONV_L2(l41,l42,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			_mm_store_ps(&wao4[n3p- 8],r1);
			_mm_store_ps(&wao4[n3p- 4],r0);
		}
		else if(n3t<13) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			CONV_L3(l41,l42,l43,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			_mm_store_ps(&wao4[n3p-12],r2);
			_mm_store_ps(&wao4[n3p- 8],r1);
			_mm_store_ps(&wao4[n3p- 4],r0);
		}
		else if(n3t<17) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			CONV_L4(l41,l42,l43,l44,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao4[n3p-16],r3);
			_mm_store_ps(&wao4[n3p-12],r2);
			_mm_store_ps(&wao4[n3p- 8],r1);
			_mm_store_ps(&wao4[n3p- 4],r0);
		}
		else if(n3t<21) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			CONV_L5(l41,l42,l43,l44,l45,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao4[n3p-20],r4);
			_mm_store_ps(&wao4[n3p-16],r3);
			_mm_store_ps(&wao4[n3p-12],r2);
			_mm_store_ps(&wao4[n3p- 8],r1);
			_mm_store_ps(&wao4[n3p- 4],r0);
		}
		else if(n3t<25) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			l46=_mm_castps_si128(_mm_load_ps(wai4+20));
			CONV_L6(l41,l42,l43,l44,l45,l46,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao4[n3p-24],r5);
			_mm_store_ps(&wao4[n3p-20],r4);
			_mm_store_ps(&wao4[n3p-16],r3);
			_mm_store_ps(&wao4[n3p-12],r2);
			_mm_store_ps(&wao4[n3p- 8],r1);
			_mm_store_ps(&wao4[n3p- 4],r0);
		}
		else if(n3t<29) {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			l46=_mm_castps_si128(_mm_load_ps(wai4+20));
			l47=_mm_castps_si128(_mm_load_ps(wai4+24));
			CONV_L7(l41,l42,l43,l44,l45,l46,l47,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao4[n3p-28],r6);
			_mm_store_ps(&wao4[n3p-24],r5);
			_mm_store_ps(&wao4[n3p-20],r4);
			_mm_store_ps(&wao4[n3p-16],r3);
			_mm_store_ps(&wao4[n3p-12],r2);
			_mm_store_ps(&wao4[n3p- 8],r1);
			_mm_store_ps(&wao4[n3p- 4],r0);
		}
		else {
			l41=_mm_castps_si128(_mm_load_ps(wai4));
			l42=_mm_castps_si128(_mm_load_ps(wai4+4));
			l43=_mm_castps_si128(_mm_load_ps(wai4+8));
			l44=_mm_castps_si128(_mm_load_ps(wai4+12));
			l45=_mm_castps_si128(_mm_load_ps(wai4+16));
			CONV_L8_PART1(l41,l42,l43,l44,l45,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r41,r42,r43,r44,r45,r46,r47);
			for(i3=28;i3<n3p-4;i3+=4) {
				l46=_mm_castps_si128(_mm_load_ps(wai4+i3-8));
				CONV_L8_PART2(l44,l45,l46,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r40,r41,r42,r43,r44,r45,r46,r47);
				_mm_store_ps(&wao4[i3-28],r47);
				l44=l45;l45=l46;
				UP_L8(r40,r41,r42,r43,r44,r45,r46,r47);
			}
			l46=_mm_castps_si128(_mm_load_ps(wai4+n3p-12));
			l47=_mm_castps_si128(_mm_load_ps(wai4+n3p- 8));
			l48=_mm_castps_si128(_mm_load_ps(wai4+n3p- 4));
			CONV_L8_PART3(l44,l45,l46,l47,l48,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r40,r41,r42,r43,r44,r45,r46,r47);
			_mm_store_ps(&wao4[n3p-32],r47);
			_mm_store_ps(&wao4[n3p-28],r46);
			_mm_store_ps(&wao4[n3p-24],r45);
			_mm_store_ps(&wao4[n3p-20],r44);
			_mm_store_ps(&wao4[n3p-16],r43);
			_mm_store_ps(&wao4[n3p-12],r42);
			_mm_store_ps(&wao4[n3p- 8],r41);
			_mm_store_ps(&wao4[n3p- 4],r40);
		}
		//----------------------------------------------------------------------
		for(i3=0;i3<(n3p-4);i3+=4) {
			y_f[loffset1+(i3+0)*nseg3]=y_f[loffset1+(i3+0)*nseg3]+wao4[i3+0];
			y_f[loffset1+(i3+1)*nseg3]=y_f[loffset1+(i3+1)*nseg3]+wao4[i3+1];
			y_f[loffset1+(i3+2)*nseg3]=y_f[loffset1+(i3+2)*nseg3]+wao4[i3+2];
			y_f[loffset1+(i3+3)*nseg3]=y_f[loffset1+(i3+3)*nseg3]+wao4[i3+3];
		}
		for(;i3<n3t;i3++) {
			y_f[loffset1+i3*nseg3]=y_f[loffset1+i3*nseg3]+wao4[i3];
		}
		} //end of if for n3t>0
		//----------------------------------------------------------------------
	} //end of loop over i1
	} //end of loop over i2
} // end of function convolut_ib_sse_fz_ssw_
//****************************************************************************************
void convolut_ib_sse_fz_wss_wsw_(int *n1,int *n2,int *n3,int *ibxy_f,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n3p,nseg1,nseg2,nseg3,n3t,istart;
	__m128 am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01,ap00,ap01;
	__m128 ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14;
	__m128 bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01;
	__m128 bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14;
	__m128 cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01;
	__m128 cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y;
	__m128i lz;
	__m128i l11,l12,l13,l14,l15,l16,l17,l18;
	__m128i l51,l52,l53,l54,l55,l56,l57,l58;
	__m128 r10,r11,r12,r13,r14,r15,r16,r17;
	__m128 r50,r51,r52,r53,r54,r55,r56,r57;
	int i1,i2,i3,loffset1,loffset2;
	float wai1[1024] __attribute__ ((aligned (16)));
	float wai5[1024] __attribute__ ((aligned (16)));
	float wao1[1024] __attribute__ ((aligned (16)));
	float wao5[1024] __attribute__ ((aligned (16)));
	//----------------------------------
	FILTER_A; FILTER_B; FILTER_C; FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	nseg3=nseg1*nseg2*7;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along Z )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n3t,loffset1,loffset2,mm,n3p,wai1,wai5,wao1,wao5)\
           private(l11,l12,l13,l14,l15,l16,l17,l18,l51,l52,l53,l54,l55,l56,l57,l58)\
           private(r10,r11,r12,r13,r14,r15,r16,r17,r50,r51,r52,r53,r54,r55,r56,r57)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i2=*nfl2;i2<=*nfu2;i2++) {
	for(i1=*nfl1;i1<=*nfu1;i1++) {
		istart=ibxy_f[0+2*(i1+(*n1+1)*i2)];
		n3t=ibxy_f[1+2*(i1+(*n1+1)*i2)]-istart+1;
		if(n3t>0) {
		loffset1=0+7*(i1-*nfl1+nseg1*(i2-*nfl2+(istart-*nfl3)*nseg2));
		loffset2=4+7*(i1-*nfl1+nseg1*(i2-*nfl2+(istart-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i3=0;i3<n3t-4;i3+=4) {
			wai1[i3+0]=x_f[loffset1+(i3+0)*nseg3];
			wai5[i3+0]=x_f[loffset2+(i3+0)*nseg3];
			wai1[i3+1]=x_f[loffset1+(i3+1)*nseg3];
			wai5[i3+1]=x_f[loffset2+(i3+1)*nseg3];
			wai1[i3+2]=x_f[loffset1+(i3+2)*nseg3];
			wai5[i3+2]=x_f[loffset2+(i3+2)*nseg3];
			wai1[i3+3]=x_f[loffset1+(i3+3)*nseg3];
			wai5[i3+3]=x_f[loffset2+(i3+3)*nseg3];
		}
		mm=n3t%4;
		if(mm==0) {
			n3p=n3t;
			wai1[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai5[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai1[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai5[n3p-3]=x_f[loffset2+(n3p-3)*nseg3];
			wai1[n3p-2]=x_f[loffset1+(n3p-2)*nseg3];
			wai5[n3p-2]=x_f[loffset2+(n3p-2)*nseg3];
			wai1[n3p-1]=x_f[loffset1+(n3p-1)*nseg3];
			wai5[n3p-1]=x_f[loffset2+(n3p-1)*nseg3];
		}
		else if(mm==1) {
			n3p=n3t+3;
			wai1[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai5[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai1[n3p-3]=0.0;
			wai5[n3p-3]=0.0;
			wai1[n3p-2]=0.0;
			wai5[n3p-2]=0.0;
			wai1[n3p-1]=0.0;
			wai5[n3p-1]=0.0;
		}
		else if(mm==2) {
			n3p=n3t+2;
			wai1[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai5[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai1[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai5[n3p-3]=x_f[loffset2+(n3p-3)*nseg3];
			wai1[n3p-2]=0.0;
			wai5[n3p-2]=0.0;
			wai1[n3p-1]=0.0;
			wai5[n3p-1]=0.0;
		}
		else if(mm==3) {
			n3p=n3t+1;
			wai1[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai5[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai1[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai5[n3p-3]=x_f[loffset2+(n3p-3)*nseg3];
			wai1[n3p-2]=x_f[loffset1+(n3p-2)*nseg3];
			wai5[n3p-2]=x_f[loffset2+(n3p-2)*nseg3];
			wai1[n3p-1]=0.0;
			wai5[n3p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n3t<5) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			CONV_L1(l11,am03,am02,am01,ap00,ap01,ap02,ap03);
			CONV_L1C(l51,bm03,bm02,bm01,bp00,bp01,bp02,bp03);
			_mm_store_ps(&wao1[n3p- 4],r0);
			CONV_L1(l11,cm03,cm02,cm01,cp00,cp01,cp02,cp03);
			CONV_L1C(l51,em03,em02,em01,ep00,ep01,ep02,ep03);
			_mm_store_ps(&wao5[n3p- 4],r0);
		}
		else if(n3t<9) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			CONV_L2(l11,l12,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			CONV_L2C(l51,l52,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07);
			_mm_store_ps(&wao1[n3p- 8],r1);
			_mm_store_ps(&wao1[n3p- 4],r0);
			CONV_L2(l11,l12,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07);
			CONV_L2C(l51,l52,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			_mm_store_ps(&wao5[n3p- 8],r1);
			_mm_store_ps(&wao5[n3p- 4],r0);
		}
		else if(n3t<13) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			CONV_L3(l11,l12,l13,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			CONV_L3C(l51,l52,l53,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11);
			_mm_store_ps(&wao1[n3p-12],r2);
			_mm_store_ps(&wao1[n3p- 8],r1);
			_mm_store_ps(&wao1[n3p- 4],r0);
			CONV_L3(l11,l12,l13,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11);
			CONV_L3C(l51,l52,l53,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			_mm_store_ps(&wao5[n3p-12],r2);
			_mm_store_ps(&wao5[n3p- 8],r1);
			_mm_store_ps(&wao5[n3p- 4],r0);
		}
		else if(n3t<17) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			CONV_L4(l11,l12,l13,l14,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L4C(l51,l52,l53,l54,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao1[n3p-16],r3);
			_mm_store_ps(&wao1[n3p-12],r2);
			_mm_store_ps(&wao1[n3p- 8],r1);
			_mm_store_ps(&wao1[n3p- 4],r0);
			CONV_L4(l11,l12,l13,l14,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L4C(l51,l52,l53,l54,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao5[n3p-16],r3);
			_mm_store_ps(&wao5[n3p-12],r2);
			_mm_store_ps(&wao5[n3p- 8],r1);
			_mm_store_ps(&wao5[n3p- 4],r0);
		}
		else if(n3t<21) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			CONV_L5(l11,l12,l13,l14,l15,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L5C(l51,l52,l53,l54,l55,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao1[n3p-20],r4);
			_mm_store_ps(&wao1[n3p-16],r3);
			_mm_store_ps(&wao1[n3p-12],r2);
			_mm_store_ps(&wao1[n3p- 8],r1);
			_mm_store_ps(&wao1[n3p- 4],r0);
			CONV_L5(l11,l12,l13,l14,l15,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L5C(l51,l52,l53,l54,l55,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao5[n3p-20],r4);
			_mm_store_ps(&wao5[n3p-16],r3);
			_mm_store_ps(&wao5[n3p-12],r2);
			_mm_store_ps(&wao5[n3p- 8],r1);
			_mm_store_ps(&wao5[n3p- 4],r0);
		}
		else if(n3t<25) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));
			l16=_mm_castps_si128(_mm_load_ps(wai1+20));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			l56=_mm_castps_si128(_mm_load_ps(wai5+20));
			CONV_L6(l11,l12,l13,l14,l15,l16,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L6C(l51,l52,l53,l54,l55,l56,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao1[n3p-24],r5);
			_mm_store_ps(&wao1[n3p-20],r4);
			_mm_store_ps(&wao1[n3p-16],r3);
			_mm_store_ps(&wao1[n3p-12],r2);
			_mm_store_ps(&wao1[n3p- 8],r1);
			_mm_store_ps(&wao1[n3p- 4],r0);
			CONV_L6(l11,l12,l13,l14,l15,l16,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L6C(l51,l52,l53,l54,l55,l56,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao5[n3p-24],r5);
			_mm_store_ps(&wao5[n3p-20],r4);
			_mm_store_ps(&wao5[n3p-16],r3);
			_mm_store_ps(&wao5[n3p-12],r2);
			_mm_store_ps(&wao5[n3p- 8],r1);
			_mm_store_ps(&wao5[n3p- 4],r0);
		}
		else if(n3t<29) {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));
			l16=_mm_castps_si128(_mm_load_ps(wai1+20));
			l17=_mm_castps_si128(_mm_load_ps(wai1+24));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			l56=_mm_castps_si128(_mm_load_ps(wai5+20));
			l57=_mm_castps_si128(_mm_load_ps(wai5+24));
			CONV_L7(l11,l12,l13,l14,l15,l16,l17,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L7C(l51,l52,l53,l54,l55,l56,l57,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao1[n3p-28],r6);
			_mm_store_ps(&wao1[n3p-24],r5);
			_mm_store_ps(&wao1[n3p-20],r4);
			_mm_store_ps(&wao1[n3p-16],r3);
			_mm_store_ps(&wao1[n3p-12],r2);
			_mm_store_ps(&wao1[n3p- 8],r1);
			_mm_store_ps(&wao1[n3p- 4],r0);
			CONV_L7(l11,l12,l13,l14,l15,l16,l17,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L7C(l51,l52,l53,l54,l55,l56,l57,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao5[n3p-28],r6);
			_mm_store_ps(&wao5[n3p-24],r5);
			_mm_store_ps(&wao5[n3p-20],r4);
			_mm_store_ps(&wao5[n3p-16],r3);
			_mm_store_ps(&wao5[n3p-12],r2);
			_mm_store_ps(&wao5[n3p- 8],r1);
			_mm_store_ps(&wao5[n3p- 4],r0);
		}
		else {
			l11=_mm_castps_si128(_mm_load_ps(wai1));
			l12=_mm_castps_si128(_mm_load_ps(wai1+4));
			l13=_mm_castps_si128(_mm_load_ps(wai1+8));
			l14=_mm_castps_si128(_mm_load_ps(wai1+12));
			l15=_mm_castps_si128(_mm_load_ps(wai1+16));
			l51=_mm_castps_si128(_mm_load_ps(wai5));
			l52=_mm_castps_si128(_mm_load_ps(wai5+4));
			l53=_mm_castps_si128(_mm_load_ps(wai5+8));
			l54=_mm_castps_si128(_mm_load_ps(wai5+12));
			l55=_mm_castps_si128(_mm_load_ps(wai5+16));
			CONV_L8_PART1(l11,l12,l13,l14,l15,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r11,r12,r13,r14,r15,r16,r17);
			CONV_L8_PART1C(l51,l52,l53,l54,l55,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,r11,r12,r13,r14,r15,r16,r17);
			CONV_L8_PART1(l11,l12,l13,l14,l15,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,r51,r52,r53,r54,r55,r56,r57);
			CONV_L8_PART1C(l51,l52,l53,l54,l55,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r51,r52,r53,r54,r55,r56,r57);
			for(i3=28;i3<n3p-4;i3+=4) {
				l16=_mm_castps_si128(_mm_load_ps(wai1+i3-8));
				l56=_mm_castps_si128(_mm_load_ps(wai5+i3-8));
				CONV_L8_PART2(l14,l15,l16,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r10,r11,r12,r13,r14,r15,r16,r17);
				CONV_L8_PART2C(l54,l55,l56,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
					bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r10,r11,r12,r13,r14,r15,r16,r17);
				_mm_store_ps(&wao1[i3-28],r17);
				CONV_L8_PART2(l14,l15,l16,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
					cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r50,r51,r52,r53,r54,r55,r56,r57);
				CONV_L8_PART2C(l54,l55,l56,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r50,r51,r52,r53,r54,r55,r56,r57);
				_mm_store_ps(&wao5[i3-28],r57);
				l14=l15;l15=l16;
				UP_L8(r10,r11,r12,r13,r14,r15,r16,r17);
				l54=l55;l55=l56;
				UP_L8(r50,r51,r52,r53,r54,r55,r56,r57);
			}
			l16=_mm_castps_si128(_mm_load_ps(wai1+n3p-12));
			l17=_mm_castps_si128(_mm_load_ps(wai1+n3p- 8));
			l18=_mm_castps_si128(_mm_load_ps(wai1+n3p- 4));
			l56=_mm_castps_si128(_mm_load_ps(wai5+n3p-12));
			l57=_mm_castps_si128(_mm_load_ps(wai5+n3p- 8));
			l58=_mm_castps_si128(_mm_load_ps(wai5+n3p- 4));
			CONV_L8_PART3(l14,l15,l16,l17,l18,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r10,r11,r12,r13,r14,r15,r16,r17);
			CONV_L8_PART3C(l54,l55,l56,l57,l58,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r10,r11,r12,r13,r14,r15,r16,r17);
			_mm_store_ps(&wao1[n3p-32],r17);
			_mm_store_ps(&wao1[n3p-28],r16);
			_mm_store_ps(&wao1[n3p-24],r15);
			_mm_store_ps(&wao1[n3p-20],r14);
			_mm_store_ps(&wao1[n3p-16],r13);
			_mm_store_ps(&wao1[n3p-12],r12);
			_mm_store_ps(&wao1[n3p- 8],r11);
			_mm_store_ps(&wao1[n3p- 4],r10);
			CONV_L8_PART3(l14,l15,l16,l17,l18,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r50,r51,r52,r53,r54,r55,r56,r57);
			CONV_L8_PART3C(l54,l55,l56,l57,l58,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r50,r51,r52,r53,r54,r55,r56,r57);
			_mm_store_ps(&wao5[n3p-32],r57);
			_mm_store_ps(&wao5[n3p-28],r56);
			_mm_store_ps(&wao5[n3p-24],r55);
			_mm_store_ps(&wao5[n3p-20],r54);
			_mm_store_ps(&wao5[n3p-16],r53);
			_mm_store_ps(&wao5[n3p-12],r52);
			_mm_store_ps(&wao5[n3p- 8],r51);
			_mm_store_ps(&wao5[n3p- 4],r50);
		}
		//----------------------------------------------------------------------
		for(i3=0;i3<(n3p-4);i3+=4) {
			y_f[loffset1+(i3+0)*nseg3]=y_f[loffset1+(i3+0)*nseg3]+wao1[i3+0];
			y_f[loffset2+(i3+0)*nseg3]=y_f[loffset2+(i3+0)*nseg3]+wao5[i3+0];
			y_f[loffset1+(i3+1)*nseg3]=y_f[loffset1+(i3+1)*nseg3]+wao1[i3+1];
			y_f[loffset2+(i3+1)*nseg3]=y_f[loffset2+(i3+1)*nseg3]+wao5[i3+1];
			y_f[loffset1+(i3+2)*nseg3]=y_f[loffset1+(i3+2)*nseg3]+wao1[i3+2];
			y_f[loffset2+(i3+2)*nseg3]=y_f[loffset2+(i3+2)*nseg3]+wao5[i3+2];
			y_f[loffset1+(i3+3)*nseg3]=y_f[loffset1+(i3+3)*nseg3]+wao1[i3+3];
			y_f[loffset2+(i3+3)*nseg3]=y_f[loffset2+(i3+3)*nseg3]+wao5[i3+3];
		}
		for(;i3<n3t;i3++) {
			y_f[loffset1+i3*nseg3]=y_f[loffset1+i3*nseg3]+wao1[i3];
			y_f[loffset2+i3*nseg3]=y_f[loffset2+i3*nseg3]+wao5[i3];
		}
		} //end of if for n3t>0
		//----------------------------------------------------------------------
	} //end of loop over i1
	} //end of loop over i2
} // end of function convolut_ib_sse_fz_wss_wsw_
//****************************************************************************************
void convolut_ib_sse_fz_sws_sww_(int *n1,int *n2,int *n3,int *ibxy_f,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n3p,nseg1,nseg2,nseg3,n3t,istart;
	__m128 am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01,ap00,ap01;
	__m128 ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14;
	__m128 bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01;
	__m128 bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14;
	__m128 cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01;
	__m128 cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y;
	__m128i lz;
	__m128i l21,l22,l23,l24,l25,l26,l27,l28;
	__m128i l61,l62,l63,l64,l65,l66,l67,l68;
	__m128 r20,r21,r22,r23,r24,r25,r26,r27;
	__m128 r60,r61,r62,r63,r64,r65,r66,r67;
	int i1,i2,i3,loffset1,loffset2;
	float wai2[1024] __attribute__ ((aligned (16)));
	float wai6[1024] __attribute__ ((aligned (16)));
	float wao2[1024] __attribute__ ((aligned (16)));
	float wao6[1024] __attribute__ ((aligned (16)));
	//----------------------------------
	FILTER_A; FILTER_B; FILTER_C; FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	nseg3=nseg1*nseg2*7;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along Z )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n3t,loffset1,loffset2,mm,n3p,wai2,wai6,wao2,wao6)\
           private(l21,l22,l23,l24,l25,l26,l27,l28,l61,l62,l63,l64,l65,l66,l67,l68)\
           private(r20,r21,r22,r23,r24,r25,r26,r27,r60,r61,r62,r63,r64,r65,r66,r67)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i2=*nfl2;i2<=*nfu2;i2++) {
	for(i1=*nfl1;i1<=*nfu1;i1++) {
		istart=ibxy_f[0+2*(i1+(*n1+1)*i2)];
		n3t=ibxy_f[1+2*(i1+(*n1+1)*i2)]-istart+1;
		if(n3t>0) {
		loffset1=1+7*(i1-*nfl1+nseg1*(i2-*nfl2+(istart-*nfl3)*nseg2));
		loffset2=5+7*(i1-*nfl1+nseg1*(i2-*nfl2+(istart-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i3=0;i3<n3t-4;i3+=4) {
			wai2[i3+0]=x_f[loffset1+(i3+0)*nseg3];
			wai6[i3+0]=x_f[loffset2+(i3+0)*nseg3];
			wai2[i3+1]=x_f[loffset1+(i3+1)*nseg3];
			wai6[i3+1]=x_f[loffset2+(i3+1)*nseg3];
			wai2[i3+2]=x_f[loffset1+(i3+2)*nseg3];
			wai6[i3+2]=x_f[loffset2+(i3+2)*nseg3];
			wai2[i3+3]=x_f[loffset1+(i3+3)*nseg3];
			wai6[i3+3]=x_f[loffset2+(i3+3)*nseg3];
		}
		mm=n3t%4;
		if(mm==0) {
			n3p=n3t;
			wai2[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai6[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai2[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai6[n3p-3]=x_f[loffset2+(n3p-3)*nseg3];
			wai2[n3p-2]=x_f[loffset1+(n3p-2)*nseg3];
			wai6[n3p-2]=x_f[loffset2+(n3p-2)*nseg3];
			wai2[n3p-1]=x_f[loffset1+(n3p-1)*nseg3];
			wai6[n3p-1]=x_f[loffset2+(n3p-1)*nseg3];
		}
		else if(mm==1) {
			n3p=n3t+3;
			wai2[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai6[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai2[n3p-3]=0.0;
			wai6[n3p-3]=0.0;
			wai2[n3p-2]=0.0;
			wai6[n3p-2]=0.0;
			wai2[n3p-1]=0.0;
			wai6[n3p-1]=0.0;
		}
		else if(mm==2) {
			n3p=n3t+2;
			wai2[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai6[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai2[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai6[n3p-3]=x_f[loffset2+(n3p-3)*nseg3];
			wai2[n3p-2]=0.0;
			wai6[n3p-2]=0.0;
			wai2[n3p-1]=0.0;
			wai6[n3p-1]=0.0;
		}
		else if(mm==3) {
			n3p=n3t+1;
			wai2[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai6[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai2[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai6[n3p-3]=x_f[loffset2+(n3p-3)*nseg3];
			wai2[n3p-2]=x_f[loffset1+(n3p-2)*nseg3];
			wai6[n3p-2]=x_f[loffset2+(n3p-2)*nseg3];
			wai2[n3p-1]=0.0;
			wai6[n3p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n3t<5) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			CONV_L1(l21,am03,am02,am01,ap00,ap01,ap02,ap03);
			CONV_L1C(l61,bm03,bm02,bm01,bp00,bp01,bp02,bp03);
			_mm_store_ps(&wao2[n3p- 4],r0);
			CONV_L1(l21,cm03,cm02,cm01,cp00,cp01,cp02,cp03);
			CONV_L1C(l61,em03,em02,em01,ep00,ep01,ep02,ep03);
			_mm_store_ps(&wao6[n3p- 4],r0);
		}
		else if(n3t<9) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			CONV_L2(l21,l22,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			CONV_L2C(l61,l62,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07);
			_mm_store_ps(&wao2[n3p- 8],r1);
			_mm_store_ps(&wao2[n3p- 4],r0);
			CONV_L2(l21,l22,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07);
			CONV_L2C(l61,l62,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			_mm_store_ps(&wao6[n3p- 8],r1);
			_mm_store_ps(&wao6[n3p- 4],r0);
		}
		else if(n3t<13) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			CONV_L3(l21,l22,l23,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			CONV_L3C(l61,l62,l63,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11);
			_mm_store_ps(&wao2[n3p-12],r2);
			_mm_store_ps(&wao2[n3p- 8],r1);
			_mm_store_ps(&wao2[n3p- 4],r0);
			CONV_L3(l21,l22,l23,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11);
			CONV_L3C(l61,l62,l63,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			_mm_store_ps(&wao6[n3p-12],r2);
			_mm_store_ps(&wao6[n3p- 8],r1);
			_mm_store_ps(&wao6[n3p- 4],r0);
		}
		else if(n3t<17) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			CONV_L4(l21,l22,l23,l24,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L4C(l61,l62,l63,l64,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao2[n3p-16],r3);
			_mm_store_ps(&wao2[n3p-12],r2);
			_mm_store_ps(&wao2[n3p- 8],r1);
			_mm_store_ps(&wao2[n3p- 4],r0);
			CONV_L4(l21,l22,l23,l24,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L4C(l61,l62,l63,l64,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao6[n3p-16],r3);
			_mm_store_ps(&wao6[n3p-12],r2);
			_mm_store_ps(&wao6[n3p- 8],r1);
			_mm_store_ps(&wao6[n3p- 4],r0);
		}
		else if(n3t<21) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			CONV_L5(l21,l22,l23,l24,l25,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L5C(l61,l62,l63,l64,l65,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao2[n3p-20],r4);
			_mm_store_ps(&wao2[n3p-16],r3);
			_mm_store_ps(&wao2[n3p-12],r2);
			_mm_store_ps(&wao2[n3p- 8],r1);
			_mm_store_ps(&wao2[n3p- 4],r0);
			CONV_L5(l21,l22,l23,l24,l25,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L5C(l61,l62,l63,l64,l65,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao6[n3p-20],r4);
			_mm_store_ps(&wao6[n3p-16],r3);
			_mm_store_ps(&wao6[n3p-12],r2);
			_mm_store_ps(&wao6[n3p- 8],r1);
			_mm_store_ps(&wao6[n3p- 4],r0);
		}
		else if(n3t<25) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));
			l26=_mm_castps_si128(_mm_load_ps(wai2+20));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			l66=_mm_castps_si128(_mm_load_ps(wai6+20));
			CONV_L6(l21,l22,l23,l24,l25,l26,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L6C(l61,l62,l63,l64,l65,l66,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao2[n3p-24],r5);
			_mm_store_ps(&wao2[n3p-20],r4);
			_mm_store_ps(&wao2[n3p-16],r3);
			_mm_store_ps(&wao2[n3p-12],r2);
			_mm_store_ps(&wao2[n3p- 8],r1);
			_mm_store_ps(&wao2[n3p- 4],r0);
			CONV_L6(l21,l22,l23,l24,l25,l26,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L6C(l61,l62,l63,l64,l65,l66,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao6[n3p-24],r5);
			_mm_store_ps(&wao6[n3p-20],r4);
			_mm_store_ps(&wao6[n3p-16],r3);
			_mm_store_ps(&wao6[n3p-12],r2);
			_mm_store_ps(&wao6[n3p- 8],r1);
			_mm_store_ps(&wao6[n3p- 4],r0);
		}
		else if(n3t<29) {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));
			l26=_mm_castps_si128(_mm_load_ps(wai2+20));
			l27=_mm_castps_si128(_mm_load_ps(wai2+24));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			l66=_mm_castps_si128(_mm_load_ps(wai6+20));
			l67=_mm_castps_si128(_mm_load_ps(wai6+24));
			CONV_L7(l21,l22,l23,l24,l25,l26,l27,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L7C(l61,l62,l63,l64,l65,l66,l67,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao2[n3p-28],r6);
			_mm_store_ps(&wao2[n3p-24],r5);
			_mm_store_ps(&wao2[n3p-20],r4);
			_mm_store_ps(&wao2[n3p-16],r3);
			_mm_store_ps(&wao2[n3p-12],r2);
			_mm_store_ps(&wao2[n3p- 8],r1);
			_mm_store_ps(&wao2[n3p- 4],r0);
			CONV_L7(l21,l22,l23,l24,l25,l26,l27,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L7C(l61,l62,l63,l64,l65,l66,l67,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao6[n3p-28],r6);
			_mm_store_ps(&wao6[n3p-24],r5);
			_mm_store_ps(&wao6[n3p-20],r4);
			_mm_store_ps(&wao6[n3p-16],r3);
			_mm_store_ps(&wao6[n3p-12],r2);
			_mm_store_ps(&wao6[n3p- 8],r1);
			_mm_store_ps(&wao6[n3p- 4],r0);
		}
		else {
			l21=_mm_castps_si128(_mm_load_ps(wai2));
			l22=_mm_castps_si128(_mm_load_ps(wai2+4));
			l23=_mm_castps_si128(_mm_load_ps(wai2+8));
			l24=_mm_castps_si128(_mm_load_ps(wai2+12));
			l25=_mm_castps_si128(_mm_load_ps(wai2+16));
			l61=_mm_castps_si128(_mm_load_ps(wai6));
			l62=_mm_castps_si128(_mm_load_ps(wai6+4));
			l63=_mm_castps_si128(_mm_load_ps(wai6+8));
			l64=_mm_castps_si128(_mm_load_ps(wai6+12));
			l65=_mm_castps_si128(_mm_load_ps(wai6+16));
			CONV_L8_PART1(l21,l22,l23,l24,l25,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r21,r22,r23,r24,r25,r26,r27);
			CONV_L8_PART1C(l61,l62,l63,l64,l65,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,r21,r22,r23,r24,r25,r26,r27);
			CONV_L8_PART1(l21,l22,l23,l24,l25,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,r61,r62,r63,r64,r65,r66,r67);
			CONV_L8_PART1C(l61,l62,l63,l64,l65,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r61,r62,r63,r64,r65,r66,r67);
			for(i3=28;i3<n3p-4;i3+=4) {
				l26=_mm_castps_si128(_mm_load_ps(wai2+i3-8));
				l66=_mm_castps_si128(_mm_load_ps(wai6+i3-8));
				CONV_L8_PART2(l24,l25,l26,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r20,r21,r22,r23,r24,r25,r26,r27);
				CONV_L8_PART2C(l64,l65,l66,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
					bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r20,r21,r22,r23,r24,r25,r26,r27);
				_mm_store_ps(&wao2[i3-28],r27);
				CONV_L8_PART2(l24,l25,l26,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
					cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r60,r61,r62,r63,r64,r65,r66,r67);
				CONV_L8_PART2C(l64,l65,l66,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r60,r61,r62,r63,r64,r65,r66,r67);
				_mm_store_ps(&wao6[i3-28],r67);
				l24=l25;l25=l26;
				UP_L8(r20,r21,r22,r23,r24,r25,r26,r27);
				l64=l65;l65=l66;
				UP_L8(r60,r61,r62,r63,r64,r65,r66,r67);
			}
			l26=_mm_castps_si128(_mm_load_ps(wai2+n3p-12));
			l27=_mm_castps_si128(_mm_load_ps(wai2+n3p- 8));
			l28=_mm_castps_si128(_mm_load_ps(wai2+n3p- 4));
			l66=_mm_castps_si128(_mm_load_ps(wai6+n3p-12));
			l67=_mm_castps_si128(_mm_load_ps(wai6+n3p- 8));
			l68=_mm_castps_si128(_mm_load_ps(wai6+n3p- 4));
			CONV_L8_PART3(l24,l25,l26,l27,l28,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r20,r21,r22,r23,r24,r25,r26,r27);
			CONV_L8_PART3C(l64,l65,l66,l67,l68,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r20,r21,r22,r23,r24,r25,r26,r27);
			_mm_store_ps(&wao2[n3p-32],r27);
			_mm_store_ps(&wao2[n3p-28],r26);
			_mm_store_ps(&wao2[n3p-24],r25);
			_mm_store_ps(&wao2[n3p-20],r24);
			_mm_store_ps(&wao2[n3p-16],r23);
			_mm_store_ps(&wao2[n3p-12],r22);
			_mm_store_ps(&wao2[n3p- 8],r21);
			_mm_store_ps(&wao2[n3p- 4],r20);
			CONV_L8_PART3(l24,l25,l26,l27,l28,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r60,r61,r62,r63,r64,r65,r66,r67);
			CONV_L8_PART3C(l64,l65,l66,l67,l68,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r60,r61,r62,r63,r64,r65,r66,r67);
			_mm_store_ps(&wao6[n3p-32],r67);
			_mm_store_ps(&wao6[n3p-28],r66);
			_mm_store_ps(&wao6[n3p-24],r65);
			_mm_store_ps(&wao6[n3p-20],r64);
			_mm_store_ps(&wao6[n3p-16],r63);
			_mm_store_ps(&wao6[n3p-12],r62);
			_mm_store_ps(&wao6[n3p- 8],r61);
			_mm_store_ps(&wao6[n3p- 4],r60);
		}
		//----------------------------------------------------------------------
		for(i3=0;i3<(n3p-4);i3+=4) {
			y_f[loffset1+(i3+0)*nseg3]=y_f[loffset1+(i3+0)*nseg3]+wao2[i3+0];
			y_f[loffset2+(i3+0)*nseg3]=y_f[loffset2+(i3+0)*nseg3]+wao6[i3+0];
			y_f[loffset1+(i3+1)*nseg3]=y_f[loffset1+(i3+1)*nseg3]+wao2[i3+1];
			y_f[loffset2+(i3+1)*nseg3]=y_f[loffset2+(i3+1)*nseg3]+wao6[i3+1];
			y_f[loffset1+(i3+2)*nseg3]=y_f[loffset1+(i3+2)*nseg3]+wao2[i3+2];
			y_f[loffset2+(i3+2)*nseg3]=y_f[loffset2+(i3+2)*nseg3]+wao6[i3+2];
			y_f[loffset1+(i3+3)*nseg3]=y_f[loffset1+(i3+3)*nseg3]+wao2[i3+3];
			y_f[loffset2+(i3+3)*nseg3]=y_f[loffset2+(i3+3)*nseg3]+wao6[i3+3];
		}
		for(;i3<n3t;i3++) {
			y_f[loffset1+i3*nseg3]=y_f[loffset1+i3*nseg3]+wao2[i3];
			y_f[loffset2+i3*nseg3]=y_f[loffset2+i3*nseg3]+wao6[i3];
		}
		} //end of if for n3t>0
		//----------------------------------------------------------------------
	} //end of loop over i1
	} //end of loop over i2
} // end of function convolut_ib_sse_fz_sws_sww_
//****************************************************************************************
void convolut_ib_sse_fz_wws_www_(int *n1,int *n2,int *n3,int *ibxy_f,
	int *nfl1,int *nfu1,int *nfl2,int *nfu2,int *nfl3,int *nfu3,
	float *x_f,float *y_f,float *fa,float *fb,float *fc,float *fe) {
	int mm,n3p,nseg1,nseg2,nseg3,n3t,istart;
	__m128 am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01,ap00,ap01;
	__m128 ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14;
	__m128 bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01;
	__m128 bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14;
	__m128 cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01;
	__m128 cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14;
	__m128 em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01,ep00,ep01;
	__m128 ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14;
	__m128 r0,r1,r2,r3,r4,r5,r6,r7,x,y;
	__m128i lz;
	__m128i l31,l32,l33,l34,l35,l36,l37,l38;
	__m128i l71,l72,l73,l74,l75,l76,l77,l78;
	__m128 r30,r31,r32,r33,r34,r35,r36,r37;
	__m128 r70,r71,r72,r73,r74,r75,r76,r77;
	int i1,i2,i3,loffset1,loffset2;
	float wai3[1024] __attribute__ ((aligned (16)));
	float wai7[1024] __attribute__ ((aligned (16)));
	float wao3[1024] __attribute__ ((aligned (16)));
	float wao7[1024] __attribute__ ((aligned (16)));
	//----------------------------------
	FILTER_A; FILTER_B; FILTER_C; FILTER_E;
	//----------------------------------
	lz=_mm_castps_si128(_mm_setzero_ps());
	nseg1=(*nfu1-*nfl1+1);
	nseg2=(*nfu2-*nfl2+1);
	nseg3=nseg1*nseg2*7;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//FINE ---> FINE ( along Z )
        #pragma omp parallel for schedule(static,1)\
           default(shared)\
           private(i3,i2,i1,istart,n3t,loffset1,loffset2,mm,n3p,wai3,wai7,wao3,wao7)\
           private(l31,l32,l33,l34,l35,l36,l37,l38,l71,l72,l73,l74,l75,l76,l77,l78)\
           private(r30,r31,r32,r33,r34,r35,r36,r37,r70,r71,r72,r73,r74,r75,r76,r77)\
           private(r0,r1,r2,r3,r4,r5,r6,r7,x,y,lz)
	for(i2=*nfl2;i2<=*nfu2;i2++) {
	for(i1=*nfl1;i1<=*nfu1;i1++) {
		istart=ibxy_f[0+2*(i1+(*n1+1)*i2)];
		n3t=ibxy_f[1+2*(i1+(*n1+1)*i2)]-istart+1;
		if(n3t>0) {
		loffset1=2+7*(i1-*nfl1+nseg1*(i2-*nfl2+(istart-*nfl3)*nseg2));
		loffset2=6+7*(i1-*nfl1+nseg1*(i2-*nfl2+(istart-*nfl3)*nseg2));
		//----------------------------------------------------------------------
		for(i3=0;i3<n3t-4;i3+=4) {
			wai3[i3+0]=x_f[loffset1+(i3+0)*nseg3];
			wai7[i3+0]=x_f[loffset2+(i3+0)*nseg3];
			wai3[i3+1]=x_f[loffset1+(i3+1)*nseg3];
			wai7[i3+1]=x_f[loffset2+(i3+1)*nseg3];
			wai3[i3+2]=x_f[loffset1+(i3+2)*nseg3];
			wai7[i3+2]=x_f[loffset2+(i3+2)*nseg3];
			wai3[i3+3]=x_f[loffset1+(i3+3)*nseg3];
			wai7[i3+3]=x_f[loffset2+(i3+3)*nseg3];
		}
		mm=n3t%4;
		if(mm==0) {
			n3p=n3t;
			wai3[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai7[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai3[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai7[n3p-3]=x_f[loffset2+(n3p-3)*nseg3];
			wai3[n3p-2]=x_f[loffset1+(n3p-2)*nseg3];
			wai7[n3p-2]=x_f[loffset2+(n3p-2)*nseg3];
			wai3[n3p-1]=x_f[loffset1+(n3p-1)*nseg3];
			wai7[n3p-1]=x_f[loffset2+(n3p-1)*nseg3];
		}
		else if(mm==1) {
			n3p=n3t+3;
			wai3[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai7[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai3[n3p-3]=0.0;
			wai7[n3p-3]=0.0;
			wai3[n3p-2]=0.0;
			wai7[n3p-2]=0.0;
			wai3[n3p-1]=0.0;
			wai7[n3p-1]=0.0;
		}
		else if(mm==2) {
			n3p=n3t+2;
			wai3[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai7[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai3[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai7[n3p-3]=x_f[loffset2+(n3p-3)*nseg3];
			wai3[n3p-2]=0.0;
			wai7[n3p-2]=0.0;
			wai3[n3p-1]=0.0;
			wai7[n3p-1]=0.0;
		}
		else if(mm==3) {
			n3p=n3t+1;
			wai3[n3p-4]=x_f[loffset1+(n3p-4)*nseg3];
			wai7[n3p-4]=x_f[loffset2+(n3p-4)*nseg3];
			wai3[n3p-3]=x_f[loffset1+(n3p-3)*nseg3];
			wai7[n3p-3]=x_f[loffset2+(n3p-3)*nseg3];
			wai3[n3p-2]=x_f[loffset1+(n3p-2)*nseg3];
			wai7[n3p-2]=x_f[loffset2+(n3p-2)*nseg3];
			wai3[n3p-1]=0.0;
			wai7[n3p-1]=0.0;
		}
		//----------------------------------------------------------------------
		if(n3t<5) {
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			CONV_L1(l31,am03,am02,am01,ap00,ap01,ap02,ap03);
			CONV_L1C(l71,bm03,bm02,bm01,bp00,bp01,bp02,bp03);
			_mm_store_ps(&wao3[n3p- 4],r0);
			CONV_L1(l31,cm03,cm02,cm01,cp00,cp01,cp02,cp03);
			CONV_L1C(l71,em03,em02,em01,ep00,ep01,ep02,ep03);
			_mm_store_ps(&wao7[n3p- 4],r0);
		}
		else if(n3t<9) {
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			CONV_L2(l31,l32,am07,am06,am05,am04,am03,am02,am01,ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07);
			CONV_L2C(l71,l72,bm07,bm06,bm05,bm04,bm03,bm02,bm01,bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07);
			_mm_store_ps(&wao3[n3p- 8],r1);
			_mm_store_ps(&wao3[n3p- 4],r0);
			CONV_L2(l31,l32,cm07,cm06,cm05,cm04,cm03,cm02,cm01,cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07);
			CONV_L2C(l71,l72,em07,em06,em05,em04,em03,em02,em01,ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07);
			_mm_store_ps(&wao7[n3p- 8],r1);
			_mm_store_ps(&wao7[n3p- 4],r0);
		}
		else if(n3t<13) {
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			CONV_L3(l31,l32,l33,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11);
			CONV_L3C(l71,l72,l73,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11);
			_mm_store_ps(&wao3[n3p-12],r2);
			_mm_store_ps(&wao3[n3p- 8],r1);
			_mm_store_ps(&wao3[n3p- 4],r0);
			CONV_L3(l31,l32,l33,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11);
			CONV_L3C(l71,l72,l73,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11);
			_mm_store_ps(&wao7[n3p-12],r2);
			_mm_store_ps(&wao7[n3p- 8],r1);
			_mm_store_ps(&wao7[n3p- 4],r0);
		}
		else if(n3t<17) {
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			CONV_L4(l31,l32,l33,l34,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L4C(l71,l72,l73,l74,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao3[n3p-16],r3);
			_mm_store_ps(&wao3[n3p-12],r2);
			_mm_store_ps(&wao3[n3p- 8],r1);
			_mm_store_ps(&wao3[n3p- 4],r0);
			CONV_L4(l31,l32,l33,l34,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L4C(l71,l72,l73,l74,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao7[n3p-16],r3);
			_mm_store_ps(&wao7[n3p-12],r2);
			_mm_store_ps(&wao7[n3p- 8],r1);
			_mm_store_ps(&wao7[n3p- 4],r0);
		}
		else if(n3t<21) {
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			CONV_L5(l31,l32,l33,l34,l35,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L5C(l71,l72,l73,l74,l75,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao3[n3p-20],r4);
			_mm_store_ps(&wao3[n3p-16],r3);
			_mm_store_ps(&wao3[n3p-12],r2);
			_mm_store_ps(&wao3[n3p- 8],r1);
			_mm_store_ps(&wao3[n3p- 4],r0);
			CONV_L5(l31,l32,l33,l34,l35,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L5C(l71,l72,l73,l74,l75,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao7[n3p-20],r4);
			_mm_store_ps(&wao7[n3p-16],r3);
			_mm_store_ps(&wao7[n3p-12],r2);
			_mm_store_ps(&wao7[n3p- 8],r1);
			_mm_store_ps(&wao7[n3p- 4],r0);
		}
		else if(n3t<25) {
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));
			l36=_mm_castps_si128(_mm_load_ps(wai3+20));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			l76=_mm_castps_si128(_mm_load_ps(wai7+20));
			CONV_L6(l31,l32,l33,l34,l35,l36,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L6C(l71,l72,l73,l74,l75,l76,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao3[n3p-24],r5);
			_mm_store_ps(&wao3[n3p-20],r4);
			_mm_store_ps(&wao3[n3p-16],r3);
			_mm_store_ps(&wao3[n3p-12],r2);
			_mm_store_ps(&wao3[n3p- 8],r1);
			_mm_store_ps(&wao3[n3p- 4],r0);
			CONV_L6(l31,l32,l33,l34,l35,l36,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L6C(l71,l72,l73,l74,l75,l76,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao7[n3p-24],r5);
			_mm_store_ps(&wao7[n3p-20],r4);
			_mm_store_ps(&wao7[n3p-16],r3);
			_mm_store_ps(&wao7[n3p-12],r2);
			_mm_store_ps(&wao7[n3p- 8],r1);
			_mm_store_ps(&wao7[n3p- 4],r0);
		}
		else if(n3t<29) {
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));
			l36=_mm_castps_si128(_mm_load_ps(wai3+20));
			l37=_mm_castps_si128(_mm_load_ps(wai3+24));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			l76=_mm_castps_si128(_mm_load_ps(wai7+20));
			l77=_mm_castps_si128(_mm_load_ps(wai7+24));
			CONV_L7(l31,l32,l33,l34,l35,l36,l37,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14);
			CONV_L7C(l71,l72,l73,l74,l75,l76,l77,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14);
			_mm_store_ps(&wao3[n3p-28],r6);
			_mm_store_ps(&wao3[n3p-24],r5);
			_mm_store_ps(&wao3[n3p-20],r4);
			_mm_store_ps(&wao3[n3p-16],r3);
			_mm_store_ps(&wao3[n3p-12],r2);
			_mm_store_ps(&wao3[n3p- 8],r1);
			_mm_store_ps(&wao3[n3p- 4],r0);
			CONV_L7(l31,l32,l33,l34,l35,l36,l37,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14);
			CONV_L7C(l71,l72,l73,l74,l75,l76,l77,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14);
			_mm_store_ps(&wao7[n3p-28],r6);
			_mm_store_ps(&wao7[n3p-24],r5);
			_mm_store_ps(&wao7[n3p-20],r4);
			_mm_store_ps(&wao7[n3p-16],r3);
			_mm_store_ps(&wao7[n3p-12],r2);
			_mm_store_ps(&wao7[n3p- 8],r1);
			_mm_store_ps(&wao7[n3p- 4],r0);
		}
		else {
			l31=_mm_castps_si128(_mm_load_ps(wai3));
			l32=_mm_castps_si128(_mm_load_ps(wai3+4));
			l33=_mm_castps_si128(_mm_load_ps(wai3+8));
			l34=_mm_castps_si128(_mm_load_ps(wai3+12));
			l35=_mm_castps_si128(_mm_load_ps(wai3+16));
			l71=_mm_castps_si128(_mm_load_ps(wai7));
			l72=_mm_castps_si128(_mm_load_ps(wai7+4));
			l73=_mm_castps_si128(_mm_load_ps(wai7+8));
			l74=_mm_castps_si128(_mm_load_ps(wai7+12));
			l75=_mm_castps_si128(_mm_load_ps(wai7+16));
			CONV_L8_PART1(l31,l32,l33,l34,l35,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,r31,r32,r33,r34,r35,r36,r37);
			CONV_L8_PART1C(l71,l72,l73,l74,l75,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,r31,r32,r33,r34,r35,r36,r37);
			CONV_L8_PART1(l31,l32,l33,l34,l35,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,r71,r72,r73,r74,r75,r76,r77);
			CONV_L8_PART1C(l71,l72,l73,l74,l75,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,r71,r72,r73,r74,r75,r76,r77);
			for(i3=28;i3<n3p-4;i3+=4) {
				l36=_mm_castps_si128(_mm_load_ps(wai3+i3-8));
				l76=_mm_castps_si128(_mm_load_ps(wai7+i3-8));
				CONV_L8_PART2(l34,l35,l36,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
					ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r30,r31,r32,r33,r34,r35,r36,r37);
				CONV_L8_PART2C(l74,l75,l76,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
					bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r30,r31,r32,r33,r34,r35,r36,r37);
				_mm_store_ps(&wao3[i3-28],r37);
				CONV_L8_PART2(l34,l35,l36,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
					cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r70,r71,r72,r73,r74,r75,r76,r77);
				CONV_L8_PART2C(l74,l75,l76,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
					ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r70,r71,r72,r73,r74,r75,r76,r77);
				_mm_store_ps(&wao7[i3-28],r77);
				l34=l35;l35=l36;
				UP_L8(r30,r31,r32,r33,r34,r35,r36,r37);
				l74=l75;l75=l76;
				UP_L8(r70,r71,r72,r73,r74,r75,r76,r77);
			}
			l36=_mm_castps_si128(_mm_load_ps(wai3+n3p-12));
			l37=_mm_castps_si128(_mm_load_ps(wai3+n3p- 8));
			l38=_mm_castps_si128(_mm_load_ps(wai3+n3p- 4));
			l76=_mm_castps_si128(_mm_load_ps(wai7+n3p-12));
			l77=_mm_castps_si128(_mm_load_ps(wai7+n3p- 8));
			l78=_mm_castps_si128(_mm_load_ps(wai7+n3p- 4));
			CONV_L8_PART3(l34,l35,l36,l37,l38,am14,am13,am12,am11,am10,am09,am08,am07,am06,am05,am04,am03,am02,am01, \
				ap00,ap01,ap02,ap03,ap04,ap05,ap06,ap07,ap08,ap09,ap10,ap11,ap12,ap13,ap14,r30,r31,r32,r33,r34,r35,r36,r37);
			CONV_L8_PART3C(l74,l75,l76,l77,l78,bm14,bm13,bm12,bm11,bm10,bm09,bm08,bm07,bm06,bm05,bm04,bm03,bm02,bm01, \
				bp00,bp01,bp02,bp03,bp04,bp05,bp06,bp07,bp08,bp09,bp10,bp11,bp12,bp13,bp14,r30,r31,r32,r33,r34,r35,r36,r37);
			_mm_store_ps(&wao3[n3p-32],r37);
			_mm_store_ps(&wao3[n3p-28],r36);
			_mm_store_ps(&wao3[n3p-24],r35);
			_mm_store_ps(&wao3[n3p-20],r34);
			_mm_store_ps(&wao3[n3p-16],r33);
			_mm_store_ps(&wao3[n3p-12],r32);
			_mm_store_ps(&wao3[n3p- 8],r31);
			_mm_store_ps(&wao3[n3p- 4],r30);
			CONV_L8_PART3(l34,l35,l36,l37,l38,cm14,cm13,cm12,cm11,cm10,cm09,cm08,cm07,cm06,cm05,cm04,cm03,cm02,cm01, \
				cp00,cp01,cp02,cp03,cp04,cp05,cp06,cp07,cp08,cp09,cp10,cp11,cp12,cp13,cp14,r70,r71,r72,r73,r74,r75,r76,r77);
			CONV_L8_PART3C(l74,l75,l76,l77,l78,em14,em13,em12,em11,em10,em09,em08,em07,em06,em05,em04,em03,em02,em01, \
				ep00,ep01,ep02,ep03,ep04,ep05,ep06,ep07,ep08,ep09,ep10,ep11,ep12,ep13,ep14,r70,r71,r72,r73,r74,r75,r76,r77);
			_mm_store_ps(&wao7[n3p-32],r77);
			_mm_store_ps(&wao7[n3p-28],r76);
			_mm_store_ps(&wao7[n3p-24],r75);
			_mm_store_ps(&wao7[n3p-20],r74);
			_mm_store_ps(&wao7[n3p-16],r73);
			_mm_store_ps(&wao7[n3p-12],r72);
			_mm_store_ps(&wao7[n3p- 8],r71);
			_mm_store_ps(&wao7[n3p- 4],r70);
		}
		//----------------------------------------------------------------------
		for(i3=0;i3<(n3p-4);i3+=4) {
			y_f[loffset1+(i3+0)*nseg3]=y_f[loffset1+(i3+0)*nseg3]+wao3[i3+0];
			y_f[loffset2+(i3+0)*nseg3]=y_f[loffset2+(i3+0)*nseg3]+wao7[i3+0];
			y_f[loffset1+(i3+1)*nseg3]=y_f[loffset1+(i3+1)*nseg3]+wao3[i3+1];
			y_f[loffset2+(i3+1)*nseg3]=y_f[loffset2+(i3+1)*nseg3]+wao7[i3+1];
			y_f[loffset1+(i3+2)*nseg3]=y_f[loffset1+(i3+2)*nseg3]+wao3[i3+2];
			y_f[loffset2+(i3+2)*nseg3]=y_f[loffset2+(i3+2)*nseg3]+wao7[i3+2];
			y_f[loffset1+(i3+3)*nseg3]=y_f[loffset1+(i3+3)*nseg3]+wao3[i3+3];
			y_f[loffset2+(i3+3)*nseg3]=y_f[loffset2+(i3+3)*nseg3]+wao7[i3+3];
		}
		for(;i3<n3t;i3++) {
			y_f[loffset1+i3*nseg3]=y_f[loffset1+i3*nseg3]+wao3[i3];
			y_f[loffset2+i3*nseg3]=y_f[loffset2+i3*nseg3]+wao7[i3];
		}
		} //end of if for n3t>0
		//----------------------------------------------------------------------
	} //end of loop over i1
	} //end of loop over i3
} // end of function convolut_ib_sse_fz_wws_www_
//****************************************************************************************

//|  1|  2|  3|  4|  5|  6|  7|  8|  9| 10| 11| 12| 13| 14| 15| 16| 17| 18| 19| 20| 21| 22| 23| 24|
//        |-14|-13|-12|-11|-10| -9| -8| -7| -6| -5| -4| -3| -2| -1|  0|  1|  2|  3|  4|  5|  6|  7|  8|  9| 10| 11| 12| 13| 14|

//|  1|  2|  3|  4|  5|  6|  7|  8|  9| 10| 11| 12| 13| 14| 15| 16| 17| 18| 19| 20| 21| 22| 23| 24|
//            |-14|-13|-12|-11|-10| -9| -8| -7| -6| -5| -4| -3| -2| -1|  0|  1|  2|  3|  4|  5|  6|  7|  8|  9| 10| 11| 12| 13| 14|

//|  1|  2|  3|  4|  5|  6|  7|  8|  9| 10| 11| 12| 13| 14| 15| 16| 17| 18| 19| 20| 21| 22| 23| 24|
//                |-14|-13|-12|-11|-10| -9| -8| -7| -6| -5| -4| -3| -2| -1|  0|  1|  2|  3|  4|  5|  6|  7|  8|  9| 10| 11| 12| 13| 14|

//|  1|  2|  3|  4|  5|  6|  7|  8|  9| 10| 11| 12| 13| 14| 15| 16| 17| 18| 19| 20| 21| 22| 23| 24|
//                    |-14|-13|-12|-11|-10| -9| -8| -7| -6| -5| -4| -3| -2| -1|  0|  1|  2|  3|  4|  5|  6|  7|  8|  9| 10| 11| 12| 13| 14|
