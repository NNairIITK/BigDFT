/*
 *  Brute force symmetry analyzer.
 *  This is actually C++ program, masquerading as a C one!
 *
 *  (C) 1996, 2003 S. Patchkovskii, Serguei.Patchkovskii@sympatico.ca
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * $Log: symmetry.c,v $
 * Revision 1.16  2003/04/04  13:05:03  patchkov
 * Revision 1.15  2000/01/25  16:47:17  patchkov
 * Revision 1.14  2000/01/25  16:39:08  patchkov
 * Revision 1.13  1996/05/24  12:32:08  ps
 * Revision 1.12  1996/05/23  16:10:47  ps
 * First reasonably stable version.
 *
 *
 * Modification for BigDFT by TD 2014/02/07
 */

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971694
#endif

#define	DIMENSION 3
#define MAXPARAM  7

typedef struct {
        int     type ;
        double  x[ DIMENSION ] ;
    } ATOM ;

/*
 *  All specific structures should have corresponding elements in the
 *  same position generic structure does.
 *
 *  Planes are characterized by the surface normal direction
 *  (taken in the direction *from* the coordinate origin)
 *  and distance from the coordinate origin to the plane
 *  in the direction of the surface normal.
 *
 *  Inversion is characterized by location of the inversion center.
 *
 *  Rotation is characterized by a vector (distance+direction) from the origin
 *  to the rotation axis, axis direction and rotation order. Rotations
 *  are in the clockwise direction looking opposite to the direction
 *  of the axis. Note that this definition of the rotation axis
 *  is *not* unique, since an arbitrary multiple of the axis direction
 *  can be added to the position vector without changing actual operation.
 *
 *  Mirror rotation is defined by the same parameters as normal rotation,
 *  but the origin is now unambiguous since it defines the position of the
 *  plane associated with the axis.
 *
 */

typedef struct _SYMMETRY_ELEMENT_ {
        void    (*transform_atom)( struct _SYMMETRY_ELEMENT_ *el, ATOM *from, ATOM *to ) ;
        int *   transform ;     /*   Correspondence table for the transformation         */
        int     order ;         /*   Applying transformation this many times is identity */
        int     nparam ;        /*   4 for inversion and planes, 7 for axes              */
        double  maxdev ;        /*   Larges error associated with the element            */
        double  distance ;
        double  normal[ DIMENSION ] ;
        double  direction[ DIMENSION ] ;
    } SYMMETRY_ELEMENT ;

typedef struct {
        char *  group_name ;        /* Canonical group name                              */
        char *  symmetry_code ;     /* Group symmetry code                               */
        int     (*check)( void ) ;  /* Additional verification routine, not used         */
    } POINT_GROUP ;

/* Global variables */
typedef struct {
        double                 ToleranceSame         ;
        double                 TolerancePrimary      ;
        double                 ToleranceFinal        ;
        double                 MaxOptStep            ;
        double                 MinOptStep            ;
        double                 GradientStep          ;
        double                 OptChangeThreshold    ;
        double                 CenterOfSomething[ DIMENSION ] ;
        double *               DistanceFromCenter    ;
        int                    verbose               ;
        int                    MaxOptCycles          ;
        int                    OptChangeHits         ;
        int                    MaxAxisOrder          ;
        int                    AtomsCount            ;
        ATOM *                 Atoms                 ;
        int                    PlanesCount           ;
        SYMMETRY_ELEMENT **    Planes                ;
        SYMMETRY_ELEMENT *     MolecularPlane        ;
        int                    InversionCentersCount ;
        SYMMETRY_ELEMENT **    InversionCenters      ;
        int                    NormalAxesCount       ;
        SYMMETRY_ELEMENT **    NormalAxes            ;
        int                    ImproperAxesCount     ;
        SYMMETRY_ELEMENT **    ImproperAxes          ;
        int *                  NormalAxesCounts      ;
        int *                  ImproperAxesCounts    ;
        int                    BadOptimization       ;
        char *                 SymmetryCode          ;
} SYMMETRIES ;


/* Statistics */
typedef struct {
        long                   StatTotal             ;
        long                   StatEarly             ;
        long                   StatPairs             ;
        long                   StatDups              ;
        long                   StatOrder             ;
        long                   StatOpt               ;
        long                   StatAccept            ;
} STATISTICS ;

