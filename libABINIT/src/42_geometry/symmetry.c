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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "symmetry.h"

/*
 *    Point groups I know about
 */
int true(void){ return 1 ; }
POINT_GROUP            PointGroups[]         = {
    {  "C1",    "",                                                          true  },
    {  "Cs",    "(sigma) ",                                                  true  },
    {  "Ci",    "(i) ",                                                      true  },
    {  "C2",    "(C2) ",                                                     true  },
    {  "C3",    "(C3) ",                                                     true  },
    {  "C4",    "(C4) (C2) ",                                                true  },
    {  "C5",    "(C5) ",                                                     true  },
    {  "C6",    "(C6) (C3) (C2) ",                                           true  },
    {  "C7",    "(C7) ",                                                     true  },
    {  "C8",    "(C8) (C4) (C2) ",                                           true  },
    {  "D2",    "3*(C2) ",                                                   true  },
    {  "D3",    "(C3) 3*(C2) ",                                              true  },
    {  "D4",    "(C4) 5*(C2) ",                                              true  },
    {  "D5",    "(C5) 5*(C2) ",                                              true  },
    {  "D6",    "(C6) (C3) 7*(C2) ",                                         true  },
    {  "D7",    "(C7) 7*(C2) ",                                              true  },
    {  "D8",    "(C8) (C4) 9*(C2) ",                                         true  },
    {  "C2v",   "(C2) 2*(sigma) ",                                           true  },
    {  "C3v",   "(C3) 3*(sigma) ",                                           true  },
    {  "C4v",   "(C4) (C2) 4*(sigma) ",                                      true  },
    {  "C5v",   "(C5) 5*(sigma) ",                                           true  },
    {  "C6v",   "(C6) (C3) (C2) 6*(sigma) ",                                 true  },
    {  "C7v",   "(C7) 7*(sigma) ",                                           true  },
    {  "C8v",   "(C8) (C4) (C2) 8*(sigma) ",                                 true  },
    {  "C2h",   "(i) (C2) (sigma) ",                                         true  },
    {  "C3h",   "(C3) (S3) (sigma) ",                                        true  },
    {  "C4h",   "(i) (C4) (C2) (S4) (sigma) ",                               true  },
    {  "C5h",   "(C5) (S5) (sigma) ",                                        true  },
    {  "C6h",   "(i) (C6) (C3) (C2) (S6) (S3) (sigma) ",                     true  },
    {  "C7h",   "(C7) (S7) (sigma) ",                                        true  },
    {  "C8h",   "(i) (C8) (C4) (C2) (S8) (S4) (sigma) ",                     true  },
    {  "D2h",   "(i) 3*(C2) 3*(sigma) ",                                     true  },
    {  "D3h",   "(C3) 3*(C2) (S3) 4*(sigma) ",                               true  },
    {  "D4h",   "(i) (C4) 5*(C2) (S4) 5*(sigma) ",                           true  },
    {  "D5h",   "(C5) 5*(C2) (S5) 6*(sigma) ",                               true  },
    {  "D6h",   "(i) (C6) (C3) 7*(C2) (S6) (S3) 7*(sigma) ",                 true  },
    {  "D7h",   "(C7) 7*(C2) (S7) 8*(sigma) ",                               true  },
    {  "D8h",   "(i) (C8) (C4) 9*(C2) (S8) (S4) 9*(sigma) ",                 true  },
    {  "D2d",   "3*(C2) (S4) 2*(sigma) ",                                    true  },
    {  "D3d",   "(i) (C3) 3*(C2) (S6) 3*(sigma) ",                           true  },
    {  "D4d",   "(C4) 5*(C2) (S8) 4*(sigma) ",                               true  },
    {  "D5d",   "(i) (C5) 5*(C2) (S10) 5*(sigma) ",                          true  },
    {  "D6d",   "(C6) (C3) 7*(C2) (S12) (S4) 6*(sigma) ",                    true  },
    {  "D7d",   "(i) (C7) 7*(C2) (S14) 7*(sigma) ",                          true  },
    {  "D8d",   "(C8) (C4) 9*(C2) (S16) 8*(sigma) ",                         true  },
    {  "S4",    "(C2) (S4) ",                                                true  },
    {  "S6",    "(i) (C3) (S6) ",                                            true  },
    {  "S8",    "(C4) (C2) (S8) ",                                           true  },
    {  "T",     "4*(C3) 3*(C2) ",                                            true  },
    {  "Th",    "(i) 4*(C3) 3*(C2) 4*(S6) 3*(sigma) ",                       true  },
    {  "Td",    "4*(C3) 3*(C2) 3*(S4) 6*(sigma) ",                           true  },
    {  "O",     "3*(C4) 4*(C3) 9*(C2) ",                                     true  },
    {  "Oh",    "(i) 3*(C4) 4*(C3) 9*(C2) 4*(S6) 3*(S4) 9*(sigma) ",         true  },
    {  "Cinfv", "(Cinf) (sigma) ",                                           true  },
    {  "Dinfh", "(i) (Cinf) (C2) 2*(sigma) ",                                true  },
    {  "I",     "6*(C5) 10*(C3) 15*(C2) ",                                   true  },
    {  "Ih",    "(i) 6*(C5) 10*(C3) 15*(C2) 6*(S10) 10*(S6) 15*(sigma) ",    true  },
    {  "Kh",    "(i) (Cinf) (sigma) ",                                       true  },
    } ;

#define PointGroupsCount (sizeof(PointGroups)/sizeof(POINT_GROUP))
char *                 PointGroupRejectionReason = NULL ;

/* Initialize the struct SYMMETRIES */
void
init_symmetries( SYMMETRIES *gsym ) {
        gsym->ToleranceSame         = 1e-3 ;
        gsym->TolerancePrimary      = 5e-2 ;
        gsym->ToleranceFinal        = 1e-4 ;
        gsym->MaxOptStep            = 5e-1 ;
        gsym->MinOptStep            = 1e-7 ;
        gsym->GradientStep          = 1e-7 ;
        gsym->OptChangeThreshold    = 1e-10 ;
        gsym->DistanceFromCenter    = NULL ;
        gsym->verbose               = 0 ;
        gsym->MaxOptCycles          = 200 ;
        gsym->OptChangeHits         = 5 ;
        gsym->MaxAxisOrder          = 20 ;
        gsym->AtomsCount            = 0 ;
        gsym->Atoms                 = NULL ;
        gsym->PlanesCount           = 0 ;
        gsym->Planes                = NULL ;
        gsym->MolecularPlane        = NULL ;
        gsym->InversionCentersCount = 0 ;
        gsym->InversionCenters      = NULL ;
        gsym->NormalAxesCount       = 0 ;
        gsym->NormalAxes            = NULL ;
        gsym->ImproperAxesCount     = 0 ;
        gsym->ImproperAxes          = NULL ;
        gsym->NormalAxesCounts      = NULL ;
        gsym->ImproperAxesCounts    = NULL ;
        gsym->BadOptimization       = 0 ;
        gsym->SymmetryCode          = "" ;
} ;


/* Initialize the Statistics */
void
init_statistics( STATISTICS *stat ) {
        stat->StatTotal             = 0 ;
        stat->StatEarly             = 0 ;
        stat->StatPairs             = 0 ;
        stat->StatDups              = 0 ;
        stat->StatOrder             = 0 ;
        stat->StatOpt               = 0 ;
        stat->StatAccept            = 0 ;
} ;

/*
 *   Generic functions
 */

double
pow2( double x )
{
return x * x ;
}

int
establish_pairs( SYMMETRY_ELEMENT *elem, SYMMETRIES *gsym )
{
        int               i, j, k, best_j ;
        char *            atom_used = calloc( gsym->AtomsCount, 1 ) ;
        double            distance, best_distance ;
        ATOM              symmetric ;

if( atom_used == NULL ){
    fprintf( stderr, "Out of memory for tagging array in establish_pairs()\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 0 ; i < gsym->AtomsCount ; i++ ){
    if( elem->transform[i] >= gsym->AtomsCount ){ /* No symmetric atom yet          */
        if( gsym->verbose > 2 ) printf( "        looking for a pair for %d\n", i ) ;
        elem->transform_atom( elem, gsym->Atoms+i, &symmetric ) ;
        if( gsym->verbose > 2 ) printf( "        new coordinates are: (%g,%g,%g)\n", 
                              symmetric.x[0], symmetric.x[1], symmetric.x[2] ) ;
        best_j        = i ;
        best_distance = 2*gsym->TolerancePrimary ;/* Performance value we'll reject */
        for( j = 0 ; j < gsym->AtomsCount ; j++ ){
            if( gsym->Atoms[j].type != symmetric.type || atom_used[j] )
                continue ;
            for( k = 0, distance = 0 ; k < DIMENSION ; k++ ){
                distance += pow2( symmetric.x[k] - gsym->Atoms[j].x[k] ) ;
                }
            distance = sqrt( distance ) ;
            if( gsym->verbose > 2 ) printf( "        distance to %d is %g\n", j, distance ) ;
            if( distance < best_distance ){
                best_j        = j ;
                best_distance = distance ;
                }
            }
        if( best_distance > gsym->TolerancePrimary ){ /* Too bad, there is no symmetric atom */
            if( gsym->verbose > 0 ) 
                printf( "        no pair for atom %d - best was %d with err = %g\n", i, best_j, best_distance ) ;
            free( atom_used ) ;
            return -1 ;
            }
        elem->transform[i] = best_j ;
        atom_used[best_j]  = 1 ;
        if( gsym->verbose > 1 ) printf( "        atom %d transforms to the atom %d, err = %g\n", i, best_j, best_distance ) ;
        }
    }
free( atom_used ) ;
return 0 ;
}

int
check_transform_order( SYMMETRY_ELEMENT *elem, SYMMETRIES *gsym )
{
        int             i, j, k ;
        void            rotate_reflect_atom( SYMMETRY_ELEMENT *, ATOM *, ATOM *) ;

for( i = 0 ; i < gsym->AtomsCount ; i++ ){
    if( elem->transform[i] == i )   /* Identity transform is Ok for any order */
        continue ;
    if( elem->transform_atom == rotate_reflect_atom ){
        j = elem->transform[i] ;
        if( elem->transform[j] == i )
            continue ; /* Second-order transform is Ok for improper axis */
        }
    for( j = elem->order - 1, k = elem->transform[i] ; j > 0 ; j--, k = elem->transform[k] ){
        if( k == i ){
            if( gsym->verbose > 0 ) printf( "        transform looped %d steps too early from atom %d\n", j, i ) ;
            return -1 ;
            }
        }
    if( k != i && elem->transform_atom == rotate_reflect_atom ){
        /* For improper axes, the complete loop may also take twice the order */
        for( j = elem->order ; j > 0 ; j--, k = elem->transform[k] ){
            if( k == i ){
                if( gsym->verbose > 0 ) printf( "        (improper) transform looped %d steps too early from atom %d\n", j, i ) ;
                return -1 ;
                }
            }
        }
    if( k != i ){
        if( gsym->verbose > 0 ) printf( "        transform failed to loop after %d steps from atom %d\n", elem->order, i ) ;
        return -1 ;
        }
    }
return 0 ;
}

int
same_transform( SYMMETRY_ELEMENT *a, SYMMETRY_ELEMENT *b, SYMMETRIES *gsym)
{
        int               i, j ;
        int               code ;

if( ( a->order != b->order ) || ( a->nparam != b->nparam ) || ( a->transform_atom != b->transform_atom ) )
    return 0 ;
for( i = 0, code = 1 ; i < gsym->AtomsCount ; i++ ){
    if( a->transform[i] != b->transform[i] ){
        code = 0 ;
        break ;
        }
    }
if( code == 0 && a->order > 2 ){  /* b can also be a reverse transformation for a */
    for( i = 0 ; i < gsym->AtomsCount ; i++ ){
        j = a->transform[i] ;
        if( b->transform[j] != i )
            return 0 ;
        }
    return 1 ;
    }
return code ;
}

SYMMETRY_ELEMENT *
alloc_symmetry_element( SYMMETRIES *gsym )
{
        SYMMETRY_ELEMENT * elem = calloc( 1, sizeof( SYMMETRY_ELEMENT ) ) ;
        int                i ;

if( elem == NULL ){
    fprintf( stderr, "Out of memory allocating symmetry element\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
elem->transform = calloc( gsym->AtomsCount, sizeof( int ) ) ;
if( elem->transform == NULL ){
    fprintf( stderr, "Out of memory allocating transform table for symmetry element\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 0 ; i < gsym->AtomsCount ; i++ ){
    elem->transform[i] = gsym->AtomsCount + 1 ; /* An impossible value */
    }
return elem ;
}

/* Free SYMMETRY_ELEMENT */
void
destroy_symmetry_element( SYMMETRY_ELEMENT *elem )
{
if( elem != NULL ){
    if( elem->transform != NULL )
        free( elem->transform ) ;
    free( elem ) ;
    }
}

int
check_transform_quality( SYMMETRY_ELEMENT *elem, SYMMETRIES *gsym )
{
        int               i, j, k ;
        ATOM              symmetric ;
        double            r, max_r ;

for( i = 0, max_r = 0 ; i < gsym->AtomsCount ; i++ ){
    j = elem->transform[i] ;
    elem->transform_atom( elem, gsym->Atoms + i, &symmetric ) ;
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += pow2( symmetric.x[k] - gsym->Atoms[j].x[k] ) ;
        }
    r = sqrt( r ) ;
    if( r > gsym->ToleranceFinal ){
        if( gsym->verbose > 0 ) printf( "        distance to symmetric atom (%g) is too big for %d\n", r, i ) ;
        return -1 ;
        }
    if( r > max_r ) max_r = r ;
    }
elem->maxdev = max_r ;
return 0 ;
}

double
eval_optimization_target_function( SYMMETRY_ELEMENT *elem, int *finish, SYMMETRIES *gsym )
{
        int               i, j, k ;
        ATOM              symmetric ;
        double            target, r, maxr ;

if( elem->nparam >= 4 ){
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += elem->normal[k]*elem->normal[k] ;
        }
    r = sqrt( r ) ;
    if( r < gsym->ToleranceSame ){
        fprintf( stderr, "Normal collapced!\n" ) ;
        exit( EXIT_FAILURE ) ;
        }
    for( k = 0 ; k < DIMENSION ; k++ ){
        elem->normal[k] /= r ;
        }
    if( elem->distance < 0 ){
        elem->distance = -elem->distance ;
        for( k = 0 ; k < DIMENSION ; k++ ){
            elem->normal[k] = -elem->normal[k] ;
            }
        }
    }
if( elem->nparam >= 7 ){
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += elem->direction[k]*elem->direction[k] ;
        }
    r = sqrt( r ) ;
    if( r < gsym->ToleranceSame ){
        fprintf( stderr, "Direction collapced!\n" ) ;
        exit( EXIT_FAILURE ) ;
        }
    for( k = 0 ; k < DIMENSION ; k++ ){
        elem->direction[k] /= r ;
        }
    }
for( i = 0, target = maxr = 0 ; i < gsym->AtomsCount ; i++ ){
    elem->transform_atom( elem, gsym->Atoms + i, &symmetric ) ;
    j = elem->transform[i] ;
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += pow2( gsym->Atoms[j].x[k] - symmetric.x[k] ) ;
        }
    if( r > maxr ) maxr = r ;
    target += r ;
    }
if( finish != NULL ){
    *finish = 0 ;
    if( sqrt( maxr ) < gsym->ToleranceFinal )
        *finish = 1 ;
    }
return target ;
}

void
get_params( SYMMETRY_ELEMENT *elem, double values[] )
{
memcpy( values, &elem->distance, elem->nparam * sizeof( double ) ) ;
}

void
set_params( SYMMETRY_ELEMENT *elem, double values[] )
{
memcpy( &elem->distance, values, elem->nparam * sizeof( double ) ) ;
}

void
optimize_transformation_params( SYMMETRY_ELEMENT *elem, SYMMETRIES *gsym )
{
        double            values[ MAXPARAM ] ;
        double            grad  [ MAXPARAM ] ;
        double            force [ MAXPARAM ] ;
        double            step  [ MAXPARAM ] ;
        double            f, fold, fnew, fnew2, fdn, fup, snorm ;
        double            a, b, x ;
        int               vars  = elem->nparam ;
        int               cycle = 0 ;
        int               i, finish ;
        int               hits = 0 ;

if( vars > MAXPARAM ){
    fprintf( stderr, "Catastrophe in optimize_transformation_params()!\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
f = 0 ;
do {
    fold = f ;
    f    = eval_optimization_target_function( elem, &finish, gsym ) ;
    /* Evaluate function, gradient and diagonal force constants */
    if( gsym->verbose > 1 ) printf( "            function value = %g\n", f ) ;
    if( finish ){
        if( gsym->verbose > 1 ) printf( "        function value is small enough\n" ) ;
        break ;
        }
    if( cycle > 0 ){
        if( fabs( f-fold ) > gsym->OptChangeThreshold )
             hits = 0 ;
        else hits++ ;
        if( hits >= gsym->OptChangeHits ){
            if( gsym->verbose > 1 ) printf( "        no progress is made, stop optimization\n" ) ;
            break ;
            }
        }
    get_params( elem, values ) ;
    for( i = 0 ; i < vars ; i++ ){
        values[i] -= gsym->GradientStep ;
        set_params( elem, values ) ;
        fdn        = eval_optimization_target_function( elem, NULL, gsym ) ;
        values[i] += 2*gsym->GradientStep ;
        set_params( elem, values ) ;
        fup        = eval_optimization_target_function( elem, NULL, gsym ) ;
        values[i] -= gsym->GradientStep ;
        grad[i]    = ( fup - fdn ) / ( 2 * gsym->GradientStep ) ;
        force[i]   = ( fup + fdn - 2*f ) / ( gsym->GradientStep * gsym->GradientStep ) ;
        if( gsym->verbose > 1 ) printf( "        i = %d, grad = %12.6e, force = %12.6e\n", i, grad[i], force[i] ) ;
        }
    /* Do a quasy-Newton step */
    for( i = 0, snorm = 0 ; i < vars ; i++ ){
        if( force[i] <  0   ) force[i] = -force[i] ;
        if( force[i] < 1e-3 ) force[i] = 1e-3 ;
        if( force[i] > 1e3  ) force[i] = 1e3 ;
        step[i] = - grad[i]/force[i] ;
        snorm += step[i] * step[i] ;
        }
    snorm = sqrt( snorm ) ;
    if( snorm > gsym->MaxOptStep ){ /* Renormalize step */
        for( i = 0 ; i < vars ; i++ )
            step[i] *= gsym->MaxOptStep/snorm ;
        snorm = gsym->MaxOptStep ;
        }
    do {
        for( i = 0 ; i < vars ; i++ ){
            values[i] += step[i] ;
            }
        set_params( elem, values ) ;
        fnew = eval_optimization_target_function( elem, NULL, gsym ) ;
        if( fnew < f )
            break ;
        for( i = 0 ; i < vars ; i++ ){
            values[i] -= step[i] ;
            step  [i] /= 2 ;
            }
        set_params( elem, values ) ;
        snorm /= 2 ;
        } while( snorm > gsym->MinOptStep ) ;
        if( (snorm > gsym->MinOptStep) && (snorm < gsym->MaxOptStep / 2) ){  /* try to do quadratic interpolation */
            for( i = 0 ; i < vars ; i++ )
                values[i] += step[i] ;
            set_params( elem, values ) ;
            fnew2 = eval_optimization_target_function( elem, NULL, gsym ) ;
            if( gsym->verbose > 1 ) printf( "        interpolation base points: %g, %g, %g\n", f, fnew, fnew2 ) ;
            for( i = 0 ; i < vars ; i++ )
                values[i] -= 2*step[i] ;
            a     = ( 4*f - fnew2 - 3*fnew ) / 2 ;
            b     = ( f + fnew2 - 2*fnew ) / 2 ;
            if( gsym->verbose > 1 ) printf( "        linear interpolation coefficients %g, %g\n", a, b ) ;
            if( b > 0 ){
                x = -a/(2*b) ;
                if( x > 0.2 && x < 1.8 ){
                    if( gsym->verbose > 1 ) printf( "        interpolated: %g\n", x ) ;
                    for( i = 0 ; i < vars ; i++ )
                        values[i] += x*step[i] ;
                    }
                else b = 0 ;
                }
            if( b <= 0 ){
                if( fnew2 < fnew ){
                    for( i = 0 ; i < vars ; i++ )
                        values[i] += 2*step[i] ;
                    }
                else {
                    for( i = 0 ; i < vars ; i++ )
                        values[i] += step[i] ;
                    }
                }
            set_params( elem, values ) ;
            }
    } while( snorm > gsym->MinOptStep && ++cycle < gsym->MaxOptCycles ) ;
f = eval_optimization_target_function( elem, NULL, gsym ) ;
if( cycle >= gsym->MaxOptCycles ) gsym->BadOptimization = 1 ;
if( gsym->verbose > 0 ) {
    if( cycle >= gsym->MaxOptCycles )
        printf( "        maximum number of optimization cycles made\n" ) ;
        printf( "        optimization completed after %d cycles with f = %g\n", cycle, f ) ;
    }
}

int
refine_symmetry_element( SYMMETRY_ELEMENT *elem, int build_table, SYMMETRIES *gsym, STATISTICS *stat )
{
        int               i ;


if( build_table && (establish_pairs( elem, gsym ) < 0) ){
    stat->StatPairs++ ;
    if( gsym->verbose > 0 ) printf( "        no transformation correspondence table can be constructed\n" ) ;
    return -1 ;
    }
for( i = 0 ; i < gsym->PlanesCount ; i++ ){
    if( same_transform( gsym->Planes[i], elem, gsym ) ){
        stat->StatDups++ ;
        if( gsym->verbose > 0 ) printf( "        transformation is identical to plane %d\n", i ) ;
        return -1 ;
        }
    }
for( i = 0 ; i < gsym->InversionCentersCount ; i++ ){
    if( same_transform( gsym->InversionCenters[i], elem, gsym ) ){
        stat->StatDups++ ;
        if( gsym->verbose > 0 ) printf( "        transformation is identical to inversion center %d\n", i ) ;
        return -1 ;
        }
    }
for( i = 0 ; i < gsym->NormalAxesCount ; i++ ){
    if( same_transform( gsym->NormalAxes[i], elem, gsym ) ){
        stat->StatDups++ ;
        if( gsym->verbose > 0 ) printf( "        transformation is identical to normal axis %d\n", i ) ;
        return -1 ;
        }
    }
for( i = 0 ; i < gsym->ImproperAxesCount ; i++ ){
    if( same_transform( gsym->ImproperAxes[i], elem, gsym ) ){
        stat->StatDups++ ;
        if( gsym->verbose > 0 ) printf( "        transformation is identical to improper axis %d\n", i ) ;
        return -1 ;
        }
    }
if( check_transform_order( elem, gsym ) < 0 ){
    stat->StatOrder++ ;
    if( gsym->verbose > 0 ) printf( "        incorrect transformation order\n" ) ;
    return -1 ;
    }
optimize_transformation_params( elem, gsym ) ;
if( check_transform_quality( elem, gsym ) < 0 ){
    stat->StatOpt++ ;
    if( gsym->verbose > 0 ) printf( "        refined transformation does not pass the numeric threshold\n" ) ;
    return -1 ;
    }
stat->StatAccept++ ;
return 0 ;
}

/*
 *   Plane-specific functions
 */

void
mirror_atom( SYMMETRY_ELEMENT *plane, ATOM *from, ATOM *to )
{
        int                i ;
        double             r ;

for( i = 0, r = plane->distance ; i < DIMENSION ; i++ ){
    r -= from->x[i] * plane->normal[i] ;
    }
to->type = from->type ;
for( i = 0 ; i < DIMENSION ; i++ ){
    to->x[i] = from->x[i] + 2*r*plane->normal[i] ;
    }
}

SYMMETRY_ELEMENT * 
init_mirror_plane( int i, int j, SYMMETRIES *gsym, STATISTICS *stat )
{
        SYMMETRY_ELEMENT * plane = alloc_symmetry_element(gsym) ;
        double             dx[ DIMENSION ], midpoint[ DIMENSION ], rab, r ;
        int                k ;

if( gsym->verbose > 0 ) printf( "Trying mirror plane for atoms %d,%d\n", i, j ) ;
stat->StatTotal++ ;
plane->transform_atom = mirror_atom ;
plane->order          = 2 ;
plane->nparam         = 4 ;
for( k = 0, rab = 0 ; k < DIMENSION ; k++ ){
    dx[k]       = gsym->Atoms[i].x[k] - gsym->Atoms[j].x[k] ;
    midpoint[k] = ( gsym->Atoms[i].x[k] + gsym->Atoms[j].x[k] ) / 2.0 ;
    rab        += dx[k]*dx[k] ;
    }
rab = sqrt(rab) ;
if( rab < gsym->ToleranceSame ){
    fprintf( stderr, "gsym->Atoms %d and %d coincide (r = %g)\n", i, j, rab ) ;
    exit( EXIT_FAILURE ) ;
    }
for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
    plane->normal[k] = dx[k]/rab ;
    r += midpoint[k]*plane->normal[k] ;
    }
if( r < 0 ){  /* Reverce normal direction, distance is always positive! */
    r = -r ;
    for( k = 0 ; k < DIMENSION ; k++ ){
        plane->normal[k] = -plane->normal[k] ;
        }
    }
plane->distance = r ;
if( gsym->verbose > 0 ) printf( "    initial plane is at %g from the origin\n", r ) ;
if( refine_symmetry_element( plane, 1, gsym, stat ) < 0 ){
    if( gsym->verbose > 0 ) printf( "    refinement failed for the plane\n" ) ;
    destroy_symmetry_element( plane ) ;
    return NULL ;
    }
return plane ;
}

SYMMETRY_ELEMENT *
init_ultimate_plane( SYMMETRIES *gsym, STATISTICS *stat )
{
        SYMMETRY_ELEMENT * plane = alloc_symmetry_element(gsym) ;
        double             d0[ DIMENSION ], d1[ DIMENSION ], d2[ DIMENSION ] ;
        double             p[ DIMENSION ] ;
        double             r, s0, s1, s2 ;
        double *           d ;
        int                i, j, k ;

if( gsym->verbose > 0 ) printf( "Trying whole-molecule mirror plane\n" ) ;
stat->StatTotal++ ;
plane->transform_atom = mirror_atom ;
plane->order          = 1 ;
plane->nparam         = 4 ;
for( k = 0 ; k < DIMENSION ; k++ )
    d0[k] = d1[k] = d2[k] = 0 ;
d0[0] = 1 ; d1[1] = 1 ; d2[2] = 1 ;
for( i = 1 ; i < gsym->AtomsCount ; i++ ){
    for( j = 0 ; j < i ; j++ ){
        for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
            p[k] = gsym->Atoms[i].x[k] - gsym->Atoms[j].x[k] ;
            r   += p[k]*p[k] ;
            }
        r = sqrt(r) ;
        for( k = 0, s0=s1=s2=0 ; k < DIMENSION ; k++ ){
            p[k] /= r ;
            s0   += p[k]*d0[k] ;
            s1   += p[k]*d1[k] ;
            s2   += p[k]*d2[k] ;
            }
        for( k = 0 ; k < DIMENSION ; k++ ){
            d0[k] -= s0*p[k] ;
            d1[k] -= s1*p[k] ;
            d2[k] -= s2*p[k] ;
            }
        }
    }
for( k = 0, s0=s1=s2=0 ; k < DIMENSION ; k++ ){
    s0 += d0[k] ;
    s1 += d1[k] ;
    s2 += d2[k] ;
    }
d = NULL ;
if( s0 >= s1 && s0 >= s2 ) d = d0 ;
if( s1 >= s0 && s1 >= s2 ) d = d1 ;
if( s2 >= s0 && s2 >= s1 ) d = d2 ;
if( d == NULL ){
    fprintf( stderr, "Catastrophe in init_ultimate_plane(): %g, %g and %g have no ordering!\n", s0, s1, s2 ) ;
    exit( EXIT_FAILURE ) ;
    }
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += d[k]*d[k] ;
r = sqrt(r) ;
if( r > 0 ){
    for( k = 0 ; k < DIMENSION ; k++ )
        plane->normal[k] = d[k]/r ;
    }
else {
    for( k = 1 ; k < DIMENSION ; k++ )
        plane->normal[k] = 0 ;
    plane->normal[0] = 1 ;
    }
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += gsym->CenterOfSomething[k]*plane->normal[k] ;
plane->distance = r ;
for( k = 0 ; k < gsym->AtomsCount ; k++ )
    plane->transform[k] = k ;
if( refine_symmetry_element( plane, 0, gsym, stat ) < 0 ){
    if( gsym->verbose > 0 ) printf( "    refinement failed for the plane\n" ) ;
    destroy_symmetry_element( plane ) ;
    return NULL ;
    }
return plane ;
}
/*
 *   Inversion-center specific functions
 */
void
invert_atom( SYMMETRY_ELEMENT *center, ATOM *from, ATOM *to )
{
        int                i ;

to->type = from->type ;
for( i = 0 ; i < DIMENSION ; i++ ){
    to->x[i] = 2*center->distance*center->normal[i] - from->x[i] ;
    }
}

SYMMETRY_ELEMENT *
init_inversion_center( SYMMETRIES *gsym, STATISTICS *stat )
{
        SYMMETRY_ELEMENT * center = alloc_symmetry_element(gsym) ;
        int                k ;
        double             r ;

if( gsym->verbose > 0 ) printf( "Trying inversion center at the center of something\n" ) ;
stat->StatTotal++ ;
center->transform_atom = invert_atom ;
center->order          = 2 ;
center->nparam         = 4 ;
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += gsym->CenterOfSomething[k]*gsym->CenterOfSomething[k] ;
r = sqrt(r) ;
if( r > 0 ){
    for( k = 0 ; k < DIMENSION ; k++ )
        center->normal[k] = gsym->CenterOfSomething[k]/r ;
    }
else {
    center->normal[0] = 1 ;
    for( k = 1 ; k < DIMENSION ; k++ )
        center->normal[k] = 0 ;
    }
center->distance = r ;
if( gsym->verbose > 0 ) printf( "    initial inversion center is at %g from the origin\n", r ) ;
if( refine_symmetry_element( center, 1, gsym, stat ) < 0 ){
    if( gsym->verbose > 0 ) printf( "    refinement failed for the inversion center\n" ) ;
    destroy_symmetry_element( center ) ;
    return NULL ;
    }
return center ;
}

/*
 *   Normal rotation axis-specific routines.
 */
void
rotate_atom( SYMMETRY_ELEMENT *axis, ATOM *from, ATOM *to )
{
        double             x[3], y[3], a[3], b[3], c[3] ;
        double             angle = axis->order ? 2*M_PI/axis->order : 1.0 ;
        double             a_sin = sin( angle ) ;
        double             a_cos = cos( angle ) ;
        double             dot ;
        int                i ;

if( DIMENSION != 3 ){
    fprintf( stderr, "Catastrophe in rotate_atom!\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 0 ; i < 3 ; i++ )
    x[i] = from->x[i] - axis->distance * axis->normal[i] ;
for( i = 0, dot = 0 ; i < 3 ; i++ )
    dot += x[i] * axis->direction[i] ;
for( i = 0 ; i < 3 ; i++ )
    a[i] = axis->direction[i] * dot ;
for( i = 0 ; i < 3 ; i++ )
    b[i] = x[i] - a[i] ;
c[0] = b[1]*axis->direction[2] - b[2]*axis->direction[1] ;
c[1] = b[2]*axis->direction[0] - b[0]*axis->direction[2] ;
c[2] = b[0]*axis->direction[1] - b[1]*axis->direction[0] ;
for( i = 0 ; i < 3 ; i++ )
    y[i] = a[i] + b[i]*a_cos + c[i]*a_sin ;
for( i = 0 ; i < 3 ; i++ )
    to->x[i] = y[i] + axis->distance * axis->normal[i] ;
to->type = from->type ;
}

SYMMETRY_ELEMENT *
init_ultimate_axis( SYMMETRIES *gsym, STATISTICS *stat )
{
        SYMMETRY_ELEMENT * axis = alloc_symmetry_element(gsym) ;
        double             dir[ DIMENSION ], rel[ DIMENSION ] ;
        double             s ;
        int                i, k ;

if( gsym->verbose > 0 ) printf( "Trying infinity axis\n" ) ;
stat->StatTotal++ ;
axis->transform_atom = rotate_atom ;
axis->order          = 0 ;
axis->nparam         = 7 ;
for( k = 0 ; k < DIMENSION ; k++ )
    dir[k] = 0 ;
for( i = 0 ; i < gsym->AtomsCount ; i++ ){
    for( k = 0, s = 0 ; k < DIMENSION ; k++ ){
        rel[k] = gsym->Atoms[i].x[k] - gsym->CenterOfSomething[k] ;
        s     += rel[k]*dir[k] ;
        }
    if( s >= 0 )
         for( k = 0 ; k < DIMENSION ; k++ )
             dir[k] += rel[k] ;
    else for( k = 0 ; k < DIMENSION ; k++ )
             dir[k] -= rel[k] ;
    }
for( k = 0, s = 0 ; k < DIMENSION ; k++ )
    s += pow2( dir[k] ) ;
s = sqrt(s) ;
if( s > 0 )
     for( k = 0 ; k < DIMENSION ; k++ )
         dir[k] /= s ;
else dir[0] = 1 ;
for( k = 0 ; k < DIMENSION ; k++ )
    axis->direction[k] = dir[k] ;
for( k = 0, s = 0 ; k < DIMENSION ; k++ )
    s += pow2( gsym->CenterOfSomething[k] ) ;
s = sqrt(s) ;
if( s > 0 )
    for( k = 0 ; k < DIMENSION ; k++ )
        axis->normal[k] = gsym->CenterOfSomething[k]/s ;
else {
    for( k = 1 ; k < DIMENSION ; k++ )
        axis->normal[k] = 0 ;
    axis->normal[0] = 1 ;
    }
axis->distance = s ;
for( k = 0 ; k < gsym->AtomsCount ; k++ )
    axis->transform[k] = k ;
if( refine_symmetry_element( axis, 0, gsym, stat ) < 0 ){
    if( gsym->verbose > 0 ) printf( "    refinement failed for the infinity axis\n" ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}


SYMMETRY_ELEMENT *
init_c2_axis( int i, int j, double support[ DIMENSION ], SYMMETRIES *gsym, STATISTICS *stat )
{
        SYMMETRY_ELEMENT * axis ;
        int                k ;
        double             ris, rjs ;
        double             r, center[ DIMENSION ] ;

if( gsym->verbose > 0 ) 
    printf( "Trying c2 axis for the pair (%d,%d) with the support (%g,%g,%g)\n", 
             i, j, support[0], support[1], support[2] ) ;
stat->StatTotal++ ;
/* First, do a quick sanity check */
for( k = 0, ris = rjs = 0 ; k < DIMENSION ; k++ ){
    ris += pow2( gsym->Atoms[i].x[k] - support[k] ) ;
    rjs += pow2( gsym->Atoms[j].x[k] - support[k] ) ;
    }
ris = sqrt( ris ) ;
rjs = sqrt( rjs ) ;
if( fabs( ris - rjs ) > gsym->TolerancePrimary ){
    stat->StatEarly++ ;
    if( gsym->verbose > 0 ) printf( "    Support can't actually define a rotation axis\n" ) ;
    return NULL ;
    }
axis                 = alloc_symmetry_element(gsym) ;
axis->transform_atom = rotate_atom ;
axis->order          = 2 ;
axis->nparam         = 7 ;
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += gsym->CenterOfSomething[k]*gsym->CenterOfSomething[k] ;
r = sqrt(r) ;
if( r > 0 ){
    for( k = 0 ; k < DIMENSION ; k++ )
        axis->normal[k] = gsym->CenterOfSomething[k]/r ;
    }
else {
    axis->normal[0] = 1 ;
    for( k = 1 ; k < DIMENSION ; k++ )
        axis->normal[k] = 0 ;
    }
axis->distance = r ;
for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
    center[k] = ( gsym->Atoms[i].x[k] + gsym->Atoms[j].x[k] ) / 2 - support[k] ;
    r        += center[k]*center[k] ;
    }
r = sqrt(r) ;
if( r <= gsym->TolerancePrimary ){ /* c2 is underdefined, let's do something special */
    if( gsym->MolecularPlane != NULL ){
        if( gsym->verbose > 0 ) printf( "    c2 is underdefined, but there is a molecular plane\n" ) ;
        for( k = 0 ; k < DIMENSION ; k++ )
            axis->direction[k] = gsym->MolecularPlane->normal[k] ;
        }
    else {
        if( gsym->verbose > 0 ) printf( "    c2 is underdefined, trying random direction\n" ) ;
        for( k = 0 ; k < DIMENSION ; k++ )
            center[k] = gsym->Atoms[i].x[k] - gsym->Atoms[j].x[k] ;
        if( fabs( center[2] ) + fabs( center[1] ) > gsym->ToleranceSame ){
            axis->direction[0] =  0 ;
            axis->direction[1] =  center[2] ;
            axis->direction[2] = -center[1] ;
            }
        else {
            axis->direction[0] = -center[2] ;
            axis->direction[1] =  0 ;
            axis->direction[2] =  center[0] ;
            }
        for( k = 0, r = 0 ; k < DIMENSION ; k++ )
            r += axis->direction[k] * axis->direction[k] ;
        r = sqrt(r) ;
        for( k = 0 ; k < DIMENSION ; k++ )
            axis->direction[k] /= r ;
        }
    }
else { /* direction is Ok, renormalize it */
    for( k = 0 ; k < DIMENSION ; k++ )
        axis->direction[k] = center[k]/r ;
    }
if( refine_symmetry_element( axis, 1, gsym, stat ) < 0 ){
    if( gsym->verbose > 0 ) printf( "    refinement failed for the c2 axis\n" ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}

SYMMETRY_ELEMENT *
init_axis_parameters( double a[3], double b[3], double c[3], SYMMETRIES *gsym, STATISTICS *stat )
{
        SYMMETRY_ELEMENT * axis ;
        int                i, order, sign ;
        double             ra, rb, rc, rab, rbc, rac, r ;
        double             angle ;

ra = rb = rc = rab = rbc = rac = 0 ;
for( i = 0 ; i < DIMENSION ; i++ ){
    ra  += a[i]*a[i] ;
    rb  += b[i]*b[i] ;
    rc  += c[i]*c[i] ;
    }
ra = sqrt(ra) ; rb  = sqrt(rb) ; rc  = sqrt(rc) ;
if( fabs( ra - rb ) > gsym->TolerancePrimary || fabs( ra - rc ) > gsym->TolerancePrimary || fabs( rb - rc ) > gsym->TolerancePrimary ){
    stat->StatEarly++ ;
    if( gsym->verbose > 0 ) printf( "    points are not on a sphere\n" ) ;
    return NULL ;
    }
for( i = 0 ; i < DIMENSION ; i++ ){
    rab += (a[i]-b[i])*(a[i]-b[i]) ;
    rac += (a[i]-c[i])*(a[i]-c[i]) ;
    rbc += (c[i]-b[i])*(c[i]-b[i]) ;
    }
rab = sqrt(rab) ;
rac = sqrt(rac) ;
rbc = sqrt(rbc) ;
if( fabs( rab - rbc ) > gsym->TolerancePrimary ){
    stat->StatEarly++ ;
    if( gsym->verbose > 0 ) printf( "    points can't be rotation-equivalent\n" ) ;
    return NULL ;
    }
if( rab <= gsym->ToleranceSame || rbc <= gsym->ToleranceSame || rac <= gsym->ToleranceSame ){
    stat->StatEarly++ ;
    if( gsym->verbose > 0 ) printf( "    rotation is underdefined by these points\n" ) ;
    return NULL ;
    }
rab   = (rab+rbc)/2 ;
angle = M_PI - 2*asin( rac/(2*rab) ) ;
if( gsym->verbose > 1 ) printf( "    rotation angle is %f\n", angle ) ;
if( fabs(angle) <= M_PI/(gsym->MaxAxisOrder+1) ){
    stat->StatEarly++ ;
    if( gsym->verbose > 0 ) printf( "    atoms are too close to a straight line\n" ) ;
    return NULL ;
    }
order = floor( (2*M_PI)/angle + 0.5 ) ;
if( order <= 2 || order > gsym->MaxAxisOrder ){
    stat->StatEarly++ ;
    if( gsym->verbose > 0 ) printf( "    rotation axis order (%d) is not from 3 to %d\n", order, gsym->MaxAxisOrder ) ;
    return NULL ;
    }
axis = alloc_symmetry_element(gsym) ;
axis->order          = order ;
axis->nparam         = 7 ;
for( i = 0, r = 0 ; i < DIMENSION ; i++ )
    r += gsym->CenterOfSomething[i]*gsym->CenterOfSomething[i] ;
r = sqrt(r) ;
if( r > 0 ){
    for( i = 0 ; i < DIMENSION ; i++ )
        axis->normal[i] = gsym->CenterOfSomething[i]/r ;
    }
else {
    axis->normal[0] = 1 ;
    for( i = 1 ; i < DIMENSION ; i++ )
        axis->normal[i] = 0 ;
    }
axis->distance = r ;
axis->direction[0] = (b[1]-a[1])*(c[2]-b[2]) - (b[2]-a[2])*(c[1]-b[1]) ;
axis->direction[1] = (b[2]-a[2])*(c[0]-b[0]) - (b[0]-a[0])*(c[2]-b[2]) ;
axis->direction[2] = (b[0]-a[0])*(c[1]-b[1]) - (b[1]-a[1])*(c[0]-b[0]) ;
/*
 *  Arbitrarily select axis direction so that first non-zero component
 *  or the direction is positive.
 */
sign = 0 ;
if( axis->direction[0] <= 0 )
    if( axis->direction[0] < 0 )
         sign = 1 ;
    else if( axis->direction[1] <= 0 )
             if( axis->direction[1] < 0 )
                  sign = 1 ;
             else if( axis->direction[2] < 0 )
                      sign = 1 ;
if( sign )
    for( i = 0 ; i < DIMENSION ; i++ )
        axis->direction[i] = -axis->direction[i] ;
for( i = 0, r = 0 ; i < DIMENSION ; i++ )
    r += axis->direction[i]*axis->direction[i] ;
r = sqrt(r) ;
for( i = 0 ; i < DIMENSION ; i++ )
    axis->direction[i] /= r ;
if( gsym->verbose > 1 ){
    printf( "    axis origin is at (%g,%g,%g)\n", 
        axis->normal[0]*axis->distance, axis->normal[1]*axis->distance, axis->normal[2]*axis->distance ) ;
    printf( "    axis is in the direction (%g,%g,%g)\n", axis->direction[0], axis->direction[1], axis->direction[2] ) ;
    }
return axis ;
}

SYMMETRY_ELEMENT *
init_higher_axis( int ia, int ib, int ic, SYMMETRIES *gsym, STATISTICS *stat )
{
        SYMMETRY_ELEMENT * axis ;
        double             a[ DIMENSION ], b[ DIMENSION ], c[ DIMENSION ] ;
        int                i ;

if( gsym->verbose > 0 ) printf( "Trying cn axis for the triplet (%d,%d,%d)\n", ia, ib, ic ) ;
stat->StatTotal++ ;
/* Do a quick check of geometry validity */
for( i = 0 ; i < DIMENSION ; i++ ){
    a[i] = gsym->Atoms[ia].x[i] - gsym->CenterOfSomething[i] ;
    b[i] = gsym->Atoms[ib].x[i] - gsym->CenterOfSomething[i] ;
    c[i] = gsym->Atoms[ic].x[i] - gsym->CenterOfSomething[i] ;
    }
if( ( axis = init_axis_parameters( a, b, c, gsym, stat ) ) == NULL ){
    if( gsym->verbose > 0 ) printf( "    no coherent axis is defined by the points\n" ) ;
    return NULL ;
    }
axis->transform_atom = rotate_atom ;
if( refine_symmetry_element( axis, 1, gsym, stat ) < 0 ){
    if( gsym->verbose > 0 ) printf( "    refinement failed for the c%d axis\n", axis->order ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}

/*
 *   Improper axes-specific routines.
 *   These are obtained by slight modifications of normal rotation
 *       routines.
 */
void
rotate_reflect_atom( SYMMETRY_ELEMENT *axis, ATOM *from, ATOM *to )
{
        double             x[3], y[3], a[3], b[3], c[3] ;
        double             angle = 2*M_PI/axis->order ;
        double             a_sin = sin( angle ) ;
        double             a_cos = cos( angle ) ;
        double             dot ;
        int                i ;

if( DIMENSION != 3 ){
    fprintf( stderr, "Catastrophe in rotate_reflect_atom!\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 0 ; i < 3 ; i++ )
    x[i] = from->x[i] - axis->distance * axis->normal[i] ;
for( i = 0, dot = 0 ; i < 3 ; i++ )
    dot += x[i] * axis->direction[i] ;
for( i = 0 ; i < 3 ; i++ )
    a[i] = axis->direction[i] * dot ;
for( i = 0 ; i < 3 ; i++ )
    b[i] = x[i] - a[i] ;
c[0] = b[1]*axis->direction[2] - b[2]*axis->direction[1] ;
c[1] = b[2]*axis->direction[0] - b[0]*axis->direction[2] ;
c[2] = b[0]*axis->direction[1] - b[1]*axis->direction[0] ;
for( i = 0 ; i < 3 ; i++ )
    y[i] = -a[i] + b[i]*a_cos + c[i]*a_sin ;
for( i = 0 ; i < 3 ; i++ )
    to->x[i] = y[i] + axis->distance * axis->normal[i] ;
to->type = from->type ;
}

SYMMETRY_ELEMENT *
init_improper_axis( int ia, int ib, int ic, SYMMETRIES *gsym, STATISTICS *stat )
{
        SYMMETRY_ELEMENT * axis ;
        double             a[ DIMENSION ], b[ DIMENSION ], c[ DIMENSION ] ;
        double             centerpoint[ DIMENSION ] ;
        double             r ;
        int                i ;

if( gsym->verbose > 0 ) printf( "Trying sn axis for the triplet (%d,%d,%d)\n", ia, ib, ic ) ;
stat->StatTotal++ ;
/* First, reduce the problem to Cn case */
for( i = 0 ; i < DIMENSION ; i++ ){
    a[i] = gsym->Atoms[ia].x[i] - gsym->CenterOfSomething[i] ;
    b[i] = gsym->Atoms[ib].x[i] - gsym->CenterOfSomething[i] ;
    c[i] = gsym->Atoms[ic].x[i] - gsym->CenterOfSomething[i] ;
    }
for( i = 0, r = 0 ; i < DIMENSION ; i++ ){
    centerpoint[i] = a[i] + c[i] + 2*b[i] ;
    r             += centerpoint[i]*centerpoint[i] ;
    }
r = sqrt(r) ;
if( r <= gsym->ToleranceSame ){
    stat->StatEarly++ ;
    if( gsym->verbose > 0 ) printf( "    atoms can not define improper axis of the order more than 2\n" ) ;
    return NULL ;
    }
for( i = 0 ; i < DIMENSION ; i++ )
    centerpoint[i] /= r ;
for( i = 0, r = 0 ; i < DIMENSION ; i++ )
    r += centerpoint[i] * b[i] ;
for( i = 0 ; i < DIMENSION ; i++ )
    b[i] = 2*r*centerpoint[i] - b[i] ;
/* Do a quick check of geometry validity */
if( ( axis = init_axis_parameters( a, b, c, gsym, stat ) ) == NULL ){
    if( gsym->verbose > 0 ) printf( "    no coherrent improper axis is defined by the points\n" ) ;
    return NULL ;
    }
axis->transform_atom = rotate_reflect_atom ;
if( refine_symmetry_element( axis, 1, gsym, stat ) < 0 ){
    if( gsym->verbose > 0 ) printf( "    refinement failed for the s%d axis\n", axis->order ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}

/*
 *   Control routines
 */

void
find_center_of_something( SYMMETRIES *gsym )
{
        int                i, j ;
        double             coord_sum[ DIMENSION ] ;
        double             r ;

for( j = 0 ; j < DIMENSION ; j++ )
    coord_sum[j] = 0 ;
for( i = 0 ; i < gsym->AtomsCount ; i++ ){
    for( j = 0 ; j < DIMENSION ; j++ )
        coord_sum[j] += gsym->Atoms[i].x[j] ;
    }
for( j = 0 ; j < DIMENSION ; j++ )
    gsym->CenterOfSomething[j] = coord_sum[j]/gsym->AtomsCount ;
if( gsym->verbose > 0 )
    printf( "Center of something is at %15.10f, %15.10f, %15.10f\n", 
            gsym->CenterOfSomething[0], gsym->CenterOfSomething[1], gsym->CenterOfSomething[2] ) ;
gsym->DistanceFromCenter = (double *) calloc( gsym->AtomsCount, sizeof( double ) ) ;
if( gsym->DistanceFromCenter == NULL ){
    fprintf( stderr, "Unable to allocate array for the distances\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 0 ; i < gsym->AtomsCount ; i++ ){
    for( j = 0, r = 0 ; j < DIMENSION ; j++ )
        r += pow2( gsym->Atoms[i].x[j] - gsym->CenterOfSomething[j] ) ;
    gsym->DistanceFromCenter[i] = r ;
    }
}

void
find_planes(SYMMETRIES *gsym, STATISTICS *stat)
{
        int                i, j ;
        SYMMETRY_ELEMENT * plane ;

plane = init_ultimate_plane(gsym, stat) ;
if( plane != NULL ){
    gsym->MolecularPlane = plane ;
    gsym->PlanesCount++ ;
    gsym->Planes = (SYMMETRY_ELEMENT **) realloc( gsym->Planes, sizeof( SYMMETRY_ELEMENT* ) * gsym->PlanesCount ) ;
    if( gsym->Planes == NULL ){
        perror( "Out of memory in find_planes" ) ;
        exit( EXIT_FAILURE ) ;
        }
    gsym->Planes[ gsym->PlanesCount - 1 ] = plane ;
    }
for( i = 1 ; i < gsym->AtomsCount ; i++ ){
    for( j = 0 ; j < i ; j++ ){
        if( gsym->Atoms[i].type != gsym->Atoms[j].type )
            continue ;
        if( ( plane = init_mirror_plane( i, j, gsym, stat ) ) != NULL ){
            gsym->PlanesCount++ ;
            gsym->Planes = (SYMMETRY_ELEMENT **) realloc( gsym->Planes, sizeof( SYMMETRY_ELEMENT* ) * gsym->PlanesCount ) ;
            if( gsym->Planes == NULL ){
                perror( "Out of memory in find_planes" ) ;
                exit( EXIT_FAILURE ) ;
                }
            gsym->Planes[ gsym->PlanesCount - 1 ] = plane ;
            }
        }
    }
}

void
find_inversion_centers( SYMMETRIES *gsym, STATISTICS *stat )
{
        SYMMETRY_ELEMENT * center ;

if( ( center = init_inversion_center(gsym, stat) ) != NULL ){
    gsym->InversionCenters = (SYMMETRY_ELEMENT **) calloc( 1, sizeof( SYMMETRY_ELEMENT* ) ) ;
    gsym->InversionCenters[0]   = center ;
    gsym->InversionCentersCount = 1 ;
    }
}

void
find_infinity_axis( SYMMETRIES *gsym, STATISTICS *stat )
{
        SYMMETRY_ELEMENT * axis ;

if( ( axis = init_ultimate_axis(gsym,stat) ) != NULL ){
    gsym->NormalAxesCount++ ;
    gsym->NormalAxes = (SYMMETRY_ELEMENT **) realloc( gsym->NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * gsym->NormalAxesCount ) ;
    if( gsym->NormalAxes == NULL ){
        perror( "Out of memory in find_infinity_axes()" ) ;
        exit( EXIT_FAILURE ) ;
        }
    gsym->NormalAxes[ gsym->NormalAxesCount - 1 ] = axis ;
    }
}

void
find_c2_axes( SYMMETRIES *gsym, STATISTICS *stat )
{
        int                i, j, k, l, m ;
        double             center[ DIMENSION ] ;
        double *           distances = calloc( gsym->AtomsCount, sizeof( double ) ) ;
        double             r ;
        SYMMETRY_ELEMENT * axis ;

if( distances == NULL ){
    fprintf( stderr, "Out of memory in find_c2_axes()\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 1 ; i < gsym->AtomsCount ; i++ ){
    for( j = 0 ; j < i ; j++ ){
        if( gsym->Atoms[i].type != gsym->Atoms[j].type )
            continue ;
        if( fabs( gsym->DistanceFromCenter[i] - gsym->DistanceFromCenter[j] ) > gsym->TolerancePrimary )
            continue ; /* A very cheap, but quite effective check */
        /*
         *   First, let's try to get it cheap and use gsym->CenterOfSomething
         */
        for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
            center[k] = ( gsym->Atoms[i].x[k] + gsym->Atoms[j].x[k] ) / 2 ;
            r        += pow2( center[k] - gsym->CenterOfSomething[k] ) ;
            }
        r = sqrt(r) ;
        if( r > 5*gsym->TolerancePrimary ){ /* It's Ok to use gsym->CenterOfSomething */
            if( ( axis = init_c2_axis( i, j, gsym->CenterOfSomething, gsym, stat ) ) != NULL ){
                gsym->NormalAxesCount++ ;
                gsym->NormalAxes = (SYMMETRY_ELEMENT **) realloc( gsym->NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * gsym->NormalAxesCount ) ;
                if( gsym->NormalAxes == NULL ){
                    perror( "Out of memory in find_c2_axes" ) ;
                    exit( EXIT_FAILURE ) ;
                    }
                gsym->NormalAxes[ gsym->NormalAxesCount - 1 ] = axis ;
                }
            continue ;
            }
        /*
         *  Now, C2 axis can either pass through an atom, or through the
         *  middle of the other pair.
         */
        for( k = 0 ; k < gsym->AtomsCount ; k++ ){
            if( ( axis = init_c2_axis( i, j, gsym->Atoms[k].x, gsym, stat ) ) != NULL ){
                gsym->NormalAxesCount++ ;
                gsym->NormalAxes = (SYMMETRY_ELEMENT **) realloc( gsym->NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * gsym->NormalAxesCount ) ;
                if( gsym->NormalAxes == NULL ){
                    perror( "Out of memory in find_c2_axes" ) ;
                    exit( EXIT_FAILURE ) ;
                    }
                gsym->NormalAxes[ gsym->NormalAxesCount - 1 ] = axis ;
                }
            }
        /*
         *  Prepare data for an additional pre-screening check
         */
        for( k = 0 ; k < gsym->AtomsCount ; k++ ){
            for( l = 0, r = 0 ; l < DIMENSION ; l++ )
                r += pow2( gsym->Atoms[k].x[l] - center[l] ) ;
            distances[k] = sqrt(r) ;
            }
        for( k = 0 ; k < gsym->AtomsCount ; k++ ){
            for( l = 0 ; l < gsym->AtomsCount ; l++ ){
                if( gsym->Atoms[k].type != gsym->Atoms[l].type )
                    continue ;
                if( fabs( gsym->DistanceFromCenter[k] - gsym->DistanceFromCenter[l] ) > gsym->TolerancePrimary ||
                    fabs( distances[k] - distances[l] ) > gsym->TolerancePrimary )
                        continue ; /* We really need this one to run reasonably fast! */
                for( m = 0 ; m < DIMENSION ; m++ )
                    center[m] = ( gsym->Atoms[k].x[m] + gsym->Atoms[l].x[m] ) / 2 ;
                if( ( axis = init_c2_axis( i, j, center, gsym, stat ) ) != NULL ){
                    gsym->NormalAxesCount++ ;
                    gsym->NormalAxes = (SYMMETRY_ELEMENT **) realloc( gsym->NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * gsym->NormalAxesCount ) ;
                    if( gsym->NormalAxes == NULL ){
                        perror( "Out of memory in find_c2_axes" ) ;
                        exit( EXIT_FAILURE ) ;
                        }
                    gsym->NormalAxes[ gsym->NormalAxesCount - 1 ] = axis ;
                    }
                }
            }
        }
    }
free( distances ) ;
}

void
find_higher_axes( SYMMETRIES *gsym, STATISTICS *stat )
{
        int                i, j, k ;
        SYMMETRY_ELEMENT * axis ;

for( i = 0 ; i < gsym->AtomsCount ; i++ ){
    for( j = i + 1 ; j < gsym->AtomsCount ; j++ ){
        if( gsym->Atoms[i].type != gsym->Atoms[j].type )
            continue ;
        if( fabs( gsym->DistanceFromCenter[i] - gsym->DistanceFromCenter[j] ) > gsym->TolerancePrimary )
            continue ; /* A very cheap, but quite effective check */
        for( k = 0 ; k < gsym->AtomsCount ; k++ ){
            if( gsym->Atoms[i].type != gsym->Atoms[k].type )
                continue ;
            if( ( fabs( gsym->DistanceFromCenter[i] - gsym->DistanceFromCenter[k] ) > gsym->TolerancePrimary ) ||
                ( fabs( gsym->DistanceFromCenter[j] - gsym->DistanceFromCenter[k] ) > gsym->TolerancePrimary ) )
                    continue ;
            if( ( axis = init_higher_axis( i, j, k, gsym, stat ) ) != NULL ){
                gsym->NormalAxesCount++ ;
                gsym->NormalAxes = (SYMMETRY_ELEMENT **) realloc( gsym->NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * gsym->NormalAxesCount ) ;
                if( gsym->NormalAxes == NULL ){
                    perror( "Out of memory in find_higher_axes" ) ;
                    exit( EXIT_FAILURE ) ;
                    }
                gsym->NormalAxes[ gsym->NormalAxesCount - 1 ] = axis ;
                }
            }
        }
    }
}

void
find_improper_axes(SYMMETRIES *gsym, STATISTICS *stat)
{
        int                i, j, k ;
        SYMMETRY_ELEMENT * axis ;

for( i = 0 ; i < gsym->AtomsCount ; i++ ){
    for( j = i + 1 ; j < gsym->AtomsCount ; j++ ){
        for( k = 0 ; k < gsym->AtomsCount ; k++ ){
            if( ( axis = init_improper_axis( i, j, k, gsym, stat ) ) != NULL ){
                gsym->ImproperAxesCount++ ;
                gsym->ImproperAxes = (SYMMETRY_ELEMENT **) realloc( gsym->ImproperAxes, sizeof( SYMMETRY_ELEMENT* ) * gsym->ImproperAxesCount ) ;
                if( gsym->ImproperAxes == NULL ){
                    perror( "Out of memory in find_higher_axes" ) ;
                    exit( EXIT_FAILURE ) ;
                    }
                gsym->ImproperAxes[ gsym->ImproperAxesCount - 1 ] = axis ;
                }
            }
        }
    }
}

void
report_planes( SYMMETRIES *gsym )
{
        int           i ;

if( gsym->PlanesCount == 0 )
    printf( "There are no planes of symmetry in the molecule\n" ) ;
else {
    if( gsym->PlanesCount == 1 )
         printf( "There is a plane of symmetry in the molecule\n" ) ;
    else printf( "There are %d planes of symmetry in the molecule\n", gsym->PlanesCount ) ;
    printf( "     Residual          Direction of the normal           Distance\n" ) ;
    for( i = 0 ; i < gsym->PlanesCount ; i++ ){
        printf( "%3d %8.4e ", i, gsym->Planes[i]->maxdev ) ;
        printf( "(%11.8f,%11.8f,%11.8f) ", gsym->Planes[i]->normal[0], gsym->Planes[i]->normal[1], gsym->Planes[i]->normal[2] ) ;
        printf( "%14.8f\n", gsym->Planes[i]->distance ) ;
        }
    }
}

void
report_inversion_centers( SYMMETRIES *gsym )
{
if( gsym->InversionCentersCount == 0 )
     printf( "There is no inversion center in the molecule\n" ) ;
else {
    printf( "There is an inversion center in the molecule\n" ) ;
    printf( "     Residual                      Position\n" ) ;
    printf( "   %8.4e ", gsym->InversionCenters[0]->maxdev ) ;
    printf( "(%14.8f,%14.8f,%14.8f)\n",
        gsym->InversionCenters[0]->distance * gsym->InversionCenters[0]->normal[0],
        gsym->InversionCenters[0]->distance * gsym->InversionCenters[0]->normal[1],
        gsym->InversionCenters[0]->distance * gsym->InversionCenters[0]->normal[2] ) ;
    }
}

void
report_axes( SYMMETRIES *gsym )
{
        int           i ;

if( gsym->NormalAxesCount == 0 )
    printf( "There are no normal axes in the molecule\n" ) ;
else {
    if( gsym->NormalAxesCount == 1 )
         printf( "There is a normal axis in the molecule\n" ) ;
    else printf( "There are %d normal axes in the molecule\n", gsym->NormalAxesCount ) ;
    printf( "     Residual  Order         Direction of the axis                         Supporting point\n" ) ;
    for( i = 0 ; i < gsym->NormalAxesCount ; i++ ){
        printf( "%3d %8.4e ", i, gsym->NormalAxes[i]->maxdev ) ;
        if( gsym->NormalAxes[i]->order == 0 )
             printf( "Inf " ) ;
        else printf( "%3d ", gsym->NormalAxes[i]->order ) ;
        printf( "(%11.8f,%11.8f,%11.8f) ", 
            gsym->NormalAxes[i]->direction[0], gsym->NormalAxes[i]->direction[1], gsym->NormalAxes[i]->direction[2] ) ;
        printf( "(%14.8f,%14.8f,%14.8f)\n", 
            gsym->NormalAxes[0]->distance * gsym->NormalAxes[0]->normal[0],
            gsym->NormalAxes[0]->distance * gsym->NormalAxes[0]->normal[1],
            gsym->NormalAxes[0]->distance * gsym->NormalAxes[0]->normal[2] ) ;
        }
    }
}

void
report_improper_axes( SYMMETRIES *gsym )
{
        int           i ;

if( gsym->ImproperAxesCount == 0 )
    printf( "There are no improper axes in the molecule\n" ) ;
else {
    if( gsym->ImproperAxesCount == 1 )
         printf( "There is an improper axis in the molecule\n" ) ;
    else printf( "There are %d improper axes in the molecule\n", gsym->ImproperAxesCount ) ;
    printf( "     Residual  Order         Direction of the axis                         Supporting point\n" ) ;
    for( i = 0 ; i < gsym->ImproperAxesCount ; i++ ){
        printf( "%3d %8.4e ", i, gsym->ImproperAxes[i]->maxdev ) ;
        if( gsym->ImproperAxes[i]->order == 0 )
             printf( "Inf " ) ;
        else printf( "%3d ", gsym->ImproperAxes[i]->order ) ;
        printf( "(%11.8f,%11.8f,%11.8f) ", 
            gsym->ImproperAxes[i]->direction[0], gsym->ImproperAxes[i]->direction[1], gsym->ImproperAxes[i]->direction[2] ) ;
        printf( "(%14.8f,%14.8f,%14.8f)\n", 
            gsym->ImproperAxes[0]->distance * gsym->ImproperAxes[0]->normal[0],
            gsym->ImproperAxes[0]->distance * gsym->ImproperAxes[0]->normal[1],
            gsym->ImproperAxes[0]->distance * gsym->ImproperAxes[0]->normal[2] ) ;
        }
    }
}

/*
 *  General symmetry handling
 */
void
report_and_reset_counters( STATISTICS *stat )
{
printf( "  %10ld candidates examined\n"
        "  %10ld removed early\n"
        "  %10ld removed during initial mating stage\n"
        "  %10ld removed as duplicates\n"
        "  %10ld removed because of the wrong transformation order\n"
        "  %10ld removed after unsuccessful optimization\n"
        "  %10ld accepted\n",
    stat->StatTotal, stat->StatEarly, stat->StatPairs, stat->StatDups, stat->StatOrder, stat->StatOpt, stat->StatAccept ) ;
    init_statistics(stat) ;
}

void
find_symmetry_elements( 
                SYMMETRIES *gsym,
                STATISTICS *stat )
{
find_center_of_something(gsym) ;
if( gsym->verbose > -1 ){
    printf( "Looking for the inversion center\n" ) ;
    }
find_inversion_centers(gsym, stat) ;
if( gsym->verbose > -1 ){
    report_and_reset_counters(stat) ;
    printf( "Looking for the planes of symmetry\n" ) ;
    }
find_planes(gsym, stat) ;
if( gsym->verbose > -1 ){
    report_and_reset_counters(stat) ;
    printf( "Looking for infinity axis\n" ) ;
    }
find_infinity_axis(gsym, stat) ;
if( gsym->verbose > -1 ){
    report_and_reset_counters(stat) ;
    printf( "Looking for C2 axes\n" ) ;
    }
find_c2_axes(gsym, stat) ;
if( gsym->verbose > -1 ){
    report_and_reset_counters(stat) ;
    printf( "Looking for higher axes\n" ) ;
    }
find_higher_axes(gsym, stat) ;
if( gsym->verbose > -1 ){
    report_and_reset_counters(stat) ;
    printf( "Looking for the improper axes\n" ) ;
    }
find_improper_axes(gsym, stat) ;
if( gsym->verbose > -1 ){
    report_and_reset_counters(stat) ;
    }
}

int
compare_axes( const void *a, const void *b )
{
        SYMMETRY_ELEMENT * axis_a = *(SYMMETRY_ELEMENT**) a ;
        SYMMETRY_ELEMENT * axis_b = *(SYMMETRY_ELEMENT**) b ;
        int                i, order_a, order_b ;

order_a = axis_a->order ; if( order_a == 0 ) order_a = 10000 ;
order_b = axis_b->order ; if( order_b == 0 ) order_b = 10000 ;
if( ( i = order_b - order_a ) != 0 ) return i ;
if( axis_a->maxdev > axis_b->maxdev ) return -1 ;
if( axis_a->maxdev < axis_b->maxdev ) return  1 ;
return 0 ;
}

void
sort_symmetry_elements( SYMMETRIES *gsym )
{
if( gsym->PlanesCount > 1 ){
    qsort( gsym->Planes, gsym->PlanesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
    }
if( gsym->NormalAxesCount > 1 ){
    qsort( gsym->NormalAxes, gsym->NormalAxesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
    }
if( gsym->ImproperAxesCount > 1 ){
    qsort( gsym->ImproperAxes, gsym->ImproperAxesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
    }
}

void
report_symmetry_elements_verbose( SYMMETRIES *gsym )
{
report_inversion_centers(gsym) ;
report_axes(gsym) ;
report_improper_axes(gsym) ;
report_planes(gsym) ;
}

void
summarize_symmetry_elements( SYMMETRIES *gsym )
{
        int          i ;

gsym->NormalAxesCounts   = (int*) calloc( gsym->MaxAxisOrder+1, sizeof( int ) ) ;
gsym->ImproperAxesCounts = (int*) calloc( gsym->MaxAxisOrder+1, sizeof( int ) ) ;
for( i = 0 ; i < gsym->NormalAxesCount ; i++ )
    gsym->NormalAxesCounts[ gsym->NormalAxes[i]->order ]++ ;
for( i = 0 ; i < gsym->ImproperAxesCount ; i++ )
    gsym->ImproperAxesCounts[ gsym->ImproperAxes[i]->order ]++ ;
}

void
report_symmetry_elements_brief( SYMMETRIES *gsym )
{
        int          i ;
        char *       symmetry_code = calloc( 1, 10*(gsym->PlanesCount+gsym->NormalAxesCount+gsym->ImproperAxesCount+gsym->InversionCentersCount+2) ) ;
        char         buf[ 100 ] ;

if( symmetry_code == NULL ){
    fprintf( stderr, "Unable to allocate memory for symmetry ID code in report_symmetry_elements_brief()\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
if( gsym->PlanesCount + gsym->NormalAxesCount + gsym->ImproperAxesCount + gsym->InversionCentersCount == 0 )
    printf( "Molecule has no symmetry elements\n" ) ;
else {
    printf( "Molecule has the following symmetry elements: " ) ;
    if( gsym->InversionCentersCount > 0 ) strcat( symmetry_code, "(i) " ) ;
    if( gsym->NormalAxesCounts[0] == 1 )
         strcat( symmetry_code, "(Cinf) " ) ;
    if( gsym->NormalAxesCounts[0] >  1 ) {
        sprintf( buf, "%d*(Cinf) ", gsym->NormalAxesCounts[0] ) ;
        strcat( symmetry_code, buf ) ;
        }
    for( i = gsym->MaxAxisOrder ; i >= 2 ; i-- ){
        if( gsym->NormalAxesCounts[i] == 1 ){ sprintf( buf, "(C%d) ", i ) ; strcat( symmetry_code, buf ) ; }
        if( gsym->NormalAxesCounts[i] >  1 ){ sprintf( buf, "%d*(C%d) ", gsym->NormalAxesCounts[i], i ) ; strcat( symmetry_code, buf ) ; }
        }
    for( i = gsym->MaxAxisOrder ; i >= 2 ; i-- ){
        if( gsym->ImproperAxesCounts[i] == 1 ){ sprintf( buf, "(S%d) ", i ) ; strcat( symmetry_code, buf ) ; }
        if( gsym->ImproperAxesCounts[i] >  1 ){ sprintf( buf, "%d*(S%d) ", gsym->ImproperAxesCounts[i], i ) ; strcat( symmetry_code, buf ) ; }
        }
    if( gsym->PlanesCount == 1 ) strcat( symmetry_code, "(sigma) " ) ;
    if( gsym->PlanesCount >  1 ){ sprintf( buf, "%d*(sigma) ", gsym->PlanesCount ) ; strcat( symmetry_code, buf ) ; }
    printf( "%s\n", symmetry_code ) ;
    }
gsym->SymmetryCode = symmetry_code ;
}

void
identify_point_group( SYMMETRIES *gsym )
{
        int            i ;
        int            last_matching = -1 ;
        int            matching_count = 0 ;

for( i = 0 ; i < PointGroupsCount ; i++ ){
    if( strcmp( gsym->SymmetryCode, PointGroups[i].symmetry_code ) == 0 ){
        if( PointGroups[i].check() == 1 ){
            last_matching = i ;
            matching_count++ ;
            }
        else {
            if( gsym->verbose > -2 ){
                printf( "It looks very much like %s, but it is not since %s\n", 
                    PointGroups[i].group_name, PointGroupRejectionReason ) ;
                }
            }
        }
    }
if( matching_count == 0 ){
    printf( "These symmetry elements match no point group I know of. Sorry.\n" ) ;
    }
if( matching_count >  1 ){
    printf( "These symmetry elements match more than one group I know of.\n"
            "SOMETHING IS VERY WRONG\n" ) ;
    printf( "Matching groups are:\n" ) ;
    for( i = 0 ; i < PointGroupsCount ; i++ ){
        if( ( strcmp( gsym->SymmetryCode, PointGroups[i].symmetry_code ) == 0 ) && ( PointGroups[i].check() == 1 ) ){
            printf( "    %s\n", PointGroups[i].group_name ) ;
            }
        }
    }
if( matching_count == 1 ){
    printf( "It seems to be the %s point group\n", PointGroups[last_matching].group_name ) ;
    }
}

/* Input/Output */

int
read_coordinates( FILE *in, SYMMETRIES *gsym )
{
        int i ;

if( fscanf( in, "%d", &gsym->AtomsCount ) != 1 ){
    fprintf( stderr, "Error reading atom count\n" ) ;
    return -1 ;
    }
if( gsym->verbose > 0 ) printf( "gsym->Atoms count = %d\n", gsym->AtomsCount ) ;
gsym->Atoms = calloc( gsym->AtomsCount, sizeof( ATOM ) ) ;
if( gsym->Atoms == NULL ){
    fprintf( stderr, "Out of memory for atoms coordinates\n" ) ;
    return -1 ;
    }
for( i = 0 ; i < gsym->AtomsCount ; i++ ){
    if( fscanf( in, "%d %lg %lg %lg\n", &gsym->Atoms[i].type, &gsym->Atoms[i].x[0], &gsym->Atoms[i].x[1], &gsym->Atoms[i].x[2] ) != 4 ){
        fprintf( stderr, "Error reading description of the atom %d\n", i ) ;
        return -1 ;
        }
    }
return 0 ;
}

int
set_coordinates( 
                SYMMETRIES *gsym, 
                int * nat,
                int * typat,
                double *xat )
{
        int i ;

gsym->AtomsCount = nat[0] ;
gsym->verbose = 1 ;
if( gsym->verbose > 0 ) printf( "gsym->Atoms count = %d\n", gsym->AtomsCount ) ;
gsym->Atoms = calloc( gsym->AtomsCount, sizeof( ATOM ) ) ;
if( gsym->Atoms == NULL ){
    fprintf( stderr, "Out of memory for atoms coordinates\n" ) ;
    return -1 ;
    }
for( i = 0 ; i < gsym->AtomsCount ; i++ ){
    gsym->Atoms[i].type = typat[i] ;
    gsym->Atoms[i].x[0] = xat[i*DIMENSION] ;
    gsym->Atoms[i].x[1] = xat[i*DIMENSION+1] ;
    gsym->Atoms[i].x[2] = xat[i*DIMENSION+2] ;
    printf("%d %5.2f %5.2f %5.2f\n",gsym->Atoms[i].type,gsym->Atoms[i].x[0], gsym->Atoms[i].x[1], gsym->Atoms[i].x[2]) ;
    }
return 0 ;
}


/* The main routine */
int
find_symmetries_( 
                int *nat,
                int * typat,
                double *xat,
                SYMMETRIES * gsym
                )
{
/*        SYMMETRIES    * gsym = calloc( 1, sizeof( SYMMETRIES ) ) ; */
        gsym = (SYMMETRIES *)calloc( 1, sizeof( SYMMETRIES ) ) ;
        STATISTICS    * stat = calloc( 1, sizeof( STATISTICS ) ) ;

/* init gsym and stat */
init_symmetries(gsym);
init_statistics(stat);

set_coordinates(gsym, nat, typat, xat) ;
find_symmetry_elements(gsym, stat) ;
sort_symmetry_elements(gsym) ;
summarize_symmetry_elements(gsym) ;
if( gsym->BadOptimization )
    printf( "Refinement of some symmetry elements was terminated before convergence was reached.\n"
            "Some symmetry elements may remain unidentified.\n" ) ;
if( gsym->verbose >= 0 )
    report_symmetry_elements_verbose(gsym) ;
report_symmetry_elements_brief(gsym) ;
identify_point_group(gsym) ;
/* exit( EXIT_SUCCESS ) ; */
}
