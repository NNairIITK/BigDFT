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

#include "symmetry.h"

/* The main program */
int
main( int argc, char **argv )
{
        char          *program = *argv ;
        FILE          *in ;
        SYMMETRIES    * gsym = calloc( 1, sizeof( SYMMETRIES ) ) ; 
        STATISTICS    * stat = calloc( 1, sizeof( STATISTICS ) ) ;

/* init gsym and stat */
init_symmetries(gsym);
init_statistics(stat);

for( argc--, argv++ ; argc > 0 ; argc -= 2, argv += 2 ){
    if( **argv != '-' )
        break ;
    if( strcmp( *argv, "-help"         ) == 0 ||
        strcmp( *argv, "-h"            ) == 0 ||
        strcmp( *argv, "-?"            ) == 0 ){
        argc++ ; argv-- ;
        printf( "%s [option value ...] [filename]\n" 
                "Valid options are:\n"
                "  -verbose      (%3d) Determines verbosity level\n"
                "                      All values above 0 are intended for debugging purposes\n"
                "  -maxaxisorder (%3d) Maximum order of rotation axis to look for\n"
                "  -maxoptcycles (%3d) Maximum allowed number of cycles in symmetry element optimization\n"
                "  --                  Terminates option processing\n"
                "Defaults should be Ok for these:\n"
                "  -same         (%8g) Atoms are colliding if distance falls below this value\n"
                "  -primary      (%8g) Initial loose criterion for atom equivalence\n"
                "  -final        (%8g) Final criterion for atom equivalence\n"
                "  -maxoptstep   (%8g) Largest step allowed in symmetry element optimization\n"
                "  -minoptstep   (%8g) Termination criterion in symmetry element optimization\n"
                "  -gradstep     (%8g) Finite step used in numeric gradient evaluation\n" 
                "  -minchange    (%8g) Minimum allowed change in target function\n" 
                "  -minchgcycles (%8d)  Number of minchange cycles before optimization stops\n",
            program, gsym->verbose, gsym->MaxAxisOrder, gsym->MaxOptCycles, gsym->ToleranceSame, gsym->TolerancePrimary,
            gsym->ToleranceFinal, gsym->MaxOptStep, gsym->MinOptStep, gsym->GradientStep, gsym->OptChangeThreshold, gsym->OptChangeHits ) ;
        printf( "\n"
                "Input is expected in the following format:\n"
                "number_of_atoms\n"
                "AtomicNumber X Y Z\n"
                "...\n" ) ;
        printf( "\n"
                "Note that only primitive rotations will be reported\n" ) ;
        printf( "This is version $Revision: 1.16 $ ($Date: 2003/04/04 13:05:03 $)\n" ) ;
        exit( EXIT_SUCCESS ) ;
        }
    else
    if( strcmp( *argv, "--"            ) == 0 ){
        argc-- ; argv++ ; break ;
        }
    if( argc < 2 ){
        fprintf( stderr, "Missing argument for \"%s\"\n", *argv ) ;
        exit( EXIT_FAILURE ) ;
        }
    if( strcmp( *argv, "-minchgcycles" ) == 0 ){
        if( sscanf( argv[1], "%d", &gsym->OptChangeHits ) != 1 ){
            fprintf( stderr, "Invalid parameter for -minchgcycles: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else
    if( strcmp( *argv, "-minchange"    ) == 0 ){
        if( sscanf( argv[1], "%lg", &gsym->OptChangeThreshold ) != 1 ){
            fprintf( stderr, "Invalid parameter for -minchange: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else
    if( strcmp( *argv, "-same"         ) == 0 ){
        if( sscanf( argv[1], "%lg", &gsym->ToleranceSame ) != 1 ){
            fprintf( stderr, "Invalid parameter for -same: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else
    if( strcmp( *argv, "-primary"      ) == 0 ){
        if( sscanf( argv[1], "%lg", &gsym->TolerancePrimary ) != 1 ){
            fprintf( stderr, "Invalid parameter for -primary: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else
    if( strcmp( *argv, "-final"        ) == 0 ){
        if( sscanf( argv[1], "%lg", &gsym->ToleranceFinal ) != 1 ){
            fprintf( stderr, "Invalid parameter for -final: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else
    if( strcmp( *argv, "-maxoptstep"   ) == 0 ){
        if( sscanf( argv[1], "%lg", &gsym->MaxOptStep ) != 1 ){
            fprintf( stderr, "Invalid parameter for -maxoptstep: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else
    if( strcmp( *argv, "-minoptstep"   ) == 0 ){
        if( sscanf( argv[1], "%lg", &gsym->MinOptStep ) != 1 ){
            fprintf( stderr, "Invalid parameter for -minoptstep: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else
    if( strcmp( *argv, "-gradstep"     ) == 0 ){
        if( sscanf( argv[1], "%lg", &gsym->GradientStep ) != 1 ){
            fprintf( stderr, "Invalid parameter for -gradstep: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else
    if( strcmp( *argv, "-verbose"      ) == 0 ){
        if( sscanf( argv[1], "%d", &gsym->verbose ) != 1 ){
            fprintf( stderr, "Invalid parameter for -verbose: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else
    if( strcmp( *argv, "-maxoptcycles" ) == 0 ){
        if( sscanf( argv[1], "%d", &gsym->MaxOptCycles ) != 1 ){
            fprintf( stderr, "Invalid parameter for -maxoptcycles: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else
    if( strcmp( *argv, "-maxaxisorder" ) == 0 ){
        if( sscanf( argv[1], "%d", &gsym->MaxAxisOrder ) != 1 ){
            fprintf( stderr, "Invalid parameter for -maxaxisorder: \"%s\"\n", argv[1] ) ;
            exit( EXIT_FAILURE ) ;
            }
        }
    else {
        fprintf( stderr, "Unrecognized option \"%s\"\n", *argv ) ;
        exit( EXIT_FAILURE ) ;
        }
    }
if( argc > 0 ){
    if( ( in = fopen( *argv, "rt" ) ) == NULL ){
        perror( *argv ) ;
        exit( EXIT_FAILURE ) ;
        }
    }
else {
    in = stdin ;
    }
if( read_coordinates( in, gsym ) < 0 ){
    fprintf( stderr, "Error reading in atomic coordinates\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
fclose( in ) ;
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
exit( EXIT_SUCCESS ) ;
}
