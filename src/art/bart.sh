#!/bin/csh
#____________________________ ATOMS 
setenv NATOMS                           8    # Number of atoms in the problem
setenv type1                            C 
setenv type2                            H 
#_____________________________ ART 
setenv Setup_Initial              .false.   # for fine tune of lanczos parameters
setenv EVENT_TYPE                     NEW   # Either 
                                            # 'NEW',
                                            # 'REFINE_SADDLE' when further converging a saddle point
                                            # 'REFINE_AND_RELAX', to refine at the saddle and check			                                                 #  the final minimum
                                            # 'GUESS_DIRECTION' 
setenv ENERGY_CALC                    BIG
setenv Temperature                   -0.5   # Fictive temperature, if negative always reject the event
setenv Max_Number_Events             1000   # Maximum number of events
setenv Type_of_Events               local   # Activation: global, local, list_local, list
setenv Radius_Initial_Deformation     1.2   # Cutoff for local and list_local (in Ang)
setenv Central_Atom                     1   # Number of the atom around which the initial move takes place
setenv sym_break_dist               0.001   # Breaks the symmetry of the crystal by randomly displacing
                                            # all atoms by this distance
setenv Activation_MaxIter             300   # Maximum number of iteraction for reaching the saddle point
setenv Initial_Step_Size             0.05   # Size of initial displacement, (in Ang)
setenv Increment_Size                0.09   # Overall scale for the increment moves
setenv Force_Threshold_Perp_Rel       0.5   # Threshold for perpendicular relaxation
#_____________________________ HARMONIC WELL
setenv Basin_Factor                   2.1   # Factor multiplying Increment_Size for leaving the basin
setenv Max_Perp_Moves_Basin             3   # Maximum number of perpendicular steps leaving basin
setenv Min_Number_KSteps                3   # Min. number of ksteps before calling lanczos 
setenv Eigenvalue_Threshold         -0.05   # Eigenvalue threshold for leaving basin
setenv Max_Iter_Basin                  12   # Maximum number of iteraction for leaving the basin (kter)
#_____________________________ LANCZOS
setenv gnrm                       1.0E-05   # convergence criterion for the wavefunction optimization for Lanczos 
                                            # Vectors. Reasonable values are in most cases between 2.d-5 and 1.d-5.
setenv Lanczos_of_minimum         .False.   # Calculation of the Hessian for each minimum
setenv Lanczos_SCLoop                   1   # of iterations in the lanczos Self consistent loop
setenv Lanczos_collinear             0.73   # Exit criteria of the lanczos Self consistent loop
setenv Number_Lanczos_Vectors_H        13   # Number of vectors included in lanczos procedure in the Harmonic well
setenv Number_Lanczos_Vectors_C        12   # Number of vectors included in lanczos procedure in convergence
setenv delta_disp_Lanczos            0.01   # Step of the numerical derivative of forces in lanczos (Ang)
#_____________________________ CONVERGENCE
setenv Max_Perp_Moves_Activ             3   # Maximum number of perpendicular steps during activation
setenv Exit_Force_Threshold          0.25   # Threshold for convergence at saddle point
setenv Prefactor_Push_Over_Saddle     0.3   # Fraction of displacement over the saddle
setenv delta_threshold               0.09   # if delta_e < delta_thr .and. delr < delr_thr, in the 
setenv delr_threshold                 0.7   # convergence stage, then we kill the event
#_____________________________ DIIS
setenv Use_DIIS                   .False.   # Use DIIS for the final convergence to saddle
setenv Iterative                  .False.   # Iterative use of Lanczos & DIIS
setenv Inflection                     100   # Number of Lanczos steps after an inflection in the eigenvalue
setenv DIIS_Force_Threshold           0.4   # Force threshold for call DIIS
setenv DIIS_Memory                      5   # Number of vectors kepts in memory for algorithm
setenv DIIS_Step_size                0.02   # prefactor multiplying forces
setenv FACTOR_DIIS                    4.0   # times Increment_Size, max allowed diis step size
setenv MAX_DIIS                       100   # max diis iterations per call
setenv DIIS_Check_Eigenvector      .True.   # Check that the final state is indeed a saddle
#______________________ INPUT/OUTPUT 
setenv FILECOUNTER            filecounter   # File tracking  the file (event) number - facultative
setenv LOGFILE                   log.file   # General output for message
setenv EVENTSLIST             events.list   # list of events with success or failure
setenv Save_Conf_Int               .True.   # Save the configuration at every step?
setenv Write_restart_file          .True.   # It is useful only for ab-initio
setenv RESTART                restart.dat   # current data for restarting event
setenv Write_xyz                   .True.   # Writes min. and sad. configurations in .xyz format.
setenv REFCONFIG                refconfig   # Reference configuration for refine saddle. Without ext.
setenv INITDIR                    initdir   # Configuration for GUESS_direction. Without ext.
###### RUN THE SIMULATION #######################################################################
mpirun -np 10 ./bart > bart.out
