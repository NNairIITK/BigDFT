#!/bin/csh

######  ATOMS ############################################################################################
setenv NATOMS     20                           # Number of atoms in the problem
setenv type1      C

######  ART ##############################################################################################

setenv EVENT_TYPE  NEW     # Either 'NEW', 'REFINE_SADDLE' when further converging a saddle point
                           # Or "REFINE_AND_RELAX", to refine at the saddle
			   # and check the final minimum

setenv Temperature                        10   # Temperature in kcal/mol, if negative always reject the event
setenv Prefactor_Push_Over_Saddle        0.3   # Fraction of displacement over the saddle
setenv Max_Number_Events                  20   # Maximum number of events
setenv Type_of_Events                  local   # Initial move for events - global or local
setenv Radius_Initial_Deformation        2.2   # Cutoff for local-move (in angstroems)
setenv Central_Atom                            # Number of the atom # around which the initial move takes place

setenv Eigenvalue_Threshold             -2.0   # Eigenvalue threshold for leaving basin
setenv Exit_Force_Threshold              1.5   # Threshold for convergence at saddle point

setenv Increment_Size                   0.07   # Overall scale for the increment moves in activation
setenv Initial_Step_Size                0.05   # Size of initial displacement, in A

setenv sym_break_dist                 0.0001   # Breaks the symmetry of the
                                               # crystal by randomly displacing all atoms by this distance
setenv Max_Perp_Moves_Basin                2   # Maximum number of perpendicular steps leaving basin
setenv Min_Number_KSteps                   2   # Min. number of ksteps before calling lanczos 
setenv Lanczos_of_minimum             .True.   # Calculation of the Hessian for each minimum
setenv Basin_Factor                      3.0   # Factor multiplying Increment_Size for leaving the basin

setenv Number_Lanczos_Vectors             16   # Number of vectors included in lanczos procedure
setenv Step_Lanczos                     0.06   # Step of the numerical derivative of forces in lanczos (Ang)
setenv Max_Perp_Moves_Activ                8   # Maximum number of perpendicular steps during activation
setenv Force_Threshold_Perp_Rel          0.5   # Threshold for perpendicular relaxation
setenv Max_Iter_Basin                    200   # Maximum number of iteraction for leaving the basin (kter)
setenv Max_Iter_Activation               100   # Maximum number of iteraction during activation (iter)
setenv Save_Conf_Int                  .True.   # Save the configuration at every step?    

################## Direction inversion in iterative subspace
setenv Use_DIIS                       .True.   # Use DIIS for the final convergence to saddle
setenv DIIS_Force_Threshold              0.4   # Force threshold for convergence
setenv DIIS_Memory                         6   # Number of vectors kepts in memory for algorithm
setenv DIIS_MaxIter                      100   # Maximum number of iteractions for the DIIS scheme
setenv DIIS_Check_Eigenvector         .True.   # Check that the final state is indeed a saddle
setenv DIIS_Step_size                   0.01   # Step size for the position

############### Input ####################################################################################
setenv FILECOUNTER               filecounter   # File tracking  the file (event) number - facultative
setenv REFCONFIG                   refconfig   # Reference configuration (actual local minimum)

############### Output ###################################################################################
setenv LOGFILE                      log.file   # General output for message
setenv EVENTSLIST                events.list   # list of events with success or failure
setenv RESTART                   restart.dat   # current data for restarting event

############### Run the simulation #######################################################################

mpirun bart > bart.out

