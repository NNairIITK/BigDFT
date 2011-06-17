#!/bin/csh
#____________________________ ATOMS 
setenv gnrm          1.0E-05
setenv Inflection  4 
setenv NATOMS     20                         # Number of atoms in the problem
setenv type1      C
#_____________________________ ART 
setenv EVENT_TYPE  NEW     # Either 'NEW', 'REFINE_SADDLE' when further converging a saddle point
                           # Or "REFINE_AND_RELAX", to refine at the saddle
			   # and check the final minimum
setenv Temperature                       0.5   # Temperature in kcal/mol, if negative always reject the event
setenv Max_Number_Events                1000   # Maximum number of events
setenv Type_of_Events                  local   # Initial move for events - global or local
setenv Radius_Initial_Deformation        2.4   # Cutoff for local-move (in angstroems)
setenv Central_Atom                            # Number of the atom around which the initial move takes place
setenv sym_break_dist                  0.001   # Breaks the symmetry of the crystal by randomly displacing
                                               # all atoms by this distance
setenv Activation_MaxIter                300   # Maximum number of iteraction for reaching the saddle point
setenv Initial_Step_Size                0.05   # Size of initial displacement, in A
setenv Increment_Size                   0.08   # Overall scale for the increment moves
setenv Force_Threshold_Perp_Rel          0.5   # Threshold for perpendicular relaxation
#_____________________________ HARMONIC WELL
setenv Basin_Factor                      2.4   # Factor multiplying Increment_Size for leaving the basin
setenv Max_Perp_Moves_Basin                3   # Maximum number of perpendicular steps leaving basin
setenv Min_Number_KSteps                   2   # Min. number of ksteps before calling lanczos 
setenv Eigenvalue_Threshold             -1.5   # Eigenvalue threshold for leaving basin
setenv Max_Iter_Basin                     20   # Maximum number of iteraction for leaving the basin (kter)
#_____________________________ LANCZOS
setenv Lanczos_of_minimum            .True.    # Calculation of the Hessian for each minimum
setenv Number_Lanczos_Vectors             15   # Number of vectors included in lanczos procedure
setenv delta_disp_Lanczos              0.025   # Step of the numerical derivative of forces in lanczos (Ang)
#_____________________________ CONVERGENCE
setenv Max_Perp_Moves_Activ                8   # Maximum number of perpendicular steps during activation
setenv Exit_Force_Threshold             0.25   # Threshold for convergence at saddle point
setenv Prefactor_Push_Over_Saddle        0.1   # Fraction of displacement over the saddle
setenv Save_Conf_Int                  .False.   # Save the configuration at every step?
#_____________________________ DIIS
setenv Iterative                      .False.   # Iterative use of Lanczos & DIIS
setenv Use_DIIS                       .True.   # Use DIIS for the final convergence to saddle
setenv DIIS_Force_Threshold              2.0   # Force threshold for call DIIS
setenv DIIS_Memory                         5   # Number of vectors kepts in memory for algorithm
setenv DIIS_Check_Eigenvector         .True.   # Check that the final state is indeed a saddle
setenv DIIS_Step_size                   0.03   # prefactor multiplying forces
setenv FACTOR_DIIS                       5.0   # times Increment_Size, max allowed diis step size
setenv MAX_DIIS                          150   # max diis iterations per call
#_____________________________ INPUT 
setenv FILECOUNTER               filecounter   # File tracking  the file (event) number - facultative
setenv REFCONFIG                   refconfig   # Reference configuration (actual local minimum)
#_____________________________ OUTPUT 
setenv LOGFILE                      log.file   # General output for message
setenv EVENTSLIST                events.list   # list of events with success or failure
setenv Write_restart_file             .False.   # It is useful only for ab-initio
setenv RESTART                   restart.dat   # current data for restarting event
setenv Write_JMOL                     .False.   # Writes the configuration in jmol format.
###### RUN THE SIMULATION #######################################################################
mpirun ./bart >> bart.out
