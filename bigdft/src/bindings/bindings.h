/** @file 
     Header for the BigDFT bindings.
    @author
     Copyright (C) 2011-2015 BigDFT group (DC)
     This file is distributed under the terms of the
     GNU General Public License, see ~/COPYING file
     or http://www.gnu.org/copyleft/gpl.txt .
     For the list of contributors, see ~/AUTHORS 
*/


#ifndef BINDINGS_H
#define BINDINGS_H

#ifdef HAVE_DEBUG
#define DBG_MEM(A, T) {int i__; for (i__ = 0; i__ < sizeof(T) / sizeof(void*); i__++) \
                                  fprintf(stderr, "DBG (%2d) -> %p\n", i__, ((void**)A)[i__]); }
#else
#define DBG_MEM(A, T)
#endif

/* Constructors of C wrappers around already built Fortran objects. */
BigDFT_Dict*    bigdft_dict_new_from_fortran   (f90_dictionary_pointer dict);
BigDFT_Atoms*   bigdft_atoms_new_from_fortran  (f90_atoms_data_pointer at);
BigDFT_Inputs*  bigdft_inputs_new_from_fortran (f90_input_variables_pointer inputs);

/* Additional private methods. */
void _inputs_sync(BigDFT_Inputs *in);

/*  Generic tools. */
gchar* _get_c_string(const gchar *fstr, guint len);

#endif
