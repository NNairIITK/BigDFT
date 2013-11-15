/*
Contains only basic ABINIT macros:
Author: TRangel.
*/

#  define ABI_ALLOCATE(ARR,SIZE) \
   allocate(ARR SIZE)
#  define ABI_DEALLOCATE(ARR) \
   deallocate(ARR)
#  define ABI_DATATYPE_ALLOCATE(ARR,SIZE) \
   allocate(ARR SIZE)
#  define ABI_DATATYPE_DEALLOCATE(ARR) \
   deallocate(ARR)


#  define MSG_COMMENT(msg)      call msg_hndl(msg,"COMMENT","COLL")
#  define MSG_WARNING(msg)      call msg_hndl(msg,"WARNING","COLL")
#  define MSG_ERROR(msg)        call msg_hndl(msg,"ERROR"  ,"COLL")
#  define MSG_BUG(msg)          call msg_hndl(msg,"BUG"    ,"COLL")

