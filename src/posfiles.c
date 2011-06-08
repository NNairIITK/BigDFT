/*
!> @file
!!  Routines to read and write position files (C part).
!! @author
!!    Copyright (C) 2007-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
*/

#include <config.h>

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#ifdef HAVE_LIB_ARCHIVE
#include <archive.h>
#include <archive_entry.h>
static struct archive *_posout_;
#endif

void FC_FUNC(addtocompress, ADDTOCOMPRESS)(const char *archive, int *lgAr,
                                           const char *filename, int *lgF)
{
  static int count = 0;
  char *arFilename, *addFilename;
#ifdef HAVE_LIB_ARCHIVE
  struct archive_entry *entry;
  struct stat st;
  char buff[8192];
  int len;
  int fd;
#endif

  arFilename = (char*)malloc(sizeof(char) * (*lgAr + 1));
  memcpy(arFilename, archive, (size_t)*lgAr);
  arFilename[*lgAr] = '\0';

  addFilename = (char*)malloc(sizeof(char) * (*lgF + 1));
  memcpy(addFilename, filename, (size_t)*lgF);
  addFilename[*lgF] = '\0';

  if (!count)
    {
#ifdef HAVE_LIB_ARCHIVE
      fprintf(stdout, "Init '%s'.\n", arFilename);
      _posout_ = archive_write_new();
      if (strstr(arFilename, ".tar.bz2"))
        archive_write_set_compression_bzip2(_posout_);
      else
        archive_write_set_compression_gzip(_posout_);
      archive_write_set_format_pax_restricted(_posout_);
      archive_write_open_filename(_posout_, arFilename);
#else
      fprintf(stdout, "Warning: no possibility of archive,"
              " compile BigDFT with libarchive.\n");
#endif
    }
  count +=1;

  fprintf(stdout, "Copy '%s' to '%s'.\n", addFilename, arFilename);
#ifdef HAVE_LIB_ARCHIVE
  stat(addFilename, &st);
  if (st.st_size > 0)
    {
      entry = archive_entry_new();
      archive_entry_set_pathname(entry, addFilename);
      /* archive_entry_copy_stat(entry, &st); */
      archive_entry_set_size(entry, st.st_size);
      archive_entry_set_filetype(entry, AE_IFREG);
      fprintf(stdout, "Write header to '%s'.\n", arFilename);
      archive_write_header(_posout_, entry);
      fprintf(stdout, "Open '%s'.\n", addFilename);
      fd = open(addFilename, O_RDONLY);
      fprintf(stdout, "Read '%s'.\n", addFilename);
      len = read(fd, buff, sizeof(buff));
      while ( len > 0 )
        {
          archive_write_data(_posout_, buff, len);
          len = read(fd, buff, sizeof(buff));
        }
      close(fd);
      archive_entry_free(entry);
    }
  fprintf(stdout, "Remove '%s'.\n", addFilename);
  unlink(addFilename);
#endif

  free(addFilename);
  free(arFilename);
}

void FC_FUNC(finalisecompress, FINALISECOMPRESS)(void)
{
#ifdef HAVE_LIB_ARCHIVE
  fprintf(stdout, "Finalise\n");
  archive_write_close(_posout_);
  archive_write_finish(_posout_);
#endif
}
