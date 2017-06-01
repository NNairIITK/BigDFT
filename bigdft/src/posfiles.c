/*
*  @file Routines to read and write position files (C part).
*
*  @author
*    Copyright (C) 2007-2011 BigDFT group 
*    This file is distributed under the terms of the
*    GNU General Public License, see ~/COPYING file
*    or http://www.gnu.org/copyleft/gpl.txt .
*    For the list of contributors, see ~/AUTHORS 
*/

#include <config.h>

#define _GNU_SOURCE

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#ifdef HAVE_LIB_ARCHIVE
#include <archive.h>
#include <archive_entry.h>
static struct archive *_posout_ = NULL, *_posinp_ = NULL;
#endif

#ifndef HAVE_STRNDUP
char* strndup(const char *src, size_t len)
{
  char *out;

  out = malloc(sizeof(char) * (len + 1));
  memcpy(out, src, sizeof(char) * len);
  out[len] = '\0';

  return out;
}
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

  arFilename = strndup(archive, (size_t)*lgAr);
  addFilename = strndup(filename, (size_t)*lgF);

  if (!count)
    {
#ifdef HAVE_LIB_ARCHIVE
      /* fprintf(stdout, "Init '%s'.\n", arFilename); */
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

  /* fprintf(stdout, "Copy '%s' to '%s'.\n", addFilename, arFilename); */
#ifdef HAVE_LIB_ARCHIVE
  stat(addFilename, &st);
  if (st.st_size > 0)
    {
      entry = archive_entry_new();
      archive_entry_set_pathname(entry, addFilename);
      /* archive_entry_copy_stat(entry, &st); */
      archive_entry_set_size(entry, st.st_size);
      archive_entry_set_filetype(entry, AE_IFREG);
      archive_entry_set_mode(entry, st.st_mode);
      archive_entry_set_mtime(entry, st.st_mtime, 0);
      archive_write_header(_posout_, entry);
      fd = open(addFilename, O_RDONLY);
      len = read(fd, buff, sizeof(buff));
      while ( len > 0 )
        {
          archive_write_data(_posout_, buff, len);
          len = read(fd, buff, sizeof(buff));
        }
      close(fd);
      archive_entry_free(entry);
    }
  /* fprintf(stdout, "Remove '%s'.\n", addFilename); */
  unlink(addFilename);
#endif

  free(addFilename);
  free(arFilename);
}

void FC_FUNC(finalisecompress, FINALISECOMPRESS)(void)
{
#ifdef HAVE_LIB_ARCHIVE
  /* fprintf(stdout, "Finalise\n"); */
  if (_posout_)
    {
      archive_write_close(_posout_);
      archive_write_finish(_posout_);
    }
#endif
}

static void _openArchive_(const char *arFilename)
{
  static int count = 0;

  if (!count)
    {
#ifdef HAVE_LIB_ARCHIVE
      /* fprintf(stdout, "Init '%s'.\n", arFilename); */
      _posinp_ = archive_read_new();
      archive_read_support_compression_all(_posinp_);
      archive_read_support_format_all(_posinp_);
      archive_read_open_filename(_posinp_, arFilename, 10240);
#else
      fprintf(stdout, "Warning: no possibility of archive,"
              " compile BigDFT with libarchive.\n");
#endif
    }
  count +=1;
}

void FC_FUNC(extractnextcompress, EXTRACTNEXTCOMPRESS)(const char *archive, int *lgAr,
                                                       const char *filename, int *lgF,
                                                       int *extract, char ext[6])
{
  char *arFilename, *addFilename;
#ifdef HAVE_LIB_ARCHIVE
  int len;
  struct archive_entry *entry;
  const void *buff;
  size_t size;
  ssize_t res;
  off_t offset;
  int fd;
#endif

  arFilename = strndup(archive, (size_t)*lgAr);
  addFilename = strndup(filename, (size_t)*lgF);
#ifdef HAVE_LIB_ARCHIVE
  len = strlen(addFilename);
#endif

  *extract = 0;
  _openArchive_(arFilename);

  /* fprintf(stdout, "Extract '%s' from '%s'.\n", addFilename, arFilename); */
#ifdef HAVE_LIB_ARCHIVE
  while (archive_read_next_header(_posinp_, &entry) == ARCHIVE_OK)
    if (!strncmp(addFilename, archive_entry_pathname(entry), len))
      {
        /* We extract the file. */
        fd = creat(archive_entry_pathname(entry), 0640);
        /* fprintf(stdout, "Create '%s' (%d).\n", archive_entry_pathname(entry), fd); */
        res = 1;
        while (archive_read_data_block(_posinp_, &buff, &size, &offset) == ARCHIVE_OK &&
               res > 0)
          res = write(fd, buff, size);
        close(fd);
        *extract = 1;
        strcpy(ext, archive_entry_pathname(entry) + len + 1);
        memset(ext + strlen(ext), ' ', 6 - strlen(ext));
        break;
      }
#endif

  free(addFilename);
  free(arFilename);
}

void FC_FUNC(unlinkextract, UNLINKEXTRACT)(const char *filename, int *lgF)
{
  char *addFilename;

  addFilename = (char*)malloc(sizeof(char) * (*lgF + 1));
  memcpy(addFilename, filename, (size_t)*lgF);
  addFilename[*lgF] = '\0';

#ifdef HAVE_LIB_ARCHIVE
  unlink(addFilename);
#endif

  free(addFilename);
}

void FC_FUNC(finaliseextract, FINALISEEXTRACT)(void)
{
#ifdef HAVE_LIB_ARCHIVE
  /* fprintf(stdout, "Finalise\n"); */
  if (_posinp_)
    {
      archive_write_close(_posinp_);
      archive_write_finish(_posinp_);
    }
#endif
}

#ifdef HAVE_LIB_ARCHIVE
static const void *buff = (void*)0;
static size_t size;
static size_t pos = 0;
#endif

void FC_FUNC(opennextcompress, OPENNEXTCOMPRESS)(const char *archive, int *lgAr,
                                                 const char *filename, int *lgF,
                                                 int *extract, char ext[6])
{
  char *arFilename, *addFilename;
#ifdef HAVE_LIB_ARCHIVE
  int len;
  struct archive_entry *entry;
#endif

  arFilename = strndup(archive, (size_t)*lgAr);
  addFilename = strndup(filename, (size_t)*lgF);
#ifdef HAVE_LIB_ARCHIVE
  len = strlen(addFilename);
#endif

  *extract = 0;
  _openArchive_(arFilename);

  /* fprintf(stdout, "Open '%s' from '%s'.\n", addFilename, arFilename); */
#ifdef HAVE_LIB_ARCHIVE

  while (archive_read_next_header(_posinp_, &entry) == ARCHIVE_OK)
    if (!strncmp(addFilename, archive_entry_pathname(entry), len))
      {
        *extract = 1;
        strcpy(ext, archive_entry_pathname(entry) + len + 1);
        memset(ext + strlen(ext), ' ', 6 - strlen(ext));

        buff = (void*)0;
        size = 0;
        pos = 0;
        break;
      }
#endif

  free(addFilename);
  free(arFilename);
}

void FC_FUNC(extractnextline, EXTRACTNEXTLINE)(char line[150], int *eof)
{
#ifdef HAVE_LIB_ARCHIVE
  char *vals;
  size_t  idx, jdx;
  off_t offset;
  int r;

  *eof = 0;
  jdx = 0;
  /* fprintf(stdout, "Read a line (buff = %p, pos = %d)...\n", buff, pos); */

  do
    {
      if (!buff)
        {
          r = archive_read_data_block(_posinp_, &buff, &size, &offset);
          if (r != ARCHIVE_OK)
            {
              buff = (void*)0;
              pos = 0;
              if (!jdx)
                *eof = 1;
              return;
            }
        }
      vals = (char*)buff;
      for (idx = pos; idx < size; idx++)
        {
          if (vals[idx] != '\n')
            line[jdx++] = vals[idx];
          else
            {
              line[jdx] = '\0';
              /* fprintf(stdout, "Found '%s'.\n", line); */
              memset(line + jdx, ' ', 150 - jdx);
              pos = idx + 1;
              return;
            }
        }
      /* We eat all. */
      buff = (void*)0;
      pos = 0;
    }
  while(1);
#else
  fprintf(stdout, "Warning: no possibility of archive,"
          " compile BigDFT with libarchive.\n");
#endif
}
