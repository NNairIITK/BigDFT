#include <yaml.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"

#ifndef FC_FUNC_
#define FC_FUNC_(A,B) A ## _
#endif

typedef struct _PosinpAtoms
{
  /* The cell. */
  unsigned int BC, Units;
  double acell[3], angdeg[3];

  /* The positions. */
  unsigned int nat, ntypes, units;
  double *rxyz;
  char **atomnames;
  unsigned int *iatype, *ifrztyp;
  int *igspin;
} PosinpAtoms;
static const char *UnitsPositions_keys[] = {"angstroem", "atomic", "bohr", "reduced", NULL};
static const char *BC_keys[] = {"periodic", "free", "surface", "wire", NULL};
static const char *Units_keys[] = {"angstroem", "atomic", "bohr", NULL};

static PosinpAtoms* posinp_atoms_new()
{
  PosinpAtoms *atoms;

  atoms = malloc(sizeof(PosinpAtoms));
  memset(atoms, 0, sizeof(PosinpAtoms));

  return atoms;
}
static void posinp_atoms_free(PosinpAtoms *atoms)
{
  unsigned int i;

  /* Freeing. */
  if (atoms->rxyz)
    free(atoms->rxyz);
  if (atoms->atomnames)
    {
      for (i = 0; i < atoms->ntypes; i++)
        if (atoms->atomnames[i])
          free(atoms->atomnames[i]);
      free(atoms->atomnames);
    }
  if (atoms->iatype)
    free(atoms->iatype);
  if (atoms->ifrztyp)
    free(atoms->ifrztyp);
  if (atoms->igspin)
    free(atoms->igspin);

  free(atoms);
}
#ifdef TEST_ME
static void posinp_atoms_trace(PosinpAtoms *atoms)
{
  unsigned int i;

  fprintf(stdout, "BC is %d (%s).\n", atoms->BC, BC_keys[atoms->BC]);
  fprintf(stdout, "Units are %d (%s).\n", atoms->Units, Units_keys[atoms->Units]);
  fprintf(stdout, "acell is %g %g %g.\n", atoms->acell[0], atoms->acell[1], atoms->acell[2]);
  fprintf(stdout, "angdeg is %g %g %g.\n", atoms->angdeg[0], atoms->angdeg[1], atoms->angdeg[2]);
  fprintf(stdout, "Position units are %d (%s).\n", atoms->units, UnitsPositions_keys[atoms->units]);
  for (i = 0; i < atoms->nat; i++)
    fprintf(stdout, "Atom '%s' at %g %g %g (%d) f:%d s:%d.\n", atoms->atomnames[atoms->iatype[i]],
            atoms->rxyz[3 * i + 0], atoms->rxyz[3 * i + 1], atoms->rxyz[3 * i + 2], atoms->iatype[i],
            atoms->ifrztyp[i], atoms->igspin[i]);
}
#endif

typedef struct _PosinpList PosinpList;
struct _PosinpList
{
  PosinpList *next;
  PosinpAtoms *data;
};

static void _yaml_parser_error(const yaml_parser_t *parser)
{
  switch (parser->error)
    {
    case YAML_MEMORY_ERROR:
      fprintf(stderr, "Memory error: Not enough memory for parsing\n");
      return;
    case YAML_READER_ERROR:
      if (parser->problem_value != -1)
        {
          fprintf(stderr, "Reader error: %s: #%X at %ld\n", parser->problem,
                  parser->problem_value, parser->problem_offset);
        }
      else
        {
          fprintf(stderr, "Reader error: %s at %ld\n", parser->problem,
                  parser->problem_offset);
        }
      return;
    case YAML_SCANNER_ERROR:
      if (parser->context)
        {
          fprintf(stderr, "Scanner error: %s at line %ld, column %ld\n"
                  "%s at line %ld, column %ld\n", parser->context,
                  parser->context_mark.line+1, parser->context_mark.column+1,
                  parser->problem, parser->problem_mark.line+1,
                  parser->problem_mark.column+1);
        }
      else
        {
          fprintf(stderr, "Scanner error: %s at line %ld, column %ld\n",
                  parser->problem, parser->problem_mark.line+1,
                  parser->problem_mark.column+1);
        }
      return;
    case YAML_PARSER_ERROR:
      if (parser->context)
        {
          fprintf(stderr, "Parser error: %s at line %ld, column %ld\n"
                  "%s at line %ld, column %ld\n", parser->context,
                  parser->context_mark.line+1, parser->context_mark.column+1,
                  parser->problem, parser->problem_mark.line+1,
                  parser->problem_mark.column+1);
        }
      else
        {
          fprintf(stderr, "Parser error: %s at line %ld, column %ld\n",
                  parser->problem, parser->problem_mark.line+1,
                  parser->problem_mark.column+1);
        }
      return;
    default:
      /* Couldn't happen. */
      fprintf(stderr, "Internal error\n");
      return;
    }
}

static int _yaml_parser_read_int(yaml_parser_t *parser, int *val)
{
  yaml_event_t event;
  int done;
  char *end;

  /* Read the value. */
  done = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SCALAR_EVENT)
        {
          *val = (int)strtol((const char*)event.data.scalar.value, &end, 10);
          if (end == (char*)event.data.scalar.value)
            {
              fprintf(stderr, "Parser error: cannot convert '%s' to an int.\n", end);
              done = -1;
            }
          else
            done = 0;
        }
      else
        {
          fprintf(stderr, "Parser error: value awaited after key '%s'.\n", "BC");
          done = (event.type == YAML_STREAM_END_EVENT)?1:-1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser);
      done = -1;
    }

  return done;
}
static int _yaml_parser_read_double(yaml_parser_t *parser, double *val)
{
  yaml_event_t event;
  int done;
  char *end;

  /* Read the value. */
  done = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SCALAR_EVENT)
        {
          *val = strtod((const char*)event.data.scalar.value, &end);
          if (end == (char*)event.data.scalar.value)
            {
              fprintf(stderr, "Parser error: cannot convert '%s' to a double.\n", end);
              done = -1;
            }
          else
            done = 0;
        }
      else
        {
          fprintf(stderr, "Parser error: value awaited after key '%s'.\n", "BC");
          done = (event.type == YAML_STREAM_END_EVENT)?1:-1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser);
      done = -1;
    }

  return done;
}
static int _yaml_parser_read_double_array(yaml_parser_t *parser, const char *key,
                                          double *vals, unsigned int n)
{
  yaml_event_t event;
  int done;
  unsigned int i;

  /* Read the value. */
  done = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SEQUENCE_START_EVENT)
        {
          for (i = 0; i < n && done == 0; i++)
            done = _yaml_parser_read_double(parser, vals + i);
          yaml_event_delete(&event);
          if (yaml_parser_parse(parser, &event))
            {
              if (event.type != YAML_SEQUENCE_END_EVENT)
                {
                  fprintf(stderr, "Parser error: end sequence missing for key '%s' after %d values.\n", key, n);
                  done = -1;
                }
            }
          else
            {
              /* Error treatment. */
              _yaml_parser_error(parser);
              done = -1;
            }
        }
      else
        {
          fprintf(stderr, "Parser error: sequence awaited after key '%s'.\n", key);
          done = -1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser);
      done = -1;
    }

  return done;
}
static int _yaml_parser_read_keyword(yaml_parser_t *parser, const char *key,
                                     const char *keys[], unsigned int *id)
{
  yaml_event_t event;
  int done;

  /* Read the value. */
  done = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SCALAR_EVENT)
        {
          for (*id = 0; keys[*id]; *id += 1)
            if (!strcasecmp((const char*)event.data.scalar.value, keys[*id]))
              break;
          if (keys[*id])
            done = 0;
          else
            {
              fprintf(stderr, "Parser error: cannot find key value '%s' in:\n",
                      event.data.scalar.value);
              for (*id = 0; keys[*id]; *id += 1)
                fprintf(stderr, "              - '%s'\n", keys[*id]);
            }
        }
      else
        {
          fprintf(stderr, "Parser error: value awaited after key '%s'.\n", key);
          done = -1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser);
      done = -1;
    }

  return done;
}








static int posinp_yaml_cell(yaml_parser_t *parser, PosinpAtoms *atoms)
{
  yaml_event_t event;
  int done, count;

  /* Default values. */
  atoms->BC = 1;
  atoms->Units = 1;
  atoms->acell[0] = 5.;
  atoms->acell[1] = 5.;
  atoms->acell[2] = 5.;
  atoms->angdeg[0] = 90.;
  atoms->angdeg[1] = 90.;
  atoms->angdeg[2] = 90.;

  /* Read the event sequence. */
  done = 0;
  count = 0;
  while (!done)
    /* Get the next event. */
    if (yaml_parser_parse(parser, &event))
      {
        switch(event.type)
          {
          case YAML_SEQUENCE_START_EVENT:
          case YAML_MAPPING_START_EVENT:
            count += 1;
            done = 0;
            break;
          case YAML_SEQUENCE_END_EVENT:
          case YAML_MAPPING_END_EVENT:
            count -= 1;
            done = 0;
            break;
          case YAML_SCALAR_EVENT:
            if (!strcmp((const char*)event.data.scalar.value, "BC"))
              done = _yaml_parser_read_keyword(parser, "BC", BC_keys, &atoms->BC);
            else if (!strcmp((const char*)event.data.scalar.value, "Units"))
              done = _yaml_parser_read_keyword(parser, "Units", Units_keys, &atoms->Units);
            else if (!strcmp((const char*)event.data.scalar.value, "acell"))
              done = _yaml_parser_read_double_array(parser, "acell", atoms->acell, 3);
            else if (!strcmp((const char*)event.data.scalar.value, "angdeg"))
              done = _yaml_parser_read_double_array(parser, "angdeg", atoms->angdeg, 3);
            else
              done = 0;
            break;
          default:
            done = (event.type == YAML_STREAM_END_EVENT);
            break;
          }

        /* Are we finished? */
        if (count == 0)
          done = 2;

        /* The application is responsible for destroying the event object. */
        yaml_event_delete(&event);
      }
    else
      {
        /* Error treatment. */
        _yaml_parser_error(parser);
        done = -1;
      }

  if (done == 2)
    done = 0;

  return done;
}
static int posinp_yaml_frozen(yaml_parser_t *parser, unsigned int *ifrztyp)
{
  yaml_event_t event;
  int done;

  /* Read the value. */
  done = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SCALAR_EVENT)
        {
          if (!strcasecmp((const char*)event.data.scalar.value, "No"))
            *ifrztyp = 0;
          else if (!strcasecmp((const char*)event.data.scalar.value, "Yes") ||
                   !strcasecmp((const char*)event.data.scalar.value, "f"))
            *ifrztyp = 1;
          else if (!strcasecmp((const char*)event.data.scalar.value, "fy"))
            *ifrztyp = 2;
          else if (!strcasecmp((const char*)event.data.scalar.value, "fxz"))
            *ifrztyp = 3;
          else 
            {
              fprintf(stderr, "Parser error: frozen keyword '%s' unknown.\n", event.data.scalar.value);
              done = -1;
            }
        }
      else
        {
          fprintf(stderr, "Parser error: atom name awaited.\n");
          done = -1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser);
      done = -1;
    }

  return done;
}
static int posinp_yaml_coord(yaml_parser_t *parser, double coords[3], char **names, unsigned int *iat)
{
  yaml_event_t event;
  int done;

  /* Read the value. */
  done = 0;
  *iat = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SCALAR_EVENT)
        {
          /* Here parse the name... */
          for (*iat = 0; names[*iat] && strcmp((const char*)names[*iat],
                                               (const char*)event.data.scalar.value); *iat += 1);
          if (!names[*iat])
            names[*iat] = strdup((const char*)event.data.scalar.value);
          /* Then the coordinates. */
          done = _yaml_parser_read_double_array(parser, names[*iat], coords, 3);
        }
      else
        {
          fprintf(stderr, "Parser error: atom name awaited.\n");
          done = -1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser);
      done = -1;
    }

  return done;
}
static int posinp_yaml_coords(yaml_parser_t *parser, PosinpAtoms *atoms)
{
  yaml_event_t event;
  int done;
  unsigned int count, atom_size;
#define ATOM_INC 100

  /* Read the event sequence. */
  done = 0;
  count = 0;
  atom_size = ATOM_INC;
  atoms->rxyz = malloc(sizeof(double) * atom_size * 3);
  atoms->atomnames = malloc(sizeof(char*) * atom_size);
  memset(atoms->atomnames, 0, sizeof(char*) * atom_size);
  atoms->iatype = malloc(sizeof(unsigned int) * atom_size);
  atoms->ifrztyp = malloc(sizeof(unsigned int) * atom_size);
  memset(atoms->ifrztyp, 0, sizeof(unsigned int) * atom_size);
  atoms->igspin = malloc(sizeof(int) * atom_size);
  memset(atoms->igspin, 0, sizeof(int) * atom_size);
  
  while (!done)
    /* Get the next event. */
    if (yaml_parser_parse(parser, &event))
      {
        switch(event.type)
          {
          case YAML_MAPPING_START_EVENT:
            /* Each mapping is one atom. */
            done = posinp_yaml_coord(parser, atoms->rxyz + 3 * count,
                                     atoms->atomnames, atoms->iatype + count);
            count += 1;
            if (count >= atom_size)
              {
                atom_size += ATOM_INC;
                atoms->rxyz = realloc(atoms->rxyz, sizeof(double) * 3 * atom_size);
                atoms->atomnames = realloc(atoms->atomnames, sizeof(char*) * atom_size);
                memset(atoms->atomnames - ATOM_INC, 0, sizeof(char*) * ATOM_INC);
                atoms->iatype = realloc(atoms->iatype, sizeof(unsigned int) * atom_size);
                atoms->ifrztyp = realloc(atoms->ifrztyp, sizeof(unsigned int) * atom_size);
                memset(atoms->ifrztyp - ATOM_INC, 0, sizeof(unsigned int) * ATOM_INC);
                atoms->igspin = realloc(atoms->igspin, sizeof(int) * atom_size);
                memset(atoms->igspin - ATOM_INC, 0, sizeof(int) * ATOM_INC);
              }
            done = 0;
            break;
          case YAML_SEQUENCE_END_EVENT:
            done = 2;
            break;
          case YAML_SCALAR_EVENT:
            if (!strcmp((const char*)event.data.scalar.value, "IGSpin"))
              done = _yaml_parser_read_int(parser, atoms->igspin + count - 1);
            if (!strcmp((const char*)event.data.scalar.value, "Frozen"))
              done = posinp_yaml_frozen(parser, atoms->ifrztyp + count - 1);
            done = 0;
            break;
          default:
            done = (event.type == YAML_STREAM_END_EVENT);
            break;
          }

        /* The application is responsible for destroying the event object. */
        yaml_event_delete(&event);
      }
    else
      {
        /* Error treatment. */
        _yaml_parser_error(parser);
        done = -1;
      }
  
  if (done == 2)
    done = 0;

  atoms->nat       = count;
  atoms->rxyz      = realloc(atoms->rxyz,      sizeof(double) * 3 * atoms->nat);
  atoms->iatype    = realloc(atoms->iatype,    sizeof(unsigned int) * atoms->nat);
  atoms->ifrztyp   = realloc(atoms->ifrztyp,   sizeof(unsigned int) * atoms->nat);
  atoms->igspin    = realloc(atoms->igspin,    sizeof(int) * atoms->nat);
  for (atoms->ntypes = 0; atoms->atomnames[atoms->ntypes]; atoms->ntypes++);
  atoms->atomnames = realloc(atoms->atomnames, sizeof(char*) * atoms->ntypes);

  return done;
}
static int posinp_yaml_position(yaml_parser_t *parser, PosinpAtoms *atoms)
{
  yaml_event_t event;
  int done, count;
  unsigned int Units;

  /* Default values. */
  Units = 1;

  /* Read the event sequence. */
  done = 0;
  count = 0;
  while (!done)
    /* Get the next event. */
    if (yaml_parser_parse(parser, &event))
      {
        switch(event.type)
          {
          case YAML_SEQUENCE_START_EVENT:
            if (count == 0)
              {
                done = posinp_yaml_coords(parser, atoms);
                break;
              }
          case YAML_MAPPING_START_EVENT:
            count += 1;
            done = 0;
            break;
          case YAML_SEQUENCE_END_EVENT:
          case YAML_MAPPING_END_EVENT:
            count -= 1;
            done = 0;
            break;
          case YAML_SCALAR_EVENT:
            if (!strcmp((const char*)event.data.scalar.value, "Units"))
              done = _yaml_parser_read_keyword(parser, "Units", UnitsPositions_keys, &atoms->units);
            else if (!strcmp((const char*)event.data.scalar.value, "Values"))
              done = posinp_yaml_coords(parser, atoms);
            else
              done = 0;
            break;
          default:
            done = (event.type == YAML_STREAM_END_EVENT);
            break;
          }

        /* Are we finished? */
        if (count == 0)
          done = 2;

        /* The application is responsible for destroying the event object. */
        yaml_event_delete(&event);
      }
    else
      {
        /* Error treatment. */
        _yaml_parser_error(parser);
        done = -1;
      }

  if (done == 2)
    done = 0;

  return done;
}
PosinpList* posinp_yaml_parse(const char *filename)
{
  PosinpList *list, *tmp, *tmp2;
  FILE *input;
  yaml_parser_t parser;
  yaml_event_t event;
  int done;
  
  list = tmp = (PosinpList*)0;
  
  input = fopen(filename, "rb");
  if (!input)
    return list;

  /* Create the Parser object. */
  yaml_parser_initialize(&parser);
  yaml_parser_set_input_file(&parser, input);

  /* Read the event sequence. */
  done = 0;
  while (!done)
    /* Get the next event. */
    if (yaml_parser_parse(&parser, &event))
      {
        if (event.type == YAML_DOCUMENT_START_EVENT)
          {
            tmp2 = malloc(sizeof(PosinpList));
            memset(tmp2, 0, sizeof(PosinpList));
            tmp2->data = posinp_atoms_new();
            
            if (!tmp)
              list = tmp = tmp2;
            else
              {
                tmp->next = tmp2;
                tmp = tmp2;
              }
          }
        else if (event.type == YAML_SCALAR_EVENT &&
                 !strcmp((const char*)event.data.scalar.value, "Cell"))
          done = posinp_yaml_cell(&parser, tmp->data);
        else if (event.type == YAML_SCALAR_EVENT &&
                 !strcmp((const char*)event.data.scalar.value, "Positions"))
          done = posinp_yaml_position(&parser, tmp->data);
        else
          done = (event.type == YAML_STREAM_END_EVENT);

        /* The application is responsible for destroying the event object. */
        yaml_event_delete(&event);
      }
    else
      {
        /* Error treatment. */
        _yaml_parser_error(&parser);
        done = -1;
      }

  /* Destroy the Parser object. */
  yaml_parser_delete(&parser);
  fclose(input);

  return list;
}
void FC_FUNC_(posinp_yaml_parse, POSINP_YAML_PARSE)(PosinpList **self,
                                                    const char *filename, unsigned int *ln)
{
  char *name;

  name = malloc(sizeof(char) * (*ln + 1));
  memcpy(name, filename, sizeof(char) * *ln);
  name[*ln] = '\0';

  *self = posinp_yaml_parse(name);

  free(name);
}
static void posinp_yaml_free_list(PosinpList *lst)
{
  PosinpList *tmp;

  while(lst)
    {
      posinp_atoms_free(lst->data);
      tmp = lst;
      lst = lst->next;
      free(tmp);
    }
}
void FC_FUNC_(posinp_yaml_free_list, POSINP_YAML_FREE_LIST)(PosinpList **self)
{
  posinp_yaml_free_list(*self);
}

void FC_FUNC_(posinp_yaml_get_cell, POSINP_YAML_GET_CELL)(PosinpList **self, unsigned int *i,
                                                          unsigned int *BC, unsigned int *Units,
                                                          double acell[3], double angdeg[3])
{
  PosinpList *lst;
  unsigned int j;

  for (lst = *self, j = 0; j < *i; j++)
    if (lst)
      lst = lst->next;

  if (lst)
    {
      *BC = lst->data->BC;
      *Units = lst->data->Units;
      acell[0] = lst->data->acell[0];
      acell[1] = lst->data->acell[1];
      acell[2] = lst->data->acell[2];
      angdeg[0] = lst->data->angdeg[0];
      angdeg[1] = lst->data->angdeg[1];
      angdeg[2] = lst->data->angdeg[2];
    }
}
void FC_FUNC_(posinp_yaml_get_dims, POSINP_YAML_GET_DIMS)(PosinpList **self, unsigned int *i,
                                                          unsigned int *nat, unsigned int *ntypes)
{
  PosinpList *lst;
  unsigned int j;

  for (lst = *self, j = 0; j < *i; j++)
    if (lst)
      lst = lst->next;

  if (lst)
    {
      *nat = lst->data->nat;
      *ntypes = lst->data->ntypes;
    }
}
void FC_FUNC_(posinp_yaml_get_atoms, POSINP_YAML_GET_ATOMS)(PosinpList **self, unsigned int *i,
                                                            unsigned int *units, double *rxyz, 
                                                            unsigned int *iatype, unsigned int *ifrztyp,
                                                            int *igspin)
{
  PosinpList *lst;
  unsigned int j;

  for (lst = *self, j = 0; j < *i; j++)
    if (lst)
      lst = lst->next;

  if (lst)
    {
      *units = lst->data->units;
      memcpy(iatype,  lst->data->iatype,  sizeof(unsigned int) * lst->data->nat);
      memcpy(ifrztyp, lst->data->ifrztyp, sizeof(unsigned int) * lst->data->nat);
      memcpy(igspin,  lst->data->igspin,  sizeof(unsigned int) * lst->data->nat);
      for (j = 0; j < lst->data->nat; j++)
        {
          iatype[j] += 1;
          rxyz[0 * lst->data->nat + j] = lst->data->rxyz[j * 3 + 0];
          rxyz[1 * lst->data->nat + j] = lst->data->rxyz[j * 3 + 1];
          rxyz[2 * lst->data->nat + j] = lst->data->rxyz[j * 3 + 2];
        }
    }
}
void FC_FUNC_(posinp_yaml_get_atomname, POSINP_YAML_GET_ATOMNAME)(PosinpList **self, unsigned int *i,
                                                                  unsigned int *ityp, char name[20])
{
  PosinpList *lst;
  unsigned int j, ln;

  for (lst = *self, j = 0; j < *i; j++)
    if (lst)
      lst = lst->next;

  if (lst)
    {
      memset(name, ' ', sizeof(char) * 20);
      ln = strlen(lst->data->atomnames[*ityp]);
      memcpy(name, lst->data->atomnames[*ityp], sizeof(char) * ((ln > 20)?20:ln));
    }
}

#ifdef TEST_ME
int main(int argc, const char **argv)
{
  PosinpList *lst, *tmp;

  lst = posinp_yaml_parse(argv[1]);
  for (tmp = lst; tmp; tmp = tmp->next)
    {
      fprintf(stdout, "---\n");
      posinp_atoms_trace(tmp->data);
    }
  posinp_yaml_free_list(lst);

  return 0;
}
#endif
