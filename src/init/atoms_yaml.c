#define _ISOC99_SOURCE
#include "config.h"

#ifdef HAVE_YAML
#include <yaml.h>
#endif

#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atoms_yaml.h"

#ifndef FC_FUNC_
#define FC_FUNC_(A,B) A ## _
#endif

#ifdef HAVE_YAML
static const char *UnitsPositions_keys[] = {"bohr", "angstroem", "reduced", "atomic", NULL};
static const char *BC_keys[] = {"free", "wire", "surface", "periodic", NULL};
static const char *Units_keys[] = {"bohr", "angstroem", "atomic", NULL};
static const char *funits_keys[] = {"Ha/Bohr", "eV/Ang", NULL};
static const char *eunits_keys[] = {"Ha", "eV", "Ry", NULL};
static const char *frozen_keys[] = {"No", "Yes", "fy", "fxz",
                                    "N",  "Y",   "",   "",
                                    "false", "true", "", "",
                                    "off", "on", "", "", NULL};
static const char *bool_keys[] = {"No", "Yes",
                                  "N",  "Y",
                                  "false", "true",
                                  "off", "on", NULL};

static PosinpAtoms* posinp_atoms_new()
{
  PosinpAtoms *atoms;

  atoms = malloc(sizeof(PosinpAtoms));
  memset(atoms, 0, sizeof(PosinpAtoms));

  /* Default values. */
  atoms->BC = POSINP_BC_FREE;
  atoms->Units = POSINP_CELL_UNITS_BOHR;
  atoms->angdeg[0] = 90.;
  atoms->angdeg[1] = 90.;
  atoms->angdeg[2] = 90.;

  /* Default values. */
  atoms->units = POSINP_COORD_UNITS_BOHR;

  /* Default values. */
  atoms->funits = POSINP_FORCE_UNITS_HARTREE_PER_BOHR;

  /* Default values. */
  atoms->energy = NAN;
  atoms->eunits = POSINP_ENERG_UNITS_HARTREE;
  atoms->gnrm_wfn = NAN;

  return atoms;
}
#endif

static void posinp_atoms_free(PosinpAtoms *atoms)
{
  unsigned int i;

  /* Freeing. */
  if (atoms->comment)
    free(atoms->comment);
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
  if (atoms->igchg)
    free(atoms->igchg);
  if (atoms->fxyz)
    free(atoms->fxyz);

  free(atoms);
}
#ifdef TEST_ME
static void posinp_atoms_trace(PosinpAtoms *atoms)
{
  unsigned int i;

  fprintf(stdout, "'%s'.\n", atoms->comment);
  if (!isnan(atoms->gnrm_wfn))
    fprintf(stdout, "Converged %d (%g).\n", atoms->converged, atoms->gnrm_wfn);
  if (!isnan(atoms->energy))
    fprintf(stdout, "Energy is %g %s (%d).\n", atoms->energy, eunits_keys[atoms->eunits], atoms->eunits);
  fprintf(stdout, "BC is %d (%s).\n", atoms->BC, BC_keys[atoms->BC]);
  fprintf(stdout, "Units are %d (%s).\n", atoms->Units, Units_keys[atoms->Units]);
  fprintf(stdout, "acell is %g %g %g.\n", atoms->acell[0], atoms->acell[1], atoms->acell[2]);
  fprintf(stdout, "angdeg is %g %g %g.\n", atoms->angdeg[0], atoms->angdeg[1], atoms->angdeg[2]);
  fprintf(stdout, "Position units are %d (%s).\n", atoms->units, UnitsPositions_keys[atoms->units]);
  for (i = 0; i < atoms->nat; i++)
    fprintf(stdout, "Atom '%s' at %g %g %g (%d) f:%d s:%d c:%d.\n", atoms->atomnames[atoms->iatype[i]],
            atoms->rxyz[3 * i + 0], atoms->rxyz[3 * i + 1], atoms->rxyz[3 * i + 2], atoms->iatype[i],
            atoms->ifrztyp[i], atoms->igspin[i], atoms->igchg[i]);
  if (atoms->fxyz)
    {
      fprintf(stdout, "Forces units are %d (%s).\n", atoms->funits, funits_keys[atoms->funits]);
      fprintf(stdout, "fnrm is %g, maxval is %g.\n", atoms->fnrm, atoms->maxval);
      for (i = 0; i < atoms->nat; i++)
        fprintf(stdout, "Force on '%s' is %g %g %g.\n", atoms->atomnames[atoms->iatype[i]],
                atoms->fxyz[3 * i + 0], atoms->fxyz[3 * i + 1], atoms->fxyz[3 * i + 2]);
    }
}
#endif

#define set_error(...)                          \
  if (message && !*message)                     \
    {                                           \
      ln = snprintf(NULL, 0, __VA_ARGS__);      \
      *message = malloc(sizeof(char) * ln);     \
      sprintf(*message, __VA_ARGS__);           \
    }                                           \
  else                                          \
    fprintf(stderr, __VA_ARGS__)

#ifdef HAVE_YAML

static void _yaml_parser_error(const yaml_parser_t *parser, char **message)
{
  size_t ln;

  switch (parser->error)
    {
    case YAML_MEMORY_ERROR:
      set_error("Memory error: Not enough memory for parsing\n");
      return;
    case YAML_READER_ERROR:
      if (parser->problem_value != -1)
        {
          set_error("Reader error: %s: #%X at %ld\n", parser->problem,
                    parser->problem_value, parser->problem_offset);
        }
      else
        {
          set_error("Reader error: %s at %ld\n", parser->problem,
                    parser->problem_offset);
        }
      return;
    case YAML_SCANNER_ERROR:
      if (parser->context)
        {
          set_error("Scanner error: %s at line %ld, column %ld\n"
                    "%s at line %ld, column %ld\n", parser->context,
                    parser->context_mark.line+1, parser->context_mark.column+1,
                    parser->problem, parser->problem_mark.line+1,
                    parser->problem_mark.column+1);
        }
      else
        {
          set_error("Scanner error: %s at line %ld, column %ld\n",
                    parser->problem, parser->problem_mark.line+1,
                    parser->problem_mark.column+1);
        }
      return;
    case YAML_PARSER_ERROR:
      if (parser->context)
        {
          set_error("Parser error: %s at line %ld, column %ld\n"
                    "%s at line %ld, column %ld\n", parser->context,
                    parser->context_mark.line+1, parser->context_mark.column+1,
                    parser->problem, parser->problem_mark.line+1,
                    parser->problem_mark.column+1);
        }
      else
        {
          set_error("Parser error: %s at line %ld, column %ld\n",
                    parser->problem, parser->problem_mark.line+1,
                    parser->problem_mark.column+1);
        }
      return;
    default:
      /* Couldn't happen. */
      set_error("Internal error\n");
      return;
    }
}

static int _yaml_parser_copy_str(yaml_parser_t *parser, char **val, char **message)
{
  yaml_event_t event;
  int done;
  size_t ln;

  /* Read the value. */
  done = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SCALAR_EVENT)
        {
          ln = strlen((const char*)event.data.scalar.value);
          *val = malloc(sizeof(char) * (ln + 1));
          memcpy(*val, event.data.scalar.value, sizeof(char) * ln);
          done = 0;
        }
      else
        {
          set_error("Parser error: value awaited.\n");
          done = (event.type == YAML_STREAM_END_EVENT)?1:-1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser, message);
      done = -1;
    }

  return done;
}
static int _yaml_parser_read_int(yaml_parser_t *parser, int *val, char **message)
{
  yaml_event_t event;
  int done;
  char *end;
  size_t ln;

  /* Read the value. */
  done = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SCALAR_EVENT)
        {
          *val = (int)strtol((const char*)event.data.scalar.value, &end, 10);
          if (end == (char*)event.data.scalar.value)
            {
              set_error("Parser error: cannot convert '%s' to an int.\n", end);
              done = -1;
            }
          else
            done = 0;
        }
      else
        {
          set_error("Parser error: value awaited.\n");
          done = (event.type == YAML_STREAM_END_EVENT)?1:-1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser, message);
      done = -1;
    }

  return done;
}
static int _yaml_parser_read_double(yaml_parser_t *parser, double *val, char **message)
{
  yaml_event_t event;
  int done;
  char *end;
  size_t ln;

  /* Read the value. */
  done = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SCALAR_EVENT)
        {
          if (!event.data.scalar.value[0])
            *val = strtod("NaN", &end);
          else if (!strcasecmp((const char*)event.data.scalar.value, ".inf"))
            *val = strtod((const char*)event.data.scalar.value + 1, &end);
          else
            *val = strtod((const char*)event.data.scalar.value, &end);
          if (end == (char*)event.data.scalar.value)
            {
              set_error("Parser error: cannot convert '%s' to a double.\n", end);
              done = -1;
            }
          else
            done = 0;
        }
      else
        {
          set_error("Parser error: value awaited.\n");
          done = (event.type == YAML_STREAM_END_EVENT)?1:-1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser, message);
      done = -1;
    }

  return done;
}
static int _yaml_parser_read_double_array(yaml_parser_t *parser, const char *key,
                                          double *vals, unsigned int n, char **message)
{
  yaml_event_t event;
  int done;
  unsigned int i;
  size_t ln;

  /* Read the value. */
  done = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SEQUENCE_START_EVENT)
        {
          for (i = 0; i < n && done == 0; i++)
            done = _yaml_parser_read_double(parser, vals + i, message);
          yaml_event_delete(&event);
          if (yaml_parser_parse(parser, &event))
            {
              if (event.type != YAML_SEQUENCE_END_EVENT)
                {
                  set_error("Parser error: end sequence missing for key '%s' after %d values.\n", key, n);
                  done = -1;
                }
            }
          else
            {
              /* Error treatment. */
              _yaml_parser_error(parser, message);
              done = -1;
            }
        }
      else
        {
          set_error("Parser error: sequence awaited after key '%s'.\n", key);
          done = -1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser, message);
      done = -1;
    }

  return done;
}
static int _find_keyword(const char *keys[], const char *value, unsigned int *id,
                         unsigned int modulo, char **message)
{
  int done;
  size_t ln;

  for (*id = 0; keys[*id]; *id += 1)
    if (!strcasecmp(value, keys[*id]))
      break;
  if (keys[*id])
    {
      *id = *id % modulo;
      done = 0;
    }
  else
    {
      *id = 0;
      set_error("Parser error: cannot find key value '%s'.\n", value);
      done = -1;
    }
  return done;
}
static int _find_units(const char *keys[], const char *value, unsigned int *id,
                       unsigned int modulo, char **message)
{
  int done;
  size_t ln;
  char *start, *end, *unit;

  *id = 0;
  start = strchr(value, '(');
  if (!start)
    /* No unit specified, no error. */
    return 0;
  end = strchr(start, ')');
  if (!start)
    {
      /* Parentethis not closed, error. */
      set_error("Parser error: unit not properly written in '%s'.\n", value);
      return -1;
    }
  ln = end - start - 1;
  unit = malloc(sizeof(char) * (ln + 1));
  memcpy(unit, start + 1, ln);
  unit[ln] = '\0';

  done = _find_keyword(keys, unit, id, modulo, message);

  free(unit);
  return done;
}
static int _yaml_parser_read_keyword(yaml_parser_t *parser, const char *key,
                                     const char *keys[], unsigned int *id,
                                     unsigned int modulo, char **message)
{
  yaml_event_t event;
  int done;
  size_t ln;

  /* Read the value. */
  done = 0;
  if (yaml_parser_parse(parser, &event))
    {
      if (event.type == YAML_SCALAR_EVENT)
        done = _find_keyword(keys, (const char*)event.data.scalar.value, id, modulo, message);
      else
        {
          set_error("Parser error: value awaited after key '%s'.\n", key);
          done = -1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser, message);
      done = -1;
    }

  return done;
}








static int posinp_yaml_cell(yaml_parser_t *parser, PosinpAtoms *atoms, char **message)
{
  yaml_event_t event;
  int done, count;

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
              done = _yaml_parser_read_keyword(parser, "BC", BC_keys, &atoms->BC, POSINP_N_BC, message);
            else if (!strcmp((const char*)event.data.scalar.value, "Units"))
              done = _yaml_parser_read_keyword(parser, "Units", Units_keys,
                                               &atoms->Units, POSINP_CELL_N_UNITS, message);
            else if (!strcmp((const char*)event.data.scalar.value, "acell"))
              done = _yaml_parser_read_double_array(parser, "acell", atoms->acell, 3, message);
            else if (!strcmp((const char*)event.data.scalar.value, "angdeg"))
              done = _yaml_parser_read_double_array(parser, "angdeg", atoms->angdeg, 3, message);
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
        _yaml_parser_error(parser, message);
        done = -1;
      }

  if (done == 2)
    done = 0;

  return done;
}
static int posinp_yaml_coord(yaml_parser_t *parser, double coords[3],
                             char **names, unsigned int *iat, char **message)
{
  yaml_event_t event;
  unsigned int ln;
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
            {
              ln = strlen((const char*)event.data.scalar.value);
              names[*iat] = malloc(sizeof(char*) * (ln + 1));
              memcpy(names[*iat], (const char*)event.data.scalar.value, sizeof(char*) * ln);
              names[*iat][ln] = '\0';
            }
          /* Then the coordinates. */
          done = _yaml_parser_read_double_array(parser, names[*iat], coords, 3, message);
        }
      else
        {
          set_error("Parser error: atom name awaited.\n");
          done = -1;
        }

      /* The application is responsible for destroying the event object. */
      yaml_event_delete(&event);
    }
  else
    {
      /* Error treatment. */
      _yaml_parser_error(parser, message);
      done = -1;
    }

  return done;
}
static int posinp_yaml_coords(yaml_parser_t *parser, PosinpAtoms *atoms, char **message)
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
  atoms->igchg = malloc(sizeof(int) * atom_size);
  memset(atoms->igchg, 0, sizeof(int) * atom_size);
  
  while (!done)
    /* Get the next event. */
    if (yaml_parser_parse(parser, &event))
      {
        switch(event.type)
          {
          case YAML_MAPPING_START_EVENT:
            /* Each mapping is one atom. */
            done = posinp_yaml_coord(parser, atoms->rxyz + 3 * count,
                                     atoms->atomnames, atoms->iatype + count, message);
            count += 1;
            if (count >= atom_size)
              {
                atom_size += ATOM_INC;
                atoms->rxyz = realloc(atoms->rxyz, sizeof(double) * 3 * atom_size);
                atoms->atomnames = realloc(atoms->atomnames, sizeof(char*) * atom_size);
                memset(atoms->atomnames + atom_size - ATOM_INC, 0, sizeof(char*) * ATOM_INC);
                atoms->iatype = realloc(atoms->iatype, sizeof(unsigned int) * atom_size);
                atoms->ifrztyp = realloc(atoms->ifrztyp, sizeof(unsigned int) * atom_size);
                memset(atoms->ifrztyp + atom_size - ATOM_INC, 0, sizeof(unsigned int) * ATOM_INC);
                atoms->igspin = realloc(atoms->igspin, sizeof(int) * atom_size);
                memset(atoms->igspin + atom_size - ATOM_INC, 0, sizeof(int) * ATOM_INC);
                atoms->igchg = realloc(atoms->igchg, sizeof(int) * atom_size);
                memset(atoms->igchg + atom_size - ATOM_INC, 0, sizeof(int) * ATOM_INC);
              }
            break;
          case YAML_SEQUENCE_END_EVENT:
            done = 2;
            break;
          case YAML_SCALAR_EVENT:
            if (!strcmp((const char*)event.data.scalar.value, "IGSpin"))
              done = _yaml_parser_read_int(parser, atoms->igspin + count - 1, message);
            else if (!strcmp((const char*)event.data.scalar.value, "IGChg"))
              done = _yaml_parser_read_int(parser, atoms->igchg + count - 1, message);
            else if (!strcmp((const char*)event.data.scalar.value, "Frozen"))
              done = _yaml_parser_read_keyword(parser, "Frozen", frozen_keys,
                                               atoms->ifrztyp + count - 1, POSINP_N_FROZEN, message);
            else
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
        _yaml_parser_error(parser, message);
        done = -1;
      }
  
  if (done == 2)
    done = 0;

  atoms->nat       = count;
  atoms->rxyz      = realloc(atoms->rxyz,      sizeof(double) * 3 * atoms->nat);
  atoms->iatype    = realloc(atoms->iatype,    sizeof(unsigned int) * atoms->nat);
  atoms->ifrztyp   = realloc(atoms->ifrztyp,   sizeof(unsigned int) * atoms->nat);
  atoms->igspin    = realloc(atoms->igspin,    sizeof(int) * atoms->nat);
  atoms->igchg     = realloc(atoms->igchg,     sizeof(int) * atoms->nat);
  for (atoms->ntypes = 0; atoms->atomnames[atoms->ntypes]; atoms->ntypes++);
  atoms->atomnames = realloc(atoms->atomnames, sizeof(char*) * atoms->ntypes);

  return done;
}
static int posinp_yaml_position(yaml_parser_t *parser, PosinpAtoms *atoms, char **message)
{
  yaml_event_t event;
  int done, count;

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
                done = posinp_yaml_coords(parser, atoms, message);
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
              done = _yaml_parser_read_keyword(parser, "Units", UnitsPositions_keys,
                                               &atoms->units, POSINP_COORD_N_UNITS, message);
            else if (!strcmp((const char*)event.data.scalar.value, "Values"))
              done = posinp_yaml_coords(parser, atoms, message);
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
        _yaml_parser_error(parser, message);
        done = -1;
      }

  if (done == 2)
    done = 0;

  return done;
}
static int posinp_yaml_force(yaml_parser_t *parser, PosinpAtoms *atoms, char **message)
{
  yaml_event_t event, event2;
  int done;
  unsigned int count;
  size_t ln;

  if (atoms->nat < 1)
    {
      set_error("Parser error: forces are defined before atoms.\n");
      done = -1;
    }

  /* Read the event sequence. */
  done = 0;
  count = 0;

  atoms->fxyz = malloc(sizeof(double) * atoms->nat * 3);
  memset(atoms->fxyz, 0, sizeof(double) * atoms->nat * 3);
  
  while (!done)
    /* Get the next event. */
    if (yaml_parser_parse(parser, &event))
      {
        switch(event.type)
          {
          case YAML_MAPPING_START_EVENT:
            /* Each mapping is one atom. */
            if (count >= atoms->nat)
              {
                set_error("Parser error: there are more forces than actual atoms.\n");
                done = -1;
                break;
              }
            if (yaml_parser_parse(parser, &event2))
              {
                if (event2.type == YAML_SCALAR_EVENT)
                  {
                    /* Here parse the name... */
                    if (!strcmp(atoms->atomnames[atoms->iatype[count]],
                                (const char*)event2.data.scalar.value))
                      /* Then the coordinates. */
                      done = _yaml_parser_read_double_array(parser,
                                                            atoms->atomnames[atoms->iatype[count]],
                                                            atoms->fxyz + 3 * count, 3, message);
                    else
                      {
                        set_error("Parser error: force %d is applied on atom '%s' while atom"
                                  " %d is named '%s'.\n", count, (const char*)event2.data.scalar.value,
                                  count, atoms->atomnames[atoms->iatype[count]]);
                        done = -1;
                      }
                  }
                else
                  {
                    set_error("Parser error: atom name awaited.\n");
                    done = -1;
                  }
              }
            else
              {
                /* Error treatment. */
                _yaml_parser_error(parser, message);
                done = -1;
              }
            yaml_event_delete(&event2);
            count += 1;
            break;
          case YAML_SEQUENCE_END_EVENT:
            done = 2;
            break;
          case YAML_SCALAR_EVENT:
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
        _yaml_parser_error(parser, message);
        done = -1;
      }
  
  if (done == 2)
    done = 0;

  return done;
}
static int posinp_yaml_forces(yaml_parser_t *parser, PosinpAtoms *atoms, char **message)
{
  yaml_event_t event;
  int done, count;

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
                done = posinp_yaml_force(parser, atoms, message);
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
              done = _yaml_parser_read_keyword(parser, "Units", funits_keys,
                                               &atoms->funits, POSINP_FORCE_N_UNITS, message);
            else if (!strcmp((const char*)event.data.scalar.value, "Values"))
              done = posinp_yaml_force(parser, atoms, message);
            else if (!strcmp((const char*)event.data.scalar.value, "Fnrm"))
              done = _yaml_parser_read_double(parser, &atoms->fnrm, message);
            else if (!strcmp((const char*)event.data.scalar.value, "MaxVal"))
              done = _yaml_parser_read_double(parser, &atoms->maxval, message);
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
        _yaml_parser_error(parser, message);
        done = -1;
      }

  if (done == 2)
    done = 0;

  return done;
}
static int posinp_yaml_properties(yaml_parser_t *parser, PosinpAtoms *atoms, char **message)
{
  yaml_event_t event;
  int done, count;

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
            if (!strcmp((const char*)event.data.scalar.value, "Converged"))
              done = _yaml_parser_read_keyword(parser, "Converged", bool_keys,
                                               &atoms->converged, 2, message);
            else if (!strcmp((const char*)event.data.scalar.value, "Gnrm_wfn"))
              done = _yaml_parser_read_double(parser, &atoms->gnrm_wfn, message);
            else if (!strncmp((const char*)event.data.scalar.value, "Energy", 6))
              {
                done = _find_units(eunits_keys, (const char*)event.data.scalar.value,
                                   &atoms->eunits, POSINP_ENERG_N_UNITS, message);
                done = (!done)?_yaml_parser_read_double(parser, &atoms->energy, message):done;
              }
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
        _yaml_parser_error(parser, message);
        done = -1;
      }

  if (done == 2)
    done = 0;

  return done;
}
#endif

PosinpList* posinp_yaml_parse(const char *filename, char **message)
{
  PosinpList *list;
#ifdef HAVE_YAML
  PosinpList *tmp, *tmp2;
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
          done = posinp_yaml_cell(&parser, tmp->data, message);
        else if (event.type == YAML_SCALAR_EVENT &&
                 !strcmp((const char*)event.data.scalar.value, "Positions"))
          done = posinp_yaml_position(&parser, tmp->data, message);
        else if (event.type == YAML_SCALAR_EVENT &&
                 !strcmp((const char*)event.data.scalar.value, "Forces"))
          done = posinp_yaml_forces(&parser, tmp->data, message);
        else if (event.type == YAML_SCALAR_EVENT &&
                 !strcmp((const char*)event.data.scalar.value, "Properties"))
          done = posinp_yaml_properties(&parser, tmp->data, message);
        else if (event.type == YAML_SCALAR_EVENT &&
                 !strcmp((const char*)event.data.scalar.value, "Comment"))
          done = _yaml_parser_copy_str(&parser, &tmp->data->comment, message);
        else
          done = (event.type == YAML_STREAM_END_EVENT);

        /* The application is responsible for destroying the event object. */
        yaml_event_delete(&event);
      }
    else
      {
        /* Error treatment. */
        _yaml_parser_error(&parser, message);
        done = -1;
      }

  /* Destroy the Parser object. */
  yaml_parser_delete(&parser);
  fclose(input);
#else
  size_t ln;

  set_error("No YAML support, cannot read file '%s'.\n", filename);
  list = (PosinpList*)0;
#endif

  return list;
}
void FC_FUNC_(f90_posinp_yaml_parse, F90_POSINP_YAML_PARSE)(PosinpList **self,
                                                            const char *filename, unsigned int *ln)
{
  char *name;

  name = malloc(sizeof(char) * (*ln + 1));
  memcpy(name, filename, sizeof(char) * *ln);
  name[*ln] = '\0';

  *self = posinp_yaml_parse(name, NULL);

  free(name);
}
void posinp_yaml_free_list(PosinpList *lst)
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
void FC_FUNC_(f90_posinp_yaml_free_list, F90_POSINP_YAML_FREE_LIST)(PosinpList **self)
{
  posinp_yaml_free_list(*self);
}

void FC_FUNC_(f90_posinp_yaml_get_cell, F90_POSINP_YAML_GET_CELL)(PosinpList **self, unsigned int *i,
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
  else
    {
      angdeg[0] = 90.;
      angdeg[1] = 90.;
      angdeg[2] = 90.;
    }

}
void FC_FUNC_(f90_posinp_yaml_get_dims, F90_POSINP_YAML_GET_DIMS)(PosinpList **self, unsigned int *i,
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
void FC_FUNC_(f90_posinp_yaml_get_atoms, F90_POSINP_YAML_GET_ATOMS)(PosinpList **self, unsigned int *i,
                                                            unsigned int *units, double *rxyz, 
                                                            unsigned int *iatype, unsigned int *ifrztyp,
                                                            int *igspin, int *igchg)
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
      memcpy(igspin,  lst->data->igspin,  sizeof(int) * lst->data->nat);
      memcpy(igchg,   lst->data->igchg,   sizeof(int) * lst->data->nat);
      memcpy(rxyz,    lst->data->rxyz,    sizeof(double) * lst->data->nat * 3);
      for (j = 0; j < lst->data->nat; j++)
	iatype[j] += 1;
    }
}
void FC_FUNC_(f90_posinp_yaml_get_atomname, F90_POSINP_YAML_GET_ATOMNAME)(PosinpList **self, unsigned int *i,
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
void FC_FUNC_(f90_posinp_yaml_has_forces, F90_POSINP_YAML_HAS_FORCES)(PosinpList **self, unsigned int *i,
                                                              unsigned int *has_forces)
{
  PosinpList *lst;
  unsigned int j;

  for (lst = *self, j = 0; j < *i; j++)
    if (lst)
      lst = lst->next;

  if (lst)
    *has_forces = (lst->data->fxyz != (double*)0);
}
void FC_FUNC_(f90_posinp_yaml_get_forces, F90_POSINP_YAML_GET_FORCES)(PosinpList **self, unsigned int *i,
                                                              unsigned int *units, double *fnrm,
                                                              double *maxval, double *fxyz)
{
  PosinpList *lst;
  unsigned int j;

  for (lst = *self, j = 0; j < *i; j++)
    if (lst)
      lst = lst->next;

  if (lst)
    {
      if (lst->data->fxyz)
        memcpy(fxyz, lst->data->fxyz, sizeof(double) * lst->data->nat * 3);
      *units  = lst->data->funits;
      *fnrm   = lst->data->fnrm;
      *maxval = lst->data->maxval;
    }
}
void FC_FUNC_(f90_posinp_yaml_get_properties, F90_POSINP_YAML_GET_PROPERTIES)
     (PosinpList **self, unsigned int *i, unsigned int *eunits, double *energy,
      double *gnrm, int *converged)
{
  PosinpList *lst;
  unsigned int j;

  for (lst = *self, j = 0; j < *i; j++)
    if (lst)
      lst = lst->next;

  if (lst)
    {
      *eunits    = lst->data->eunits;
      *energy    = lst->data->energy;
      *converged = lst->data->converged;
      *gnrm      = lst->data->gnrm_wfn;
    }
}
void FC_FUNC_(f90_posinp_yaml_get_comment, F90_POSINP_YAML_GET_COMMENT)
     (PosinpList **self, unsigned int *i, char *comment, unsigned int *len)
{
  PosinpList *lst;
  unsigned int j, ln;

  for (lst = *self, j = 0; j < *i; j++)
    if (lst)
      lst = lst->next;

  if (lst)
    {
      memset(comment, ' ', sizeof(char) * *len);
      if (lst->data->comment)
        {
          ln = strlen(lst->data->comment);
          memcpy(comment, lst->data->comment, sizeof(char) * ((ln <= (*len))?ln:*len));
        }
    }
}

#ifdef TEST_ME
int main(int argc, const char **argv)
{
  PosinpList *lst, *tmp;

  lst = posinp_yaml_parse(argv[1], NULL);
  for (tmp = lst; tmp; tmp = tmp->next)
    {
      fprintf(stdout, "---\n");
      posinp_atoms_trace(tmp->data);
    }
  posinp_yaml_free_list(lst);

  return 0;
}
#endif
