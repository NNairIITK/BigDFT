#ifndef BINDINGS_SIGNALS_H
#define BINDINGS_SIGNALS_H

#include <config.h>

#ifdef HAVE_GLIB
#include <glib.h>
#include <glib-object.h>
#include <gio/gio.h>
#endif

#include "bigdft.h"

typedef struct BigDFT_Main_
{
#ifdef HAVE_GLIB
  GMainLoop *loop;
#endif

  BigDFT_SignalModes kind;
  /* DBus transport variables. */
#ifdef HAVE_GDBUS
  GDBusConnection *bus;
  GDBusObjectManagerServer *manager;
#endif
  guint busId;

  /* Inet transport variables. */
#ifdef HAVE_GLIB
  GSocket *socket, *recv;
  GSource *source;

  GDestroyNotify destroy;
  gpointer destroyData;
#endif

  BigDFT_Wf *wf;
  BigDFT_LocalFields *denspot;
  BigDFT_Energs *energs;
} BigDFT_Main;

void bigdft_signals_free_main(gpointer self);

#endif
