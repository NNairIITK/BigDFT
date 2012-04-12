#include <config.h>

#ifdef HAVE_GLIB
#include <glib.h>
#include <glib-object.h>
#include <string.h>

#include "bigdft.h"

typedef struct BigDFT_Main_
{
  GMainLoop *loop;
  guint bus;

  BigDFT_Wf *wf;
} BigDFT_Main;

static gpointer bigdft_main(gpointer data)
{
  GMainLoop *main = (GMainLoop*)data;
  
  g_main_loop_run (main);

  return (gpointer)0;
}
#endif

#ifdef HAVE_GDBUS
#include "bindings_dbus.h"

static void onPsiReady(BigDFT_Wf *wf_, guint iter, gpointer data)
{
  guint i, j;
  BigdftDBusWf *wf = BIGDFT_DBUS_WF(data);

  bigdft_dbus_wf_set_ref_psi_ready(wf, bigdft_dbus_wf_get_n_psi_ready(wf));
  bigdft_dbus_wf_emit_psi_ready(wf, iter);
  while (bigdft_dbus_wf_get_ref_psi_ready(wf) > 0)
    g_main_context_iteration(NULL, FALSE);
}
static gboolean onRegisterPsiReady(BigdftDBusWf *wf, GDBusMethodInvocation *invocation,
                                   gpointer user_data)
{
  bigdft_dbus_wf_set_n_psi_ready(wf, bigdft_dbus_wf_get_n_psi_ready(wf) + 1);
  bigdft_dbus_wf_complete_register_psi_ready(wf, invocation);
  return TRUE;
}
static gboolean onUnregisterPsiReady(BigdftDBusWf *wf, GDBusMethodInvocation *invocation,
                                     gpointer user_data)
{
  bigdft_dbus_wf_set_n_psi_ready(wf, MAX(1, bigdft_dbus_wf_get_n_psi_ready(wf)) - 1);
  bigdft_dbus_wf_complete_unregister_psi_ready(wf, invocation);
  return TRUE;
}
static gboolean onDonePsiReady(BigdftDBusWf *wf, GDBusMethodInvocation *invocation,
                               gpointer user_data)
{
  bigdft_dbus_wf_set_ref_psi_ready(wf, MAX(1, bigdft_dbus_wf_get_ref_psi_ready(wf)) - 1);
  bigdft_dbus_wf_complete_done_psi_ready(wf, invocation);
  return TRUE;
}
static gboolean onGetPsiCompress(BigdftDBusWf *wf, GDBusMethodInvocation *invocation,
                                 guint ikpt, guint iorb, guint ispin, guint ispinor,
                                 gpointer user_data)
{
  BigDFT_Wf *wf_ = BIGDFT_WF(user_data);
  /* GArray *arr; */
  const double *psic;
  guint size;
  GVariant *psi;
  double *data;

  psic = bigdft_wf_get_psi_compress(wf_, ikpt, iorb, ispin, ispinor, &size, 0);
  if (size == 0)
    {
      g_dbus_method_invocation_return_dbus_error(invocation,
                                                 "eu.etsf.BigDFT.Error.noPsi",
                                                 "Required psi is not available");
      return TRUE;
    }
  if (psic)
    {
      /* Data are on this processus, so we wrap them without copy. */
      /* arr = g_array_sized_new(FALSE, FALSE, sizeof(double), size); */
      /* memcpy(arr->data, psic, sizeof(double) * size); */
      /* arr = g_array_set_size(arr, size); */
      /* psi = g_variant_new_from_data(G_VARIANT_TYPE("ad"), arr->data, sizeof(double) * size, */
      /*                               TRUE, (GDestroyNotify)g_array_unref, arr); */
      psi = g_variant_new_from_data(G_VARIANT_TYPE("ad"), psic, sizeof(double) * size,
                                    TRUE, (GDestroyNotify)0, (gpointer)0);
    }
  else
    {
      /* Data are on processus iproc, so we copy them. */
      data = g_malloc(sizeof(double) * size);
      if (!bigdft_wf_copy_psi_compress(wf_, ikpt, iorb, ispin, ispinor, 0, data, size))
        {
          g_free(data);
          g_dbus_method_invocation_return_dbus_error(invocation,
                                                     "eu.etsf.BigDFT.Error.wrongPsi",
                                                     "Required psi cannot be copied");
          return TRUE;
        }
      psi = g_variant_new_from_data(G_VARIANT_TYPE("ad"), data, sizeof(double) * size,
                                    TRUE, (GDestroyNotify)g_free, data);
    }
  bigdft_dbus_wf_complete_get_psi_compress(wf, invocation, psi, size);

  return TRUE;
}

static GDBusObjectManagerServer *manager = NULL;
static void on_bus_acquired(GDBusConnection *connection, const gchar *name, gpointer data)
{
  BigdftDBusObjectSkeleton *obj;
  BigdftDBusWf *wf;
  BigDFT_Main *main = (BigDFT_Main*)data;

  g_print ("Acquired a message bus connection %s\n", name);

  manager = g_dbus_object_manager_server_new("/outputs");
  
  obj = bigdft_dbus_object_skeleton_new("/outputs/DFT_wavefunctions");
  wf = bigdft_dbus_wf_skeleton_new();
  bigdft_dbus_wf_set_n_psi_ready(wf, 0);
  bigdft_dbus_wf_set_ref_psi_ready(wf, 0);
  g_signal_connect(G_OBJECT(wf), "handle-register-psi-ready",
                   G_CALLBACK(onRegisterPsiReady), (gpointer)0);
  g_signal_connect(G_OBJECT(wf), "handle-unregister-psi-ready",
                   G_CALLBACK(onUnregisterPsiReady), (gpointer)0);
  g_signal_connect(G_OBJECT(wf), "handle-done-psi-ready",
                   G_CALLBACK(onDonePsiReady), (gpointer)0);
  g_signal_connect(G_OBJECT(wf), "handle-get-psi-compress",
                   G_CALLBACK(onGetPsiCompress), (gpointer)main->wf);
  g_signal_connect(G_OBJECT(main->wf), "psi-ready",
                   G_CALLBACK(onPsiReady), (gpointer)wf);
  bigdft_dbus_object_skeleton_set_wf(obj, wf);
  g_object_unref(wf);

  g_dbus_object_manager_server_export(manager, G_DBUS_OBJECT_SKELETON(obj));
  g_object_unref(obj);

  /* Export all objects */
  g_dbus_object_manager_server_set_connection(manager, connection);
}
static void on_name_acquired(GDBusConnection *connection, const gchar *name, gpointer data)
{
  g_print ("Acquired the name %s\n", name);
}
static void on_name_lost(GDBusConnection *connection, const gchar *name, gpointer data)
{
  g_print ("Lost the name %s\n", name);
}
#endif

void FC_FUNC_(bigdft_signals_init, BIGDFT_SIGNALS_INIT)()
{
#ifdef HAVE_GLIB
  g_type_init();
#endif
}

void FC_FUNC_(bigdft_signals_start, BIGDFT_SIGNALS_START)(gpointer *self, gpointer *wf)
{
#ifdef G_THREADS_ENABLED
  GThread *ld_thread;
  GError *error = (GError*)0;
  BigDFT_Main *main;

  main = g_malloc(sizeof(BigDFT_Main));
  main->loop = g_main_loop_new(NULL, FALSE);
  ld_thread = g_thread_create(bigdft_main, (gpointer)main->loop, FALSE, &error);
  main->wf = BIGDFT_WF(*wf);
  g_object_ref(G_OBJECT(*wf));

#ifdef HAVE_GDBUS
  main->bus = g_bus_own_name(G_BUS_TYPE_SESSION, "eu.etsf.BigDFT",
                             G_BUS_NAME_OWNER_FLAGS_ALLOW_REPLACEMENT |
                             G_BUS_NAME_OWNER_FLAGS_REPLACE,
                             on_bus_acquired, on_name_acquired, on_name_lost,
                             main, NULL);
#endif

  *self = main;
#endif
}

void FC_FUNC_(bigdft_signals_stop, BIGDFT_SIGNALS_STOP)(gpointer *self)
{
#ifdef G_THREADS_ENABLED
  BigDFT_Main *main;
  main = (BigDFT_Main*)(*self);
  g_main_loop_quit(main->loop);
#ifdef HAVE_GDBUS
  g_bus_unown_name(main->bus);
  g_object_unref(manager);
#endif
  g_main_loop_unref(main->loop);
  g_object_unref(G_OBJECT(main->wf));
  g_free(main);
#endif
}
