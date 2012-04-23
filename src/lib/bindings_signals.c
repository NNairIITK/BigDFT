#include <config.h>

#ifdef HAVE_GLIB
#include <glib.h>
#include <glib-object.h>
#include <gio/gio.h>
#include <string.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_signals.h"

#ifdef HAVE_GDBUS
#include "bindings_dbus.h"
#endif

/* Callbacks for Inet transport. */
gboolean onClientConnection(GSocket *socket, GIOCondition condition,
                            gpointer user_data);
void onPsiReadyInet(BigDFT_Wf *wf, guint iter, gpointer *data);
void onEKSReadyInet(BigDFT_Energs *energs, guint iter, gpointer *data);
void onDensityReadyInet(BigDFT_LocalFields *localfields, guint iter, gpointer *data);
void onVExtReadyInet(BigDFT_LocalFields *localfields, gpointer *data);
void onIterHamInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data);
void onIterSubInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data);
void onIterWfnInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data);
void onDoneHamInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data);
void onDoneSubInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data);
void onDoneWfnInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data);


static gpointer bigdft_main(gpointer data)
{
  GMainLoop *loop = (GMainLoop*)data;
  
  g_main_loop_run(loop);

  return (gpointer)0;
}
#endif

void FC_FUNC_(bigdft_signals_add_wf, BIGDFT_SIGNALS_ADD_WF)(gpointer *self, gpointer *wf_)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);
#ifdef HAVE_GDBUS
  BigdftDBusObjectSkeleton *obj;
  BigdftDBusWf *wf;
#endif

  bmain->wf = BIGDFT_WF(*wf_);
  g_object_ref(G_OBJECT(*wf_));

  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      obj = bigdft_dbus_object_skeleton_new("/outputs/DFT_wavefunctions");
      /* Wavefunctions. */
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
                       G_CALLBACK(onGetPsiCompress), (gpointer)bmain->wf);
      g_signal_connect(G_OBJECT(bmain->wf), "psi-ready",
                       G_CALLBACK(onPsiReady), (gpointer)wf);
      bigdft_dbus_object_skeleton_set_wf(obj, wf);
      g_object_unref(wf);

      g_dbus_object_manager_server_export(bmain->manager, G_DBUS_OBJECT_SKELETON(obj));
      g_object_unref(obj);
#else
      g_warning("Signals init: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      bmain->wf_id = g_signal_connect(G_OBJECT(bmain->wf), "psi-ready",
                                     G_CALLBACK(onPsiReadyInet), (gpointer)&bmain->recv);
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }
}
void FC_FUNC_(bigdft_signals_rm_wf, BIGDFT_SIGNALS_RM_WF)(gpointer *self)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);

  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      g_dbus_object_manager_server_unexport(bmain->manager, "/outputs/DFT_wavefunctions");
#else
      g_warning("Signals init: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      g_signal_handler_disconnect(G_OBJECT(bmain->wf), bmain->wf_id);
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }
}
void FC_FUNC_(bigdft_signals_add_denspot, BIGDFT_SIGNALS_ADD_DENSPOT)(gpointer *self,
                                                                      gpointer *denspot_)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);
#ifdef HAVE_GDBUS
  BigdftDBusObjectSkeleton *obj;
  BigdftDBusLocalFields *denspot;
#endif

  bmain->denspot = BIGDFT_LOCALFIELDS(*denspot_);
  g_object_ref(G_OBJECT(*denspot_));

  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      obj = bigdft_dbus_object_skeleton_new("/outputs/DFT_local_fields");
      /* Local fields. */
      denspot = bigdft_dbus_local_fields_skeleton_new();
      bigdft_dbus_local_fields_set_n_dens_pot_ready(denspot, 0);
      bigdft_dbus_local_fields_set_ref_dens_pot_ready(denspot, 0);
      g_signal_connect(G_OBJECT(denspot), "handle-register-dens-pot-ready",
                       G_CALLBACK(onRegisterDenspotReady), (gpointer)0);
      g_signal_connect(G_OBJECT(denspot), "handle-unregister-dens-pot-ready",
                       G_CALLBACK(onUnregisterDenspotReady), (gpointer)0);
      g_signal_connect(G_OBJECT(denspot), "handle-done-dens-pot-ready",
                       G_CALLBACK(onDoneDenspotReady), (gpointer)0);
      g_signal_connect(G_OBJECT(denspot), "handle-get-denspot",
                       G_CALLBACK(onGetDenspot), (gpointer)bmain->denspot);
      g_signal_connect(G_OBJECT(bmain->denspot), "density-ready",
                       G_CALLBACK(onDensityReady), (gpointer)denspot);
      g_signal_connect(G_OBJECT(bmain->denspot), "v-ext-ready",
                       G_CALLBACK(onVExtReady), (gpointer)denspot);
      bigdft_dbus_object_skeleton_set_local_fields(obj, denspot);
      g_object_unref(denspot);

      g_dbus_object_manager_server_export(bmain->manager, G_DBUS_OBJECT_SKELETON(obj));
      g_object_unref(obj);
#else
      g_warning("Signals init: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      bmain->vext_id = g_signal_connect(G_OBJECT(bmain->denspot), "v-ext-ready",
                                       G_CALLBACK(onVExtReadyInet), (gpointer)&bmain->recv);
      bmain->denspot_id = g_signal_connect(G_OBJECT(bmain->denspot), "density-ready",
                                          G_CALLBACK(onDensityReadyInet), (gpointer)&bmain->recv);
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }
}
void FC_FUNC_(bigdft_signals_rm_denspot, BIGDFT_SIGNALS_RM_DENSPOT)(gpointer *self)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);

  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      g_dbus_object_manager_server_unexport(bmain->manager, "/outputs/DFT_local_fields");
#else
      g_warning("Signals init: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      g_signal_handler_disconnect(G_OBJECT(bmain->denspot), bmain->vext_id);
      g_signal_handler_disconnect(G_OBJECT(bmain->denspot), bmain->denspot_id);
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }
}
void FC_FUNC_(bigdft_signals_add_energs, BIGDFT_SIGNALS_ADD_ENERGS)(gpointer *self,
                                                                    gpointer *energs_)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);
#ifdef HAVE_GDBUS
  BigdftDBusObjectSkeleton *obj;
  BigdftDBusEnergs *energs;
#endif

  bmain->energs = BIGDFT_ENERGS(*energs_);
  g_object_ref(G_OBJECT(*energs_));

  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      obj = bigdft_dbus_object_skeleton_new("/outputs/DFT_energies");
      /* Energies. */
      energs = bigdft_dbus_energs_skeleton_new();
      bigdft_dbus_energs_set_n_energ_ready(energs, 0);
      bigdft_dbus_energs_set_ref_energ_ready(energs, 0);
      g_signal_connect(G_OBJECT(energs), "handle-register-energ-ready",
                       G_CALLBACK(onRegisterEnergReady), (gpointer)0);
      g_signal_connect(G_OBJECT(energs), "handle-unregister-energ-ready",
                       G_CALLBACK(onUnregisterEnergReady), (gpointer)0);
      g_signal_connect(G_OBJECT(energs), "handle-done-energ-ready",
                       G_CALLBACK(onDoneEnergReady), (gpointer)0);
      g_signal_connect(G_OBJECT(bmain->energs), "eks-ready",
                       G_CALLBACK(onEKSReadyDBus), (gpointer)energs);
      bigdft_dbus_object_skeleton_set_energs(obj, energs);
      g_object_unref(energs);

      g_dbus_object_manager_server_export(bmain->manager, G_DBUS_OBJECT_SKELETON(obj));
      g_object_unref(obj);
#else
      g_warning("Signals init: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      bmain->energs_id = g_signal_connect(G_OBJECT(bmain->energs), "eks-ready",
                                         G_CALLBACK(onEKSReadyInet), (gpointer)&bmain->recv);
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }
}
void FC_FUNC_(bigdft_signals_rm_energs, BIGDFT_SIGNALS_RM_ENERGS)(gpointer *self)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);

  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      g_dbus_object_manager_server_unexport(bmain->manager, "/outputs/DFT_energies");
#else
      g_warning("Signals init: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      g_signal_handler_disconnect(G_OBJECT(bmain->energs), bmain->energs_id);
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }
}

void FC_FUNC_(bigdft_signals_add_optloop, BIGDFT_SIGNALS_ADD_OPTLOOP)(gpointer *self,
                                                                      gpointer *optloop_)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);

  bmain->optloop = BIGDFT_OPTLOOP(*optloop_);
  g_object_ref(G_OBJECT(*optloop_));

  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      g_warning("Signals init: DBus transport not implemented for optloop.");
#else
      g_warning("Signals init: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      bmain->optloop_iter_ham_id = g_signal_connect(G_OBJECT(bmain->optloop), "iter-hamiltonian",
                                                    G_CALLBACK(onIterHamInet), (gpointer)&bmain->recv);
      bmain->optloop_iter_sub_id = g_signal_connect(G_OBJECT(bmain->optloop), "iter-subspace",
                                                    G_CALLBACK(onIterSubInet), (gpointer)&bmain->recv);
      bmain->optloop_iter_wfn_id = g_signal_connect(G_OBJECT(bmain->optloop), "iter-wavefunctions",
                                                    G_CALLBACK(onIterWfnInet), (gpointer)&bmain->recv);
      bmain->optloop_done_ham_id = g_signal_connect(G_OBJECT(bmain->optloop), "done-hamiltonian",
                                                    G_CALLBACK(onDoneHamInet), (gpointer)&bmain->recv);
      bmain->optloop_done_sub_id = g_signal_connect(G_OBJECT(bmain->optloop), "done-subspace",
                                                    G_CALLBACK(onDoneSubInet), (gpointer)&bmain->recv);
      bmain->optloop_done_wfn_id = g_signal_connect(G_OBJECT(bmain->optloop), "done-wavefunctions",
                                                    G_CALLBACK(onDoneWfnInet), (gpointer)&bmain->recv);
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }
}
void FC_FUNC_(bigdft_signals_rm_optloop, BIGDFT_SIGNALS_RM_OPTLOOP)(gpointer *self)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);

  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      g_dbus_object_manager_server_unexport(bmain->manager, "/outputs/DFT_energies");
#else
      g_warning("Signals init: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      g_signal_handler_disconnect(G_OBJECT(bmain->optloop), bmain->optloop_iter_ham_id);
      g_signal_handler_disconnect(G_OBJECT(bmain->optloop), bmain->optloop_iter_sub_id);
      g_signal_handler_disconnect(G_OBJECT(bmain->optloop), bmain->optloop_iter_wfn_id);
      g_signal_handler_disconnect(G_OBJECT(bmain->optloop), bmain->optloop_done_ham_id);
      g_signal_handler_disconnect(G_OBJECT(bmain->optloop), bmain->optloop_done_sub_id);
      g_signal_handler_disconnect(G_OBJECT(bmain->optloop), bmain->optloop_done_wfn_id);
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }
}

void FC_FUNC_(bigdft_signals_init, BIGDFT_SIGNALS_INIT)(gpointer *self, guint *kind,
                                                        gchar *domain, guint *ln)
{
  BigDFT_Main *bmain;
#ifdef G_THREADS_ENABLED
  GThread *ld_thread;
#endif
  GError *error;
#ifdef HAVE_GLIB
  GResolver *dns;
  GList *lst, *tmp;
  GSocketAddress *sockaddr;
  gboolean bind;
  gchar *fqdn, *dom;
#endif

#ifdef HAVE_GLIB
  g_type_init();
#endif

  error = (GError*)0;

  bmain = g_malloc0(sizeof(BigDFT_Main));

#ifdef G_THREADS_ENABLED
  bmain->loop = g_main_loop_new(NULL, FALSE);
  ld_thread = g_thread_create(bigdft_main, (gpointer)bmain->loop, FALSE, &error);
#endif
  bmain->wf = NULL;
  bmain->denspot = NULL;
  bmain->energs = NULL;

  bmain->kind = *kind;
  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      bmain->bus = g_bus_get_sync(G_BUS_TYPE_SESSION, NULL, &error);

      bmain->manager = g_dbus_object_manager_server_new("/outputs");
      g_dbus_object_manager_server_set_connection(bmain->manager, bmain->bus);
#else
      g_warning("Signals init: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      bmain->socket = g_socket_new(G_SOCKET_FAMILY_IPV4, G_SOCKET_TYPE_STREAM,
                                  G_SOCKET_PROTOCOL_DEFAULT, &error);
      if (!bmain->socket)
        {
          g_warning("%s", error->message);
          g_error_free(error);
          bmain->kind = BIGDFT_SIGNALS_NONE;
          break;
        }
      g_socket_set_blocking(bmain->socket, FALSE);
      bind = FALSE;
      dns = g_resolver_get_default();
      if (*ln > 0)
        {
          dom = g_strndup(domain, *ln);
          fqdn = g_strdup_printf("%s.%s", g_get_host_name(), dom);
          g_free(dom);
        }
      else
        fqdn = g_strdup(g_get_host_name());
      g_print(" |  Create a socket for hostname '%s'.\n", fqdn);
      lst = g_resolver_lookup_by_name(dns, fqdn, NULL, &error);
      g_object_unref(dns);
      g_free(fqdn);
      for (tmp = lst; tmp && !bind; tmp = g_list_next(tmp))
        if (!tmp->next ||
            g_inet_address_to_bytes((GInetAddress*)tmp->data)[0] != (guint8)127)
          {
            sockaddr = g_inet_socket_address_new((GInetAddress*)tmp->data, (guint16)91691);
            bind = g_socket_bind(bmain->socket, sockaddr, TRUE, &error);
            g_print(" |   try to bind to '%s' -> %d.\n",
                    g_inet_address_to_string((GInetAddress*)tmp->data), bind);
            if (!bind)
              {
                g_warning("%s", error->message);
                g_error_free(error);
                error = (GError*)0;
              }
            g_object_unref(sockaddr);
          }
      g_resolver_free_addresses(lst);
      if (!bind)
        {
          bmain->kind = BIGDFT_SIGNALS_NONE;
          break;
        }
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }

  *self = bmain;
}

static gboolean onClientTimeout(gpointer data)
{
  g_cancellable_cancel((GCancellable*)data);
  g_object_unref(data);

  return FALSE;
}

void FC_FUNC_(bigdft_signals_start, BIGDFT_SIGNALS_START)(gpointer *self, int *timeout)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);
  GError *error;
  GCancellable *cancellable;

  error = (GError*)0;
  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      bmain->busId = g_bus_own_name_on_connection
        (bmain->bus, "eu.etsf.BigDFT",
         G_BUS_NAME_OWNER_FLAGS_ALLOW_REPLACEMENT | G_BUS_NAME_OWNER_FLAGS_REPLACE,
         on_name_acquired, on_name_lost, NULL, NULL);
#else
      g_warning("Signals start: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      /* g_print("Make the socket listen to one client max.\n"); */
      g_socket_set_listen_backlog(bmain->socket, 1);
      if (!g_socket_listen(bmain->socket, &error))
        {
          g_warning("%s", error->message);
          g_error_free(error);
          bmain->kind = BIGDFT_SIGNALS_NONE;
          break;
        }
      if (*timeout != 0)
        {
          if (*timeout > 0)
            {
              cancellable = g_cancellable_new();
              g_object_ref(cancellable);
              g_timeout_add(*timeout * 1000, onClientTimeout, cancellable);
            }
          else
            cancellable = (GCancellable*)0;
          if (g_socket_condition_wait(bmain->socket, G_IO_IN, cancellable, &error))
            onClientConnection(bmain->socket, G_IO_IN, bmain);
          else if (error->code != G_IO_ERROR_CANCELLED)
            {
              g_warning("%s", error->message);
              g_error_free(error);
              bmain->kind = BIGDFT_SIGNALS_NONE;
            }
          else
            g_error_free(error);
          if (cancellable)
            g_object_unref(cancellable);
        }
      if (bmain->kind == BIGDFT_SIGNALS_INET)
        {
          bmain->source = g_socket_create_source(bmain->socket, G_IO_IN | G_IO_HUP, NULL);
          g_source_set_callback(bmain->source, (GSourceFunc)onClientConnection, bmain, NULL);
          g_source_attach(bmain->source, NULL);
        }
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }
}

void FC_FUNC_(bigdft_signals_stop, BIGDFT_SIGNALS_STOP)(gpointer *self)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);
  GError *error;

  error = (GError*)0;
  switch (bmain->kind)
    {
    case BIGDFT_SIGNALS_DBUS:
#ifdef HAVE_GDBUS
      g_bus_unown_name(bmain->busId);
#else
      g_warning("Signals stop: DBus transport unavailable.");
#endif
      break;
    case BIGDFT_SIGNALS_INET:
#ifdef HAVE_GLIB
      /* g_print("Close the socket.\n"); */
      if (!g_socket_close(bmain->socket, &error))
        {
          g_warning("%s", error->message);
          g_error_free(error);
          bmain->kind = BIGDFT_SIGNALS_NONE;
        }
#else
      g_warning("Signals init: Inet transport unavailable.");
#endif
      break;
    default:
      break;
    }
}

void bigdft_signals_free_main(gpointer self)
{
#ifdef G_THREADS_ENABLED
  BigDFT_Main *bmain = (BigDFT_Main*)(self);

#ifdef HAVE_GDBUS
  if (bmain->manager)
    g_object_unref(bmain->manager);
  if (bmain->bus)
    g_object_unref(bmain->bus);
#endif

#ifdef HAVE_GLIB
  if (bmain->socket)
    g_object_unref(bmain->socket);
  if (bmain->source)
    g_source_unref(bmain->source);

  if (bmain->loop)
    {
      g_main_loop_quit(bmain->loop);
      g_main_loop_unref(bmain->loop);
    }
#endif

  if (bmain->wf)
    g_object_unref(bmain->wf);
  if (bmain->denspot)
    g_object_unref(bmain->denspot);
  if (bmain->energs)
    g_object_unref(bmain->energs);
  if (bmain->optloop)
    {
      if (bmain->optloop_sync)
        g_signal_handler_disconnect(G_OBJECT(bmain->optloop), bmain->optloop_sync);
      g_object_unref(bmain->optloop);
    }

  g_free(bmain);
#endif
}

void FC_FUNC_(bigdft_signals_free, BIGDFT_SIGNALS_FREE)(gpointer *self)
{
  BigDFT_Main *bmain = (BigDFT_Main*)(*self);

  bigdft_signals_free_main(bmain);
}
