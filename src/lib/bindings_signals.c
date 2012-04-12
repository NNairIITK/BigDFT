#include <config.h>

#ifdef HAVE_GLIB
#include <glib.h>
#include <glib-object.h>
#include <string.h>

#include "bigdft.h"

#ifdef HAVE_GDBUS
#include "bindings_dbus.h"
#endif

typedef struct BigDFT_Main_
{
  GMainLoop *loop;
#ifdef HAVE_GDBUS
  GDBusConnection *bus;
  GDBusObjectManagerServer *manager;
#endif
  guint busId;

  BigDFT_Wf *wf;
  BigDFT_LocalFields *denspot;
  BigDFT_Energs *energs;
} BigDFT_Main;

static gpointer bigdft_main(gpointer data)
{
  GMainLoop *main = (GMainLoop*)data;
  
  g_main_loop_run(main);

  return (gpointer)0;
}
#endif

#ifdef HAVE_GDBUS

/* Wavefunctions stuffs. */
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

/* Local field stuffs. */
static void onDensityReady(BigDFT_LocalFields *denspot_, guint iter, gpointer data)
{
  BigdftDBusLocalFields *denspot = BIGDFT_DBUS_LOCAL_FIELDS(data);

  bigdft_dbus_local_fields_set_ref_dens_pot_ready
    (denspot, bigdft_dbus_local_fields_get_n_dens_pot_ready(denspot));
  bigdft_dbus_local_fields_emit_dens_pot_ready(denspot, iter, BIGDFT_DENSPOT_DENSITY);
  while (bigdft_dbus_local_fields_get_ref_dens_pot_ready(denspot) > 0)
    g_main_context_iteration(NULL, FALSE);
}
static void onVExtReady(BigDFT_LocalFields *denspot_, gpointer data)
{
  BigdftDBusLocalFields *denspot = BIGDFT_DBUS_LOCAL_FIELDS(data);

  bigdft_dbus_local_fields_set_ref_dens_pot_ready
    (denspot, bigdft_dbus_local_fields_get_n_dens_pot_ready(denspot));
  bigdft_dbus_local_fields_emit_dens_pot_ready(denspot, 0, BIGDFT_DENSPOT_V_EXT);
  while (bigdft_dbus_local_fields_get_ref_dens_pot_ready(denspot) > 0)
    g_main_context_iteration(NULL, FALSE);
}
static gboolean onRegisterDenspotReady(BigdftDBusLocalFields *denspot,
                                       GDBusMethodInvocation *invocation,
                                       gpointer user_data)
{
  bigdft_dbus_local_fields_set_n_dens_pot_ready
    (denspot, bigdft_dbus_local_fields_get_n_dens_pot_ready(denspot) + 1);
  bigdft_dbus_local_fields_complete_register_dens_pot_ready(denspot, invocation);
  return TRUE;
}
static gboolean onUnregisterDenspotReady(BigdftDBusLocalFields *denspot,
                                         GDBusMethodInvocation *invocation,
                                         gpointer user_data)
{
  bigdft_dbus_local_fields_set_n_dens_pot_ready
    (denspot, MAX(1, bigdft_dbus_local_fields_get_n_dens_pot_ready(denspot)) - 1);
  bigdft_dbus_local_fields_complete_unregister_dens_pot_ready(denspot, invocation);
  return TRUE;
}
static gboolean onDoneDenspotReady(BigdftDBusLocalFields *denspot,
                                   GDBusMethodInvocation *invocation,
                                   gpointer user_data)
{
  bigdft_dbus_local_fields_set_ref_dens_pot_ready
    (denspot, MAX(1, bigdft_dbus_local_fields_get_ref_dens_pot_ready(denspot)) - 1);
  bigdft_dbus_local_fields_complete_done_dens_pot_ready(denspot, invocation);
  return TRUE;
}
static gboolean onGetDenspot(BigdftDBusLocalFields *denspot,
                             GDBusMethodInvocation *invocation,
                             BigDFT_DensPotIds kind, gpointer user_data)
{
  BigDFT_LocalFields *denspot_ = BIGDFT_LOCALFIELDS(user_data);
  guint size;
  GVariant *data;

  /* Data are on this processus, so we wrap them without copy. */
  switch (kind)
    {
    case BIGDFT_DENSPOT_DENSITY:
      size = denspot_->glr->ni[0] * denspot_->glr->ni[1] * denspot_->glr->ni[2];
      data = g_variant_new_from_data(G_VARIANT_TYPE("ad"), denspot_->rhov,
                                     sizeof(double) * size,
                                     TRUE, (GDestroyNotify)0, (gpointer)0);
      break;
    case BIGDFT_DENSPOT_V_EXT:
      size = denspot_->glr->ni[0] * denspot_->glr->ni[1] * denspot_->glr->ni[2];
      data = g_variant_new_from_data(G_VARIANT_TYPE("ad"), denspot_->v_ext,
                                     sizeof(double) * size,
                                     TRUE, (GDestroyNotify)0, (gpointer)0);
      break;
    }
  bigdft_dbus_local_fields_complete_get_denspot(denspot, invocation, data, size);

  return TRUE;
}

/* Energies stuffs. */
struct _energs_signal
{
  BigdftDBusEnergs *energs;
  guint iter;
  BigDFT_EnergsIds kind;
};
static gboolean _emit_energ_ready(gpointer data)
{
  struct _energs_signal *dt = (struct _energs_signal*)data;
  
  bigdft_dbus_energs_set_ref_energ_ready
    (dt->energs, bigdft_dbus_energs_get_n_energ_ready(dt->energs));
  bigdft_dbus_energs_emit_energ_ready(dt->energs, dt->iter, dt->kind);
  while (bigdft_dbus_energs_get_ref_energ_ready(dt->energs) > 0)
    g_main_context_iteration(NULL, FALSE);

  return FALSE;
}

static void onEKSReady(BigDFT_Energs *energs_, guint iter, gpointer data)
{
  struct _energs_signal *dt;
  BigdftDBusEnergs *energs = BIGDFT_DBUS_ENERGS(data);

  bigdft_dbus_energs_set_e_ks(energs, energs_->eKS);
  dt = g_malloc(sizeof(struct _energs_signal));
  dt->energs = energs;
  dt->iter = iter;
  dt->kind = BIGDFT_ENERGS_EKS;
  g_idle_add_full(G_PRIORITY_DEFAULT, _emit_energ_ready, dt, g_free);
}
static gboolean onRegisterEnergReady(BigdftDBusEnergs *energs,
                                     GDBusMethodInvocation *invocation,
                                     gpointer user_data)
{
  bigdft_dbus_energs_set_n_energ_ready
    (energs, bigdft_dbus_energs_get_n_energ_ready(energs) + 1);
  bigdft_dbus_energs_complete_register_energ_ready(energs, invocation);
  return TRUE;
}
static gboolean onUnregisterEnergReady(BigdftDBusEnergs *energs,
                                       GDBusMethodInvocation *invocation,
                                       gpointer user_data)
{
  bigdft_dbus_energs_set_n_energ_ready
    (energs, MAX(1, bigdft_dbus_energs_get_n_energ_ready(energs)) - 1);
  bigdft_dbus_energs_complete_unregister_energ_ready(energs, invocation);
  return TRUE;
}
static gboolean onDoneEnergReady(BigdftDBusEnergs *energs,
                                 GDBusMethodInvocation *invocation,
                                 gpointer user_data)
{
  bigdft_dbus_energs_set_ref_energ_ready
    (energs, MAX(1, bigdft_dbus_energs_get_ref_energ_ready(energs)) - 1);
  bigdft_dbus_energs_complete_done_energ_ready(energs, invocation);
  return TRUE;
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

void FC_FUNC_(bigdft_signals_add_wf, BIGDFT_SIGNALS_ADD_WF)(gpointer *self, gpointer *wf_)
{
  BigDFT_Main *main = (BigDFT_Main*)(*self);
  BigdftDBusObjectSkeleton *obj;
  BigdftDBusWf *wf;

  main->wf = BIGDFT_WF(*wf_);
  g_object_ref(G_OBJECT(*wf_));

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
                   G_CALLBACK(onGetPsiCompress), (gpointer)main->wf);
  g_signal_connect(G_OBJECT(main->wf), "psi-ready",
                   G_CALLBACK(onPsiReady), (gpointer)wf);
  bigdft_dbus_object_skeleton_set_wf(obj, wf);
  g_object_unref(wf);

  g_dbus_object_manager_server_export(main->manager, G_DBUS_OBJECT_SKELETON(obj));
  g_object_unref(obj);
}
void FC_FUNC_(bigdft_signals_rm_wf, BIGDFT_SIGNALS_RM_WF)(gpointer *self)
{
  BigDFT_Main *main = (BigDFT_Main*)(*self);
  g_dbus_object_manager_server_unexport(main->manager, "/outputs/DFT_wavefunctions");
}
void FC_FUNC_(bigdft_signals_add_denspot, BIGDFT_SIGNALS_ADD_DENSPOT)(gpointer *self,
                                                                      gpointer *denspot_)
{
  BigDFT_Main *main = (BigDFT_Main*)(*self);
  BigdftDBusObjectSkeleton *obj;
  BigdftDBusLocalFields *denspot;

  main->denspot = BIGDFT_LOCALFIELDS(*denspot_);
  g_object_ref(G_OBJECT(*denspot_));

  obj = bigdft_dbus_object_skeleton_new("/outputs/DFT_local_fields");
  /* Local fields. */
  denspot = bigdft_dbus_local_fields_skeleton_new();
  bigdft_dbus_local_fields_set_n_dens_pot_ready(denspot, 0);
  bigdft_dbus_local_fields_set_ref_dens_pot_ready(denspot, 0);
  g_signal_connect(G_OBJECT(denspot), "handle-register-denspot-ready",
                   G_CALLBACK(onRegisterDenspotReady), (gpointer)0);
  g_signal_connect(G_OBJECT(denspot), "handle-unregister-denspot-ready",
                   G_CALLBACK(onUnregisterDenspotReady), (gpointer)0);
  g_signal_connect(G_OBJECT(denspot), "handle-done-denspot-ready",
                   G_CALLBACK(onDoneDenspotReady), (gpointer)0);
  g_signal_connect(G_OBJECT(denspot), "handle-get-denspot",
                   G_CALLBACK(onGetDenspot), (gpointer)main->denspot);
  g_signal_connect(G_OBJECT(main->denspot), "density-ready",
                   G_CALLBACK(onDensityReady), (gpointer)denspot);
  g_signal_connect(G_OBJECT(main->denspot), "v-ext-ready",
                   G_CALLBACK(onVExtReady), (gpointer)denspot);
  bigdft_dbus_object_skeleton_set_local_fields(obj, denspot);
  g_object_unref(denspot);

  g_dbus_object_manager_server_export(main->manager, G_DBUS_OBJECT_SKELETON(obj));
  g_object_unref(obj);
}
void FC_FUNC_(bigdft_signals_rm_denspot, BIGDFT_SIGNALS_RM_DENSPOT)(gpointer *self)
{
  BigDFT_Main *main = (BigDFT_Main*)(*self);
  g_dbus_object_manager_server_unexport(main->manager, "/outputs/DFT_local_fields");
}
void FC_FUNC_(bigdft_signals_add_energs, BIGDFT_SIGNALS_ADD_ENERGS)(gpointer *self,
                                                                    gpointer *energs_)
{
  BigDFT_Main *main = (BigDFT_Main*)(*self);
  BigdftDBusObjectSkeleton *obj;
  BigdftDBusEnergs *energs;

  main->energs = BIGDFT_ENERGS(*energs_);
  g_object_ref(G_OBJECT(*energs_));

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
  g_signal_connect(G_OBJECT(main->energs), "eks-ready",
                   G_CALLBACK(onEKSReady), (gpointer)energs);
  bigdft_dbus_object_skeleton_set_energs(obj, energs);
  g_object_unref(energs);

  g_dbus_object_manager_server_export(main->manager, G_DBUS_OBJECT_SKELETON(obj));
  g_object_unref(obj);
}
void FC_FUNC_(bigdft_signals_rm_energs, BIGDFT_SIGNALS_RM_ENERGS)(gpointer *self)
{
  BigDFT_Main *main = (BigDFT_Main*)(*self);
  g_dbus_object_manager_server_unexport(main->manager, "/outputs/DFT_energies");
}

void FC_FUNC_(bigdft_signals_init, BIGDFT_SIGNALS_INIT)(gpointer *self)
{
#ifdef G_THREADS_ENABLED
  GThread *ld_thread;
  GError *error = (GError*)0;
  BigDFT_Main *main;

  g_type_init();

  main = g_malloc(sizeof(BigDFT_Main));
  main->loop = g_main_loop_new(NULL, FALSE);
  ld_thread = g_thread_create(bigdft_main, (gpointer)main->loop, FALSE, &error);
  main->wf = NULL;
  main->denspot = NULL;
  main->energs = NULL;

#ifdef HAVE_GDBUS
  main->bus = g_bus_get_sync(G_BUS_TYPE_SESSION, NULL, &error);

  main->manager = g_dbus_object_manager_server_new("/outputs");
  g_dbus_object_manager_server_set_connection(main->manager, main->bus);
#endif

  *self = main;
#endif
}

void FC_FUNC_(bigdft_signals_start, BIGDFT_SIGNALS_START)(gpointer *self)
{
#ifdef HAVE_GDBUS
  BigDFT_Main *main = (BigDFT_Main*)(*self);;

  main->busId = g_bus_own_name_on_connection
    (main->bus, "eu.etsf.BigDFT",
     G_BUS_NAME_OWNER_FLAGS_ALLOW_REPLACEMENT | G_BUS_NAME_OWNER_FLAGS_REPLACE,
     on_name_acquired, on_name_lost, NULL, NULL);
#endif
}

void FC_FUNC_(bigdft_signals_stop, BIGDFT_SIGNALS_STOP)(gpointer *self)
{
#ifdef HAVE_GDBUS
  BigDFT_Main *main = (BigDFT_Main*)(*self);

  g_bus_unown_name(main->busId);
#endif
}

void FC_FUNC_(bigdft_signals_free, BIGDFT_SIGNALS_FREE)(gpointer *self)
{
#ifdef G_THREADS_ENABLED
  BigDFT_Main *main = (BigDFT_Main*)(*self);

#ifdef HAVE_GDBUS
  g_object_unref(main->manager);
  g_object_unref(main->bus);
#endif

  g_main_loop_quit(main->loop);
  g_main_loop_unref(main->loop);
  if (main->wf)
    g_object_unref(main->wf);
  if (main->denspot)
    g_object_unref(main->denspot);
  if (main->energs)
    g_object_unref(main->energs);
  g_free(main);
#endif
}
