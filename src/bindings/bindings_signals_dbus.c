#include <config.h>

#ifdef HAVE_GDBUS
#include "bindings_signals.h"

#include <string.h>

#include "bindings.h"
#include "bindings_api.h"
#include "bindings_dbus.h"

/* Wavefunctions stuffs. */
void onPsiReady(BigDFT_Wf *wf_, guint iter, gpointer data)
{
  BigdftDBusWf *wf = BIGDFT_DBUS_WF(data);

  bigdft_dbus_wf_set_ref_psi_ready(wf, bigdft_dbus_wf_get_n_psi_ready(wf));
  bigdft_dbus_wf_emit_psi_ready(wf, iter);
  while (bigdft_dbus_wf_get_ref_psi_ready(wf) > 0)
    g_main_context_iteration(NULL, FALSE);
}
gboolean onRegisterPsiReady(BigdftDBusWf *wf, GDBusMethodInvocation *invocation,
                                   gpointer user_data)
{
  bigdft_dbus_wf_set_n_psi_ready(wf, bigdft_dbus_wf_get_n_psi_ready(wf) + 1);
  bigdft_dbus_wf_complete_register_psi_ready(wf, invocation);
  return TRUE;
}
gboolean onUnregisterPsiReady(BigdftDBusWf *wf, GDBusMethodInvocation *invocation,
                                     gpointer user_data)
{
  bigdft_dbus_wf_set_n_psi_ready(wf, MAX(1, bigdft_dbus_wf_get_n_psi_ready(wf)) - 1);
  bigdft_dbus_wf_complete_unregister_psi_ready(wf, invocation);
  return TRUE;
}
gboolean onDonePsiReady(BigdftDBusWf *wf, GDBusMethodInvocation *invocation,
                               gpointer user_data)
{
  bigdft_dbus_wf_set_ref_psi_ready(wf, MAX(1, bigdft_dbus_wf_get_ref_psi_ready(wf)) - 1);
  bigdft_dbus_wf_complete_done_psi_ready(wf, invocation);
  return TRUE;
}
gboolean onGetPsiCompress(BigdftDBusWf *wf, GDBusMethodInvocation *invocation,
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
void onDensityReady(BigDFT_LocalFields *denspot_, guint iter, gpointer data)
{
  BigdftDBusLocalFields *denspot = BIGDFT_DBUS_LOCAL_FIELDS(data);

  bigdft_dbus_local_fields_set_ref_dens_pot_ready
    (denspot, bigdft_dbus_local_fields_get_n_dens_pot_ready(denspot));
  bigdft_dbus_local_fields_emit_dens_pot_ready(denspot, iter, BIGDFT_DENSPOT_DENSITY);
  while (bigdft_dbus_local_fields_get_ref_dens_pot_ready(denspot) > 0)
    g_main_context_iteration(NULL, FALSE);
}
void onVExtReady(BigDFT_LocalFields *denspot_, gpointer data)
{
  BigdftDBusLocalFields *denspot = BIGDFT_DBUS_LOCAL_FIELDS(data);

  bigdft_dbus_local_fields_set_ref_dens_pot_ready
    (denspot, bigdft_dbus_local_fields_get_n_dens_pot_ready(denspot));
  bigdft_dbus_local_fields_emit_dens_pot_ready(denspot, 0, BIGDFT_DENSPOT_V_EXT);
  while (bigdft_dbus_local_fields_get_ref_dens_pot_ready(denspot) > 0)
    g_main_context_iteration(NULL, FALSE);
}
gboolean onRegisterDenspotReady(BigdftDBusLocalFields *denspot,
                                       GDBusMethodInvocation *invocation,
                                       gpointer user_data)
{
  bigdft_dbus_local_fields_set_n_dens_pot_ready
    (denspot, bigdft_dbus_local_fields_get_n_dens_pot_ready(denspot) + 1);
  bigdft_dbus_local_fields_complete_register_dens_pot_ready(denspot, invocation);
  return TRUE;
}
gboolean onUnregisterDenspotReady(BigdftDBusLocalFields *denspot,
                                         GDBusMethodInvocation *invocation,
                                         gpointer user_data)
{
  bigdft_dbus_local_fields_set_n_dens_pot_ready
    (denspot, MAX(1, bigdft_dbus_local_fields_get_n_dens_pot_ready(denspot)) - 1);
  bigdft_dbus_local_fields_complete_unregister_dens_pot_ready(denspot, invocation);
  return TRUE;
}
gboolean onDoneDenspotReady(BigdftDBusLocalFields *denspot,
                                   GDBusMethodInvocation *invocation,
                                   gpointer user_data)
{
  bigdft_dbus_local_fields_set_ref_dens_pot_ready
    (denspot, MAX(1, bigdft_dbus_local_fields_get_ref_dens_pot_ready(denspot)) - 1);
  bigdft_dbus_local_fields_complete_done_dens_pot_ready(denspot, invocation);
  return TRUE;
}
gboolean onGetDenspot(BigdftDBusLocalFields *denspot,
                             GDBusMethodInvocation *invocation,
                             BigDFT_DensPotIds kind, gpointer user_data)
{
  BigDFT_LocalFields *localfields = BIGDFT_LOCALFIELDS(user_data);
  guint size, iproc = 0, new;
  GVariant *data;
  f90_pointer_double tmp;

  /* We do a full_local_potential to gather all the data. */
  F90_1D_POINTER_INIT(&tmp);
  switch (kind)
    {
    case BIGDFT_DENSPOT_DENSITY:
      FC_FUNC_(denspot_full_density, DENSPOT_FULL_DENSITY)(localfields->data, &tmp, (int*)&iproc, (int*)&new);
      break;
    case BIGDFT_DENSPOT_V_EXT:
      FC_FUNC_(denspot_full_v_ext, DENSPOT_FULL_V_EXT)(localfields->data, &tmp, (int*)&iproc, (int*)&new);
      break;
    }
  size = localfields->ni[0] * localfields->ni[1] * localfields->ni[2];
  data = g_variant_new_from_data(G_VARIANT_TYPE("ad"), tmp.data,
                                 sizeof(double) * size,
                                 TRUE, (GDestroyNotify)0, (gpointer)0);
  bigdft_dbus_local_fields_complete_get_denspot(denspot, invocation, data, size);

  FC_FUNC_(deallocate_double_1d, DEALLOCATE_DOUBLE_1D)(&tmp);

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
void onEKSReadyDBus(BigDFT_Energs *energs_, guint iter, gpointer data)
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
gboolean onRegisterEnergReady(BigdftDBusEnergs *energs,
                                     GDBusMethodInvocation *invocation,
                                     gpointer user_data)
{
  bigdft_dbus_energs_set_n_energ_ready
    (energs, bigdft_dbus_energs_get_n_energ_ready(energs) + 1);
  bigdft_dbus_energs_complete_register_energ_ready(energs, invocation);
  return TRUE;
}
gboolean onUnregisterEnergReady(BigdftDBusEnergs *energs,
                                       GDBusMethodInvocation *invocation,
                                       gpointer user_data)
{
  bigdft_dbus_energs_set_n_energ_ready
    (energs, MAX(1, bigdft_dbus_energs_get_n_energ_ready(energs)) - 1);
  bigdft_dbus_energs_complete_unregister_energ_ready(energs, invocation);
  return TRUE;
}
gboolean onDoneEnergReady(BigdftDBusEnergs *energs,
                                 GDBusMethodInvocation *invocation,
                                 gpointer user_data)
{
  bigdft_dbus_energs_set_ref_energ_ready
    (energs, MAX(1, bigdft_dbus_energs_get_ref_energ_ready(energs)) - 1);
  bigdft_dbus_energs_complete_done_energ_ready(energs, invocation);
  return TRUE;
}

void on_name_acquired(GDBusConnection *connection, const gchar *name, gpointer data)
{
  g_print ("Acquired the name %s\n", name);
}
void on_name_lost(GDBusConnection *connection, const gchar *name, gpointer data)
{
  g_print ("Lost the name %s\n", name);
}



#endif
