#include "glib.h"
#include "glib-object.h"

#include "bigdft.h"
#include "bindings_dbus.h"

static void print_obj(GDBusObjectManager *manager)
{
  GList *objects, *l;
  BigdftDBusObject *object;
  GList *interfaces, *ll;
  GDBusInterface *interface;

  g_print ("Object manager at %s\n", g_dbus_object_manager_get_object_path(manager));
  objects = g_dbus_object_manager_get_objects (manager);
  for (l = objects; l != NULL; l = l->next)
    {
      object = BIGDFT_DBUS_OBJECT(l->data);
      g_print (" - Object at %s\n", g_dbus_object_get_object_path (G_DBUS_OBJECT (object)));
      
      interfaces = g_dbus_object_get_interfaces(G_DBUS_OBJECT (object));
      for (ll = interfaces; ll != NULL; ll = ll->next)
        {
          interface = G_DBUS_INTERFACE (ll->data);
          g_print ("   - Interface %s\n", g_dbus_interface_get_info(interface)->name);
        }
      g_list_free_full(interfaces, g_object_unref);
    }
  g_list_free_full(objects, g_object_unref);
}

static void onPsiReady(BigdftDBusWf *wf, guint iter, gpointer data)
{
  GError *error;
  GVariant *psi;
  guint size, i, n;
  const double *psi_;
  BigDFT_Wf *wf_ = BIGDFT_WF(data);
  double *psir, *psii;
  double minDens, maxDens;

  g_printf("Callback for 'psi-ready' signal at iter %d.\n", iter);

  /* Pulling wavefunction 4 as a test. */
  error = (GError*)0;
  if (!bigdft_dbus_wf_call_get_psi_compress_sync(wf, 1, 4, BIGDFT_SPIN_UP, 
                                                 BIGDFT_PARTIAL_DENSITY,
                                                 &psi, &size, NULL, &error))
    {
      g_warning("%s", error->message);
      g_error_free(error);
      return;
    }
  else
    g_variant_ref(psi);

  /* Finishing signal, raising lock. */
  error = (GError*)0;
  if (!bigdft_dbus_wf_call_done_psi_ready_sync(wf, NULL, &error))
    {
      g_warning("%s", error->message);
      g_error_free(error);
    }

  psi_ = (const double*)g_variant_get_data(psi);
  g_print("%g %d\n", psi_[0], size * (guint)sizeof(double));
  psir = bigdft_locreg_convert_to_isf(BIGDFT_LOCREG(wf_->lzd), psi_);
  if (BIGDFT_ORBS(wf_)->nspinor == 2)
    psii = bigdft_locreg_convert_to_isf(BIGDFT_LOCREG(wf_->lzd), psi_ + size / 2);
  g_variant_unref(psi);

  minDens = G_MAXDOUBLE;
  maxDens = 0.;
  n = BIGDFT_LOCREG(wf_->lzd)->ni[0] * 
    BIGDFT_LOCREG(wf_->lzd)->ni[1] * 
    BIGDFT_LOCREG(wf_->lzd)->ni[2];
  g_print("%d %d %d -> %d\n", BIGDFT_LOCREG(wf_->lzd)->ni[0],
          BIGDFT_LOCREG(wf_->lzd)->ni[1], BIGDFT_LOCREG(wf_->lzd)->ni[2], n);
  for (i = 0; i < n; i++)
    {
      psir[i] *= psir[i];
      if (BIGDFT_ORBS(wf_)->nspinor == 2)
        psir[i] += psii[i] * psii[i];
      minDens = MIN(minDens, psir[i]);
      maxDens = MAX(maxDens, psir[i]);
    }
  g_print(" Band 4 has min partial density %g and max %g.\n", minDens, maxDens);

  g_free(psir);
  if (BIGDFT_ORBS(wf_)->nspinor == 2)
    g_free(psii);
}

static void on_object_added(GDBusObjectManager *manager, GDBusObject *object,
                            gpointer user_data)
{
  GList *interfaces, *ll;
  GDBusInterface *interface;
  GError *error;

  interfaces = g_dbus_object_get_interfaces(G_DBUS_OBJECT(object));
  for (ll = interfaces; ll != NULL; ll = ll->next)
    {
      interface = G_DBUS_INTERFACE(ll->data);
      g_print ("   - Interface %s\n", g_dbus_interface_get_info(interface)->name);
      if (BIGDFT_DBUS_IS_WF(interface))
        {
          error = (GError*)0;
          if (!bigdft_dbus_wf_call_register_psi_ready_sync(BIGDFT_DBUS_WF(interface),
                                                           NULL, &error))
            {
              g_warning("%s", error->message);
              g_error_free(error);
            }
          g_signal_connect(interface, "psi-ready",
                           G_CALLBACK(onPsiReady), user_data);
        }
    }
  g_list_free_full(interfaces, g_object_unref);
}

int main(int argc, const char **argv)
{
  GMainLoop *loop;
  GDBusObjectManager *manager;
  GError *error;
  BigDFT_Inputs *in;
  BigDFT_Wf *wf;
  double *radii;

  g_type_init ();

  /* Load test BigDFT run. */
  in = bigdft_inputs_new("test");
  wf = bigdft_wf_new();
  bigdft_atoms_set_structure_from_file(BIGDFT_ATOMS(wf->lzd), "test.ascii");
  bigdft_atoms_set_symmetries(BIGDFT_ATOMS(wf->lzd), !in->disableSym, -1., in->elecfield);
  bigdft_inputs_parse_additional(in, BIGDFT_ATOMS(wf->lzd));
  bigdft_atoms_set_psp(BIGDFT_ATOMS(wf->lzd), in->ixc, in->nspin, (const gchar*)0);
  radii = bigdft_atoms_get_radii(BIGDFT_ATOMS(wf->lzd), in->crmult, in->frmult, 0.);
  bigdft_locreg_set_radii(BIGDFT_LOCREG(wf->lzd), radii);
  g_free(radii);
  bigdft_lzd_set_size(wf->lzd, in->h, in->crmult, in->frmult);
  bigdft_locreg_set_wave_descriptors(BIGDFT_LOCREG(wf->lzd));
  bigdft_wf_define(wf, in, 0, 1);

  loop = g_main_loop_new (NULL, FALSE);

  error = NULL;
  manager = bigdft_dbus_object_manager_client_new_for_bus_sync
    (G_BUS_TYPE_SESSION, G_DBUS_OBJECT_MANAGER_CLIENT_FLAGS_NONE,
     "eu.etsf.BigDFT", "/outputs",
     NULL, &error);

  g_signal_connect(manager, "object-added",
                   G_CALLBACK(on_object_added), (gpointer)wf);

  g_main_loop_run (loop);

  g_object_unref(in);
  g_object_unref(wf);
  g_object_unref(manager);
  g_main_loop_unref(loop);

  return 0;
}
