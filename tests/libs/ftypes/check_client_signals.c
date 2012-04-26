#include "glib.h"
#include "glib-object.h"
#include "gio/gio.h"

#include "bigdft.h"

#ifdef HAVE_GDBUS
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

  /* Finishing signal, raising lock. */
  error = (GError*)0;
  if (!bigdft_dbus_wf_call_done_psi_ready_sync(wf, NULL, &error))
    {
      g_warning("%s", error->message);
      g_error_free(error);
    }

  /* psi_ = (const double*)g_variant_get_data(psi); */
  /* g_print("%g %d\n", psi_[0], size * (guint)sizeof(double)); */
  /* psir = bigdft_locreg_convert_to_isf(BIGDFT_LOCREG(wf_->lzd), psi_); */
  /* if (BIGDFT_ORBS(wf_)->nspinor == 2) */
  /*   psii = bigdft_locreg_convert_to_isf(BIGDFT_LOCREG(wf_->lzd), psi_ + size / 2); */
  g_variant_unref(psi);

  /* minDens = G_MAXDOUBLE; */
  /* maxDens = 0.; */
  /* n = BIGDFT_LOCREG(wf_->lzd)->ni[0] *  */
  /*   BIGDFT_LOCREG(wf_->lzd)->ni[1] *  */
  /*   BIGDFT_LOCREG(wf_->lzd)->ni[2]; */
  /* g_print("%d %d %d -> %d\n", BIGDFT_LOCREG(wf_->lzd)->ni[0], */
  /*         BIGDFT_LOCREG(wf_->lzd)->ni[1], BIGDFT_LOCREG(wf_->lzd)->ni[2], n); */
  /* for (i = 0; i < n; i++) */
  /*   { */
  /*     psir[i] *= psir[i]; */
  /*     if (BIGDFT_ORBS(wf_)->nspinor == 2) */
  /*       psir[i] += psii[i] * psii[i]; */
  /*     minDens = MIN(minDens, psir[i]); */
  /*     maxDens = MAX(maxDens, psir[i]); */
  /*   } */
  /* g_print(" Band 4 has min partial density %g and max %g.\n", minDens, maxDens); */

  /* g_free(psir); */
  /* if (BIGDFT_ORBS(wf_)->nspinor == 2) */
  /*   g_free(psii); */
}
static void onEnergReady(BigdftDBusEnergs *energs, guint iter,
                         BigDFT_EnergsIds kind, gpointer data)
{
  GError *error;

  switch (kind)
    {
    case BIGDFT_ENERGS_EKS:
      g_print("Callback for 'eks-ready' signal at iter %d -> %gHt.\n", iter,
              bigdft_dbus_energs_get_e_ks(energs));
      break;
    }
  /* Finishing signal, raising lock. */
  error = (GError*)0;
  if (!bigdft_dbus_energs_call_done_energ_ready_sync(energs, NULL, &error))
    {
      g_warning("%s", error->message);
      g_error_free(error);
    }
}
static void onDenspotReceived(GObject *obj, GAsyncResult *res, gpointer user_data)
{
  GError *error;
  GVariant *data;
  guint i, size, kind;
  const double *vals;
  double charge;
  BigDFT_LocReg *glr = BIGDFT_LOCREG(BIGDFT_WF(user_data)->lzd);

  error = (GError*)0;
  if (!bigdft_dbus_local_fields_call_get_denspot_finish(BIGDFT_DBUS_LOCAL_FIELDS(obj),
                                                        &data, &size, res, &error))
    {
      g_warning("%s", error->message);
      g_error_free(error);
      return;
    }

  /* Now doing something with data. */
  kind = GPOINTER_TO_INT(g_object_get_data(obj, "kind"));
  g_print(" pulling %d elements of kind %d.\n", size, kind);
  vals = g_variant_get_data(data);
  switch (kind)
    {
    case BIGDFT_DENSPOT_DENSITY:
      for (i = 0, charge = 0.; i < size; i++)
        charge += vals[i];
      charge *= glr->h[0] * glr->h[1] * glr->h[2] * 0.125;
      g_print(" System has %16.16f electrons.\n", charge);
      break;
    case BIGDFT_DENSPOT_V_EXT:
      break;
    }
  g_variant_unref(data);
}
static void onDenspotReady(BigdftDBusLocalFields *denspot, guint iter,
                           BigDFT_DensPotIds kind, gpointer user_data)
{
  GError *error;
  GVariant *data;
  guint i, size/* , kind */;
  const double *vals;
  double charge;
  BigDFT_LocReg *glr = BIGDFT_LOCREG(BIGDFT_WF(user_data)->lzd);

  g_print("Callback for 'dens-pot-ready' signal at iter %d for denspot %d.\n", iter, kind);
  /* Pulling data as a test. */
  error = (GError*)0;
  if (!bigdft_dbus_local_fields_call_get_denspot_sync(denspot, kind, &data, &size,
                                                      NULL, &error))
    {
      g_warning("%s", error->message);
      g_error_free(error);
      return;
    }

  /* Finishing signal, raising lock. */
  error = (GError*)0;
  if (!bigdft_dbus_local_fields_call_done_dens_pot_ready_sync(denspot, NULL, &error))
    {
      g_warning("%s", error->message);
      g_error_free(error);
      return;
    }

  /* Now doing something with data. */
  /* g_print(" pulling %d elements of kind %d.\n", size, kind); */
  /* vals = g_variant_get_data(data); */
  /* switch (kind) */
  /*   { */
  /*   case BIGDFT_DENSPOT_DENSITY: */
  /*     for (i = 0, charge = 0.; i < size; i++) */
  /*       charge += vals[i]; */
  /*     charge *= glr->h[0] * glr->h[1] * glr->h[2] * 0.125; */
  /*     g_print(" System has %16.16f electrons.\n", charge); */
  /*     break; */
  /*   case BIGDFT_DENSPOT_V_EXT: */
  /*     break; */
  /*   } */
  g_variant_unref(data);
  g_mem_profile();
}

static void on_object_added(GDBusObjectManager *manager, GDBusObject *object,
                            gpointer user_data)
{
  GList *interfaces, *ll;
  GDBusInterface *interface;
  GError *error;

  g_print ("- Object at %s\n", g_dbus_object_get_object_path(object));
  interfaces = g_dbus_object_get_interfaces(G_DBUS_OBJECT(object));
  for (ll = interfaces; ll != NULL; ll = ll->next)
    {
      interface = G_DBUS_INTERFACE(ll->data);
      g_print ("   - Interface %s\n", g_dbus_interface_get_info(interface)->name);
      /* if (BIGDFT_DBUS_IS_WF(interface)) */
      /*   { */
      /*     error = (GError*)0; */
      /*     if (!bigdft_dbus_wf_call_register_psi_ready_sync(BIGDFT_DBUS_WF(interface), */
      /*                                                      NULL, &error)) */
      /*       { */
      /*         g_warning("%s", error->message); */
      /*         g_error_free(error); */
      /*       } */
      /*     g_signal_connect(interface, "psi-ready", */
      /*                      G_CALLBACK(onPsiReady), user_data); */
      /*   } */
      if (BIGDFT_DBUS_IS_ENERGS(interface))
        {
          error = (GError*)0;
          if (!bigdft_dbus_energs_call_register_energ_ready_sync
              (BIGDFT_DBUS_ENERGS(interface), NULL, &error))
            {
              g_warning("%s", error->message);
              g_error_free(error);
            }
          g_signal_connect(interface, "energ-ready",
                           G_CALLBACK(onEnergReady), NULL);
        }
      if (BIGDFT_DBUS_IS_LOCAL_FIELDS(interface))
        {
          error = (GError*)0;
          if (!bigdft_dbus_local_fields_call_register_dens_pot_ready_sync
              (BIGDFT_DBUS_LOCAL_FIELDS(interface), NULL, &error))
            {
              g_warning("%s", error->message);
              g_error_free(error);
            }
          g_signal_connect(interface, "dens-pot-ready",
                           G_CALLBACK(onDenspotReady), user_data);
        }
    }
  g_list_free_full(interfaces, g_object_unref);
}
static void on_object_removed(GDBusObjectManager *manager, GDBusObject *object,
                            gpointer user_data)
{
  GList *interfaces, *ll;
  GDBusInterface *interface;

  g_print ("- Object at %s\n", g_dbus_object_get_object_path(object));
  interfaces = g_dbus_object_get_interfaces(G_DBUS_OBJECT(object));
  for (ll = interfaces; ll != NULL; ll = ll->next)
    {
      interface = G_DBUS_INTERFACE(ll->data);
      g_print ("   - Interface %s\n", g_dbus_interface_get_info(interface)->name);
      if (BIGDFT_DBUS_IS_WF(interface))
        {
          g_main_loop_quit((GMainLoop*)user_data);
          g_mem_profile();
        }
    }
  g_list_free_full(interfaces, g_object_unref);
}
#endif

static void onEKSReady(BigDFT_Energs *energs, guint iter, gpointer data)
{
  g_print("Get eKS = %gHt at iter %d.\n", energs->eKS, iter);
}

static void onPsiReady(BigDFT_Wf *wf, guint iter, GArray *psic,
                       guint ikpt, guint iorb, guint ispin, gpointer data)
{
  double *psir, *psii;
  guint i, n;
  double minDens, maxDens, norm;

  g_print("Get one wave (%d,%d,%d) at iter %d.\n", ikpt, iorb, ispin, iter);

  psir = bigdft_locreg_convert_to_isf(BIGDFT_LOCREG(wf->lzd), (double*)psic->data);
  if (BIGDFT_ORBS(wf)->nspinor == 2)
    psii = bigdft_locreg_convert_to_isf(BIGDFT_LOCREG(wf->lzd), (double*)psic->data + psic->len / 2);

  minDens = G_MAXDOUBLE;
  maxDens = 0.;
  n = BIGDFT_LOCREG(wf->lzd)->ni[0] * 
    BIGDFT_LOCREG(wf->lzd)->ni[1] * 
    BIGDFT_LOCREG(wf->lzd)->ni[2];
  norm = 0.;
  for (i = 0; i < n; i++)
    {
      psir[i] *= psir[i];
      if (BIGDFT_ORBS(wf)->nspinor == 2)
        psir[i] += psii[i] * psii[i];
      minDens = MIN(minDens, psir[i]);
      maxDens = MAX(maxDens, psir[i]);
      norm += psir[i];
    }
  g_print(" Band has min partial density %g and max %g (nrm = %g).\n",
          minDens, maxDens, norm);

  g_free(psir);
  if (BIGDFT_ORBS(wf)->nspinor == 2)
    g_free(psii);
}
static void onDensityReady(BigDFT_LocalFields *denspot, guint iter, gpointer data)
{
  double dens;
  guint i;

  g_print("Get density at iter %d.\n", iter);

  dens = 0.;
  for (i = 0; i < denspot->ni[0] * denspot->ni[1] * denspot->ni[2]; i++)
    dens += denspot->rhov[i];
  dens *= denspot->h[0] * denspot->h[1] * denspot->h[2];
  g_print(" Density calculated by C is %16.16f.\n", dens);  
}
static void onDoneWavefunctions(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer data)
{
  g_print("Wavefunctions loop done eKS = %gHt with gnrm %g after %d iters.\n",
          energs->eKS, optloop->gnrm, optloop->iter);

  if (optloop->gnrm > 0.4)
    {
      /* Example of changing the number of sub-space diag. */
      optloop->nrepmax = 2;
      bigdft_optloop_sync_to_fortran(optloop);
    }
}

static void onClosedSocket(gpointer data)
{
  g_main_loop_quit((GMainLoop*)data);
}

int main(int argc, const char **argv)
{
#ifdef HAVE_GDBUS
  GDBusObjectManager *manager;
#endif
  GMainLoop *loop;
  GError *error;
  BigDFT_Inputs *in;
  BigDFT_Wf *wf;
  BigDFT_LocalFields *denspot;
  BigDFT_Energs *energs;
  BigDFT_OptLoop *optloop;
  double *radii;

  GSocket *socket;
  GSource *source;

  g_mem_set_vtable (glib_mem_profiler_table);
  g_type_init ();

  loop = g_main_loop_new (NULL, FALSE);

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
  g_signal_connect(G_OBJECT(wf), "one-wave-ready::1-4-up",
                   G_CALLBACK(onPsiReady), (gpointer)0);

  energs = bigdft_energs_new();
  g_signal_connect(G_OBJECT(energs), "eks-ready",
                   G_CALLBACK(onEKSReady), (gpointer)0);

  denspot = bigdft_localfields_new(BIGDFT_LZD(wf->lzd), in, 0, 1);
  g_signal_connect(G_OBJECT(denspot), "density-ready",
                   G_CALLBACK(onDensityReady), (gpointer)0);

  optloop = bigdft_optloop_new();
  g_signal_connect(G_OBJECT(optloop), "done-wavefunctions",
                   G_CALLBACK(onDoneWavefunctions), (gpointer)0);

  error = (GError*)0;
#ifdef HAVE_GDBUS
  manager = bigdft_dbus_object_manager_client_new_for_bus_sync
    (G_BUS_TYPE_SESSION, G_DBUS_OBJECT_MANAGER_CLIENT_FLAGS_NONE,
     "eu.etsf.BigDFT", "/outputs", NULL, &error);

  g_signal_connect(manager, "object-added",
                   G_CALLBACK(on_object_added), (gpointer)wf);
  g_signal_connect(manager, "object-removed",
                   G_CALLBACK(on_object_removed), (gpointer)loop);
#endif

  if (argc > 1)
    socket = bigdft_signals_client_new(argv[1], NULL, &error);
  else
    socket = bigdft_signals_client_new(g_get_host_name(), NULL, &error);
  if (socket)
    {
      source = bigdft_signals_client_create_source(socket, energs, wf, denspot, optloop,
                                                   NULL, onClosedSocket, loop);
      g_source_attach(source, NULL);

      g_main_loop_run(loop);
    }
  else
    source = (GSource*)0;

  g_object_unref(wf);
  g_object_unref(energs);
  g_object_unref(denspot);
  g_object_unref(optloop);
  bigdft_inputs_free(in);

#ifdef HAVE_GDBUS
  g_object_unref(manager);
#endif
  g_main_loop_unref(loop);
  if (source)
    g_source_unref(source);

  /* g_mem_profile(); */

  return 0;
}
