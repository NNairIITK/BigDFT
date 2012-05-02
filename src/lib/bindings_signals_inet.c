#include <config.h>

#ifdef HAVE_GLIB
#include "bindings_signals.h"

#include <string.h>

#include "bindings.h"


#define PACKET_SIZE 4096

typedef enum
  {
    BIGDFT_SIGNAL_E_READY,
    BIGDFT_SIGNAL_PSI_READY,
    BIGDFT_SIGNAL_DENSPOT_READY,
    BIGDFT_SIGNAL_OPTLOOP_READY,
  } BigDFT_SignalIds;
typedef struct _BigDFT_Signals
{
  BigDFT_SignalIds id;
  guint iter;
  guint kind;
} BigDFT_Signals;

typedef enum
  {
    BIGDFT_SIGNAL_ANSWER_DONE,
    BIGDFT_SIGNAL_ANSWER_WAIT,
    BIGDFT_SIGNAL_ANSWER_GET_E,
    BIGDFT_SIGNAL_ANSWER_GET_OPTLOOP,
    BIGDFT_SIGNAL_ANSWER_GET_DENSPOT,
    BIGDFT_SIGNAL_ANSWER_GET_PSI,
    BIGDFT_SIGNAL_ANSWER_SYNC_OPTLOOP,
  } BigDFT_SignalAnswers;
typedef struct _BigDFT_SignalReply
{
  BigDFT_SignalAnswers id;
  guint ikpt, iorb, kind;
} BigDFT_SignalReply;

static gboolean client_handle_energs(GSocket *socket, BigDFT_Energs *energs,
                                     GCancellable *cancellable, GError **error);


static void onOptLoop(BigDFT_OptLoop *optloop, BigDFT_Energs *energs,
                      GSocket **socket_, BigDFT_OptLoopIds kind)
{
  GSocket *socket = (GSocket*)(*socket_);
  GError *error;
  gssize size, psize;
  BigDFT_Signals signal = {BIGDFT_SIGNAL_OPTLOOP_READY, 0, kind};
  BigDFT_SignalReply answer;
  BigDFT_OptLoop optloop_;
  guint iproc = 0;

  error = (GError*)0;
  if (!socket || !g_socket_is_connected(socket) || g_socket_is_closed(socket))
    return;

  /* We emit the signal on the socket. */
  size = g_socket_send(socket, (const gchar*)(&signal),
                       sizeof(BigDFT_Signals), NULL, &error);
  if (size != sizeof(BigDFT_Signals))
    {
      g_warning("%s", error->message);
      g_error_free(error);
      return;
    }
  
  /* We wait for the answer. */
  do
    {
      size = g_socket_receive(socket, (gchar*)(&answer),
                              sizeof(BigDFT_SignalReply), NULL, &error);
      if (size <= 0)
        {
          /* Connection has been closed on the other side. */
          g_object_unref(socket);
          *socket_ = (gpointer)0;
          return;
        }
      if (size != sizeof(BigDFT_SignalReply))
        {
          g_warning("%s", error->message);
          g_error_free(error);
          return;
        }
      switch (answer.id)
        {
        case BIGDFT_SIGNAL_ANSWER_DONE:
          return;
        case BIGDFT_SIGNAL_ANSWER_WAIT:
          break;
        case BIGDFT_SIGNAL_ANSWER_GET_E:
          size = g_socket_send(socket, (const gchar*)energs,
                               sizeof(BigDFT_Energs), NULL, &error);
          if (size != sizeof(BigDFT_Energs))
            {
              g_warning("%s", error->message);
              g_error_free(error);
            }
          break;
        case BIGDFT_SIGNAL_ANSWER_GET_OPTLOOP:
          size = g_socket_send(socket, (const gchar*)optloop,
                               sizeof(BigDFT_OptLoop), NULL, &error);
          if (size != sizeof(BigDFT_OptLoop))
            {
              g_warning("%s", error->message);
              g_error_free(error);
            }
          break;
        case BIGDFT_SIGNAL_ANSWER_SYNC_OPTLOOP:
          size = 0;
          do
            {
              psize = g_socket_receive(socket, (gchar*)(&optloop_) + size,
                                       sizeof(BigDFT_OptLoop), NULL, &error);
              if (psize == 0)
                {
                  g_object_unref(socket);
                  *socket_ = (gpointer)0;
                  return;
                }
              if (psize <= 0)
                {
                  g_warning("%s", error->message);
                  g_error_free(error);
                  return;
                }
              size += psize;
            }
          while (size < sizeof(BigDFT_OptLoop));
          memcpy((void*)((char*)optloop + sizeof(GObject)),
                 (const void*)((char*)(&optloop_) + sizeof(GObject)),
                 sizeof(BigDFT_OptLoop) - sizeof(GObject) - sizeof(gpointer));
          bigdft_optloop_sync_to_fortran(optloop);
          FC_FUNC_(optloop_bcast, OPTLOOP_BCAST)(optloop->data, &iproc);
          break;
        default:
          g_warning("Server: wrong client answer after optloop emission.");
          return;
        }
    }
  while (1);
}
static gboolean client_handle_optloop(GSocket *socket, BigDFT_OptLoop *optloop, BigDFT_OptLoopIds kind,
                                      BigDFT_Energs *energs, GCancellable *cancellable, GError **error)
{
  BigDFT_SignalReply answer;
  gssize size, psize;
  BigDFT_OptLoop optloop_;
  gboolean ret;

  ret = client_handle_energs(socket, energs, cancellable, error);
  if (!ret)
    return FALSE;

  /* g_print("Client: send get E.\n"); */
  answer.id = BIGDFT_SIGNAL_ANSWER_GET_OPTLOOP;
  size = g_socket_send(socket, (const gchar*)(&answer),
                       sizeof(BigDFT_SignalReply), cancellable, error);
  /* g_print("Client: send %ld / %ld.\n", size, sizeof(BigDFT_SignalReply)); */
  if (size != sizeof(BigDFT_SignalReply))
    return FALSE;
  psize = 0;
  do
    {
      size = g_socket_receive(socket, (gchar*)(&optloop_) + psize,
                              sizeof(BigDFT_OptLoop), cancellable, error);
      if (size == 0)
        *error = g_error_new(G_IO_ERROR, G_IO_ERROR_CLOSED,
                             "Connection closed by peer.");
      if (size <= 0)
        return FALSE;
      psize += size;
    }
  while (psize < sizeof(BigDFT_OptLoop));
  /* g_print("Client: receive %ld / %ld.\n", psize, sizeof(BigDFT_OptLoop)); */
  memcpy((void*)((char*)optloop + sizeof(GObject)),
         (const void*)((char*)(&optloop_) + sizeof(GObject)),
         sizeof(BigDFT_OptLoop) - sizeof(GObject) - sizeof(gpointer));
  bigdft_optloop_emit(optloop, kind, energs);
  
  return TRUE;
}
static void onOptLoopSyncInet(BigDFT_OptLoop *optloop, gpointer data)
{
  GSocket *socket = (GSocket*)data;
  BigDFT_SignalReply answer;
  gssize size;
  GError *error;

  error = (GError*)0;

  answer.id = BIGDFT_SIGNAL_ANSWER_SYNC_OPTLOOP;
  size = g_socket_send(socket, (const gchar*)(&answer),
                       sizeof(BigDFT_SignalReply), NULL, &error);
  /* g_print("Client: send %ld / %ld.\n", size, sizeof(BigDFT_SignalReply)); */
  if (size != sizeof(BigDFT_SignalReply))
    {
      g_warning("%s", error->message);
      g_error_free(error);
      return;
    }

  size = g_socket_send(socket, (gchar*)optloop,
                       sizeof(BigDFT_OptLoop), NULL, &error);
  if (size != sizeof(BigDFT_OptLoop))
    {
      g_warning("%s", error->message);
      g_error_free(error);
      return;
    }
}
void onIterHamInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data)
{
  onOptLoop(optloop, energs, (GSocket**)data, BIGDFT_OPTLOOP_ITER_HAMILTONIAN);
}
void onIterSubInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data)
{
  onOptLoop(optloop, energs, (GSocket**)data, BIGDFT_OPTLOOP_ITER_SUBSPACE);
}
void onIterWfnInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data)
{
  onOptLoop(optloop, energs, (GSocket**)data, BIGDFT_OPTLOOP_ITER_WAVEFUNCTIONS);
}
void onDoneHamInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data)
{
  onOptLoop(optloop, energs, (GSocket**)data, BIGDFT_OPTLOOP_DONE_HAMILTONIAN);
}
void onDoneSubInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data)
{
  onOptLoop(optloop, energs, (GSocket**)data, BIGDFT_OPTLOOP_DONE_SUBSPACE);
}
void onDoneWfnInet(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer *data)
{
  onOptLoop(optloop, energs, (GSocket**)data, BIGDFT_OPTLOOP_DONE_WAVEFUNCTIONS);
}

static void _onDensPotReady(BigDFT_LocalFields *localfields, guint iter,
                            GSocket **socket_, BigDFT_DensPotIds kind)
{
  GSocket *socket = *socket_;
  GError *error;
  gssize size;
  BigDFT_Signals signal = {BIGDFT_SIGNAL_DENSPOT_READY, iter, kind};
  BigDFT_SignalReply answer;
  guint iproc = 0, i, s, nvects, sizeData[2], new;
  f90_pointer_double tmp;

  error = (GError*)0;
  if (!socket || !g_socket_is_connected(socket) || g_socket_is_closed(socket))
    return;

  /* We emit the signal on the socket. */
  size = g_socket_send(socket, (const gchar*)(&signal),
                       sizeof(BigDFT_Signals), NULL, &error);
  if (size != sizeof(BigDFT_Signals))
    {
      g_warning("%s", error->message);
      g_error_free(error);
      return;
    }
  
  /* We wait for the answer. */
  do
    {
      size = g_socket_receive(socket, (gchar*)(&answer),
                              sizeof(BigDFT_SignalReply), NULL, &error);
      if (size <= 0)
        {
          /* Connection has been closed on the other side. */
          g_object_unref(socket);
          *socket_ = (gpointer)0;
          return;
        }
      if (size != sizeof(BigDFT_SignalReply))
        {
          g_warning("(%ld / %ld) %s", size, sizeof(BigDFT_SignalReply), error->message);
          g_error_free(error);
          return;
        }
      if (answer.id == BIGDFT_SIGNAL_ANSWER_DONE)
        return;
  
      s = localfields->ni[0] * localfields->ni[1] * localfields->ni[2];
      /* We send the size of the data to send. */
      sizeData[0] = s;
      sizeData[1] = PACKET_SIZE;
      size = g_socket_send(socket, (const gchar*)sizeData,
                           sizeof(guint) * 2, NULL, &error);
      if (error)
        {
          g_warning("%s", error->message);
          g_error_free(error);
          return;
        }
  
      /* We do a full_local_potential to gather all the data. */
      F90_1D_POINTER_INIT(&tmp);
      switch (kind)
        {
        case BIGDFT_DENSPOT_DENSITY:
          FC_FUNC_(denspot_full_density, DENSPOT_FULL_DENSITY)(localfields->data,
                                                                   &tmp, &iproc, &new);
          break;
        case BIGDFT_DENSPOT_V_EXT:
          FC_FUNC_(denspot_full_v_ext, DENSPOT_FULL_V_EXT)(localfields->data,
                                                               &tmp, &iproc, &new);
          break;
        }

      /* We send the whole values. */
      nvects = sizeof(double) * s / PACKET_SIZE + 1;
      size = 0;
      for (i = 0; i < nvects; i++)
        {
          size += g_socket_send(socket, (const gchar*)(tmp.data) + i * PACKET_SIZE,
                                (i < nvects - 1)?PACKET_SIZE:sizeof(double) * s -
                                i * PACKET_SIZE, NULL, &error);
          if (error)
            {
              g_warning("%s", error->message);
              g_error_free(error);
              break;
            }
        }

      if (new > 0)
        FC_FUNC_(deallocate_double_1d, DEALLOCATE_DOUBLE_1D)(&tmp);
    }
  while (1);
}

void onDensityReadyInet(BigDFT_LocalFields *localfields, guint iter, gpointer *data)
{
  _onDensPotReady(localfields, iter, (GSocket**)data, BIGDFT_DENSPOT_DENSITY);
}
void onVExtReadyInet(BigDFT_LocalFields *localfields, gpointer *data)
{
  _onDensPotReady(localfields, 0, (GSocket**)data, BIGDFT_DENSPOT_V_EXT);
}

static gboolean client_handle_denspot(GSocket *socket, BigDFT_LocalFields *denspot,
                                      guint iter, BigDFT_DensPotIds kind,
                                      GCancellable *cancellable, GError **error)
{
  BigDFT_SignalReply answer;
  gssize size, psize;
  guint sizeData[2];
  gchar *target;

  /* g_print("Client: send get denspot.\n"); */
  answer.id = BIGDFT_SIGNAL_ANSWER_GET_DENSPOT;
  answer.kind = kind;
  size = g_socket_send(socket, (const gchar*)(&answer),
                       sizeof(BigDFT_SignalReply), cancellable, error);
  /* g_print("Client: send %ld / %ld.\n", size, sizeof(BigDFT_SignalReply)); */
  if (size != sizeof(BigDFT_SignalReply))
    return FALSE;

  size = g_socket_receive(socket, (gchar*)sizeData,
                          sizeof(guint) * 2, cancellable, error);
  if (size == 0)
    *error = g_error_new(G_IO_ERROR, G_IO_ERROR_CLOSED,
                         "Connection closed by peer.");
  if (size <= 0)
    return FALSE;
          
  psize = 0;
  if (kind == BIGDFT_DENSPOT_DENSITY)
    target = (gchar*)denspot->rhov;
  else if (kind == BIGDFT_DENSPOT_V_EXT)
    target = (gchar*)denspot->v_ext;
  do
    {
      size = g_socket_receive(socket, target + psize,
                              sizeData[1], cancellable, error);
      if (size == 0)
        *error = g_error_new(G_IO_ERROR, G_IO_ERROR_CLOSED,
                             "Connection closed by peer.");
      if (size <= 0)
        return FALSE;
      psize += size;
    }
  while (psize < sizeof(double) * sizeData[0]);
  /* g_print("Client: receive %ld / %ld.\n", psize, sizeof(double) * sizeData[0]); */
  if (kind == BIGDFT_DENSPOT_DENSITY)
    {
      denspot->rhov_is = BIGDFT_RHO_IS_ELECTRONIC_DENSITY;
      bigdft_localfields_emit_rhov(denspot, iter);
    }
  else if (kind == BIGDFT_DENSPOT_V_EXT)
      bigdft_localfields_emit_v_ext(denspot);

  return TRUE;
}

void onEKSReadyInet(BigDFT_Energs *energs, guint iter, gpointer *data)
{
  GSocket *socket = (GSocket*)(*data);
  GError *error;
  gssize size;
  BigDFT_Signals signal = {BIGDFT_SIGNAL_E_READY, iter, BIGDFT_ENERGS_EKS};
  BigDFT_SignalReply answer;

  error = (GError*)0;
  if (!socket || !g_socket_is_connected(socket) || g_socket_is_closed(socket))
    return;

  /* We emit the signal on the socket. */
  size = g_socket_send(socket, (const gchar*)(&signal),
                       sizeof(BigDFT_Signals), NULL, &error);
  if (size != sizeof(BigDFT_Signals))
    {
      g_warning("%s", error->message);
      g_error_free(error);
      return;
    }
  
  /* We wait for the answer. */
  do
    {
      size = g_socket_receive(socket, (gchar*)(&answer),
                              sizeof(BigDFT_SignalReply), NULL, &error);
      if (size <= 0)
        {
          /* Connection has been closed on the other side. */
          g_object_unref(socket);
          *data = (gpointer)0;
          return;
        }
      if (size != sizeof(BigDFT_SignalReply))
        {
          g_warning("%s", error->message);
          g_error_free(error);
          return;
        }
      if (answer.id == BIGDFT_SIGNAL_ANSWER_DONE)
        return;
  
      /* We send the energy values. */
      size = g_socket_send(socket, (const gchar*)energs,
                           sizeof(BigDFT_Energs), NULL, &error);
      if (size != sizeof(BigDFT_Energs))
        {
          g_warning("%s", error->message);
          g_error_free(error);
        }
    }
  while (1);
}
static gboolean client_handle_energs(GSocket *socket, BigDFT_Energs *energs,
                                     GCancellable *cancellable, GError **error)
{
  BigDFT_SignalReply answer;
  gssize size, psize;
  BigDFT_Energs energs_;

  /* g_print("Client: send get E.\n"); */
  answer.id = BIGDFT_SIGNAL_ANSWER_GET_E;
  size = g_socket_send(socket, (const gchar*)(&answer),
                       sizeof(BigDFT_SignalReply), cancellable, error);
  /* g_print("Client: send %ld / %ld.\n", size, sizeof(BigDFT_SignalReply)); */
  if (size != sizeof(BigDFT_SignalReply))
    return FALSE;
  psize = 0;
  do
    {
      size = g_socket_receive(socket, (gchar*)(&energs_) + psize,
                              sizeof(BigDFT_Energs), cancellable, error);
      if (size == 0)
        *error = g_error_new(G_IO_ERROR, G_IO_ERROR_CLOSED,
                             "Connection closed by peer.");
      if (size <= 0)
        return FALSE;
      psize += size;
    }
  while (psize < sizeof(BigDFT_Energs));
  /* g_print("Client: receive %ld / %ld.\n", psize, sizeof(BigDFT_Energs)); */
  memcpy((void*)((char*)energs + sizeof(GObject)),
         (const void*)((char*)(&energs_) + sizeof(GObject)),
         sizeof(BigDFT_Energs) - sizeof(GObject) - sizeof(gpointer));
  /* g_print("Client: emitting signal kind %d.\n", signal.kind); */
  
  return TRUE;
}

void onPsiReadyInet(BigDFT_Wf *wf, guint iter, gpointer *data)
{
  GSocket *socket = (GSocket*)(*data);
  GError *error;
  gssize size;
  BigDFT_Signals signal = {BIGDFT_SIGNAL_PSI_READY, iter, 0};
  BigDFT_SignalReply answer;
  double *psic;
  guint i, psiSize, nvects, sizeData[2];
  gboolean copy;

  error = (GError*)0;
  if (!socket || !g_socket_is_connected(socket) || g_socket_is_closed(socket))
    return;

  /* Get the wave-function size. */
  bigdft_wf_get_psi_compress(wf, 1, 1, BIGDFT_SPIN_UP, BIGDFT_PARTIAL_DENSITY,
                             &signal.kind, 0);

  /* We emit the signal on the socket. */
  size = g_socket_send(socket, (const gchar*)(&signal),
                       sizeof(BigDFT_Signals), NULL, &error);
  if (size != sizeof(BigDFT_Signals))
    {
      g_warning("%s", error->message);
      g_error_free(error);
      return;
    }
  
  /* We wait for the answer. */
  do
    {
      size = g_socket_receive(socket, (gchar*)(&answer),
                              sizeof(BigDFT_SignalReply), NULL, &error);
      if (size <= 0)
        {
          /* Connection has been closed on the other side. */
          g_object_unref(socket);
          *data = (gpointer)0;
          return;
        }
      if (size != sizeof(BigDFT_SignalReply))
        {
          g_warning("%s", error->message);
          g_error_free(error);
          return;
        }
      if (answer.id == BIGDFT_SIGNAL_ANSWER_DONE)
        return;

      /* We gather the data to send. */
      copy = FALSE;
      psic = (double*)bigdft_wf_get_psi_compress(wf, answer.ikpt, answer.iorb,
                                                 answer.kind, BIGDFT_PARTIAL_DENSITY,
                                                 &psiSize, 0);
      if (psiSize == 0)
        g_warning("Required psi is not available");
      if (!psic)
        {
          /* Data are on processus iproc, so we copy them. */
          psic = g_malloc(sizeof(double) * psiSize);
          if (!bigdft_wf_copy_psi_compress(wf, answer.ikpt, answer.iorb,
                                           answer.kind, BIGDFT_PARTIAL_DENSITY, 0,
                                           psic, psiSize))
            {
              g_free(psic);
              g_warning("Required psi cannot be copied");
              psiSize = 0;
            }
          copy = TRUE;
        }
      /* We send the size of the data to send. */
      sizeData[0] = psiSize;
      sizeData[1] = PACKET_SIZE;
      size = g_socket_send(socket, (const gchar*)sizeData,
                           sizeof(guint) * 2, NULL, &error);
      if (error)
        {
          g_warning("%s", error->message);
          g_error_free(error);
          return;
        }
  
      if (psiSize > 0)
        {
          /* We send the compressed values. */
          nvects = sizeof(double) * psiSize / PACKET_SIZE + 1;
          size = 0;
          for (i = 0; i < nvects; i++)
            {
              size += g_socket_send(socket, (const gchar*)psic + i * PACKET_SIZE,
                                    (i < nvects - 1)?PACKET_SIZE:sizeof(double) * psiSize -
                                    i * PACKET_SIZE, NULL, &error);
              if (error)
                {
                  g_warning("%s", error->message);
                  g_error_free(error);
                  break;
                }
            }
          if (copy)
            g_free(psic);
        }
    }
  while (1);
}

static gboolean client_handle_wf(GSocket *socket, BigDFT_Wf *wf, guint iter, guint psiSize,
                                 guint ikpt, guint iorb, BigDFT_Spin ispin, GQuark quark,
                                 GCancellable *cancellable, GError **error)
{
  BigDFT_SignalReply answer;
  gssize size, psize;
  GArray *psic;
  guint sizeData[2];

  psic = g_array_sized_new(FALSE, FALSE, sizeof(double), psiSize);

  /* g_print("Client: send get psi.\n"); */
  answer.id = BIGDFT_SIGNAL_ANSWER_GET_PSI;
  answer.ikpt = ikpt;
  answer.iorb = iorb;
  answer.kind = ispin;
  size = g_socket_send(socket, (const gchar*)(&answer),
                       sizeof(BigDFT_SignalReply), cancellable, error);
  /* g_print("Client: send %ld / %ld.\n", size, sizeof(BigDFT_SignalReply)); */
  if (size != sizeof(BigDFT_SignalReply))
    return FALSE;

  size = g_socket_receive(socket, (gchar*)sizeData,
                          sizeof(guint) * 2, cancellable, error);
  if (size == 0)
    *error = g_error_new(G_IO_ERROR, G_IO_ERROR_CLOSED,
                         "Connection closed by peer.");
  if (size <= 0)
    return FALSE;
  if (sizeData[0] == 0 || psiSize != sizeData[0])
    {
      *error = g_error_new(G_IO_ERROR, G_IO_ERROR_INVALID_DATA,
                           "Unable to retrieve psi.");
      g_array_free(psic, TRUE);
      return FALSE;
    }
          
  psize = 0;
  do
    {
      size = g_socket_receive(socket, psic->data + psize,
                              sizeData[1], cancellable, error);
      if (size == 0)
        *error = g_error_new(G_IO_ERROR, G_IO_ERROR_CLOSED,
                             "Connection closed by peer.");
      if (size <= 0)
        {
          g_array_free(psic, TRUE);
          return FALSE;
        }
      psize += size;
    }
  while (psize < sizeof(double) * psiSize);
  /* g_print("Client: receive %ld / %ld.\n", psize, sizeof(double) * sizeData[0]); */
  /* g_print("Client: emitting signal.\n"); */
  psic = g_array_set_size(psic, psiSize);
  bigdft_wf_emit_one_wave(wf, iter, psic, quark,
                          answer.ikpt, answer.iorb, answer.kind);
  g_array_unref(psic);

  return TRUE;
}

gboolean onClientConnection(GSocket *socket, GIOCondition condition,
                            gpointer user_data)
{
  GError *error;
  BigDFT_Main *bmain = (BigDFT_Main*)user_data;

  error = (GError*)0;

  if ((condition & G_IO_IN) > 0)
    {
      /* g_print("Server: client requesting connection.\n"); */
      bmain->recv = g_socket_accept(socket, NULL, &error);
      if (!bmain->recv)
        {
          g_warning("Server: %s", error->message);
          g_error_free(error);
        }
      g_socket_set_blocking(bmain->recv, TRUE);
      /* g_print("Server: client connected.\n"); */
    }

  if ((condition & G_IO_HUP) > 0 && bmain->recv)
    {
      g_object_unref(bmain->recv);
      bmain->recv = (GSocket*)0;
    }

  return TRUE;
}

GSocket* bigdft_signals_client_new(const gchar *hostname,
                                   GCancellable *cancellable, GError **error)
{
  GSocket *socket;
  GResolver *dns;
  GList *lst, *tmp;
  GSocketAddress *sockaddr;
  gboolean connect;

  /* g_print("Create a socket for hostname '%s'.\n", hostname); */
  socket = g_socket_new(G_SOCKET_FAMILY_IPV4, G_SOCKET_TYPE_STREAM,
                        G_SOCKET_PROTOCOL_DEFAULT, error);
  if (!socket)
    return socket;

  g_socket_set_blocking(socket, TRUE);

  connect = FALSE;
  dns = g_resolver_get_default();
  lst = g_resolver_lookup_by_name(dns, hostname, cancellable, error);
  g_object_unref(dns);
  for (tmp = lst; tmp && !connect; tmp = g_list_next(tmp))
    {
      if (*error)
        {
          g_warning("%s", (*error)->message);
          g_error_free(*error);
          *error = (GError*)0;
        }
      sockaddr = g_inet_socket_address_new((GInetAddress*)tmp->data, (guint16)91691);
      connect = g_socket_connect(socket, sockaddr, cancellable, error);
      /* g_print(" | try to connect to '%s' -> %d.\n", */
              /* g_inet_address_to_string((GInetAddress*)tmp->data), connect); */
      g_object_unref(sockaddr);
    }
  g_resolver_free_addresses(lst);
  if (!connect)
    {
      g_object_unref(socket);
      socket = (GSocket*)0;
    }
  return socket;
}


static gboolean bigdft_signals_client_handle(GSocket *socket, BigDFT_Energs *energs,
                                             BigDFT_Wf *wf, BigDFT_LocalFields *denspot,
                                             BigDFT_OptLoop *optloop,
                                             GCancellable *cancellable, GError **error)
{
  BigDFT_Signals signal;
  BigDFT_SignalReply answer;
  gssize size;
  gboolean ret;
  guint ikpt, iorb;
  GQuark quark;
  gchar *details;

  size = g_socket_receive(socket, (gchar*)(&signal),
                          sizeof(BigDFT_Signals), cancellable, error);
  if (size == 0)
    {
      *error = g_error_new(G_IO_ERROR, G_IO_ERROR_CLOSED,
                           "Connection closed by peer.");
      return FALSE;
    }
  if (size != sizeof(BigDFT_Signals))
    return FALSE;

  /* g_print("Client: get signal %d at iter %d.\n", signal.id, signal.iter); */
  ret = TRUE;
  switch (signal.id)
    {
    case BIGDFT_SIGNAL_E_READY:
      switch (signal.kind)
        {
        case BIGDFT_ENERGS_EKS:
          if (g_signal_has_handler_pending(energs,
                                           g_signal_lookup("eks-ready",
                                                           BIGDFT_ENERGS_TYPE),
                                           0, FALSE))
            ret = client_handle_energs(socket, energs, cancellable, error);
          if (ret)
            bigdft_energs_emit(energs, signal.iter, signal.kind);
          break;
        }
      break;
    case BIGDFT_SIGNAL_PSI_READY:
      /* We test event pending for all possible wavefubctions. */
      for (ikpt = 1; ikpt <= BIGDFT_ORBS(wf)->nkpts; ikpt++)
        {
          for (iorb = 1; iorb <= BIGDFT_ORBS(wf)->norbu; iorb++)
            {
              details = g_strdup_printf("%d-%d-up", ikpt, iorb);
              quark = g_quark_from_string(details);
              g_free(details);
              if (g_signal_has_handler_pending(wf,
                                               g_signal_lookup("one-wave-ready",
                                                               BIGDFT_WF_TYPE),
                                               quark, FALSE))
                ret = client_handle_wf(socket, wf, signal.iter, signal.kind,
                                       ikpt, iorb, BIGDFT_SPIN_UP, quark,
                                       cancellable, error);
              if (!ret)
                break;
            }
          if (!ret)
            break;
        }
      break;
    case BIGDFT_SIGNAL_DENSPOT_READY:
      switch (signal.kind)
        {
        case BIGDFT_DENSPOT_DENSITY:
          if (g_signal_has_handler_pending(denspot,
                                           g_signal_lookup("density-ready",
                                                           BIGDFT_LOCALFIELDS_TYPE),
                                           0, FALSE))
            ret = client_handle_denspot(socket, denspot, signal.iter, signal.kind,
                                        cancellable, error);
          break;
        case  BIGDFT_DENSPOT_V_EXT:
          if (g_signal_has_handler_pending(denspot,
                                           g_signal_lookup("v-ext-ready",
                                                           BIGDFT_LOCALFIELDS_TYPE),
                                           0, FALSE))
            ret = client_handle_denspot(socket, denspot, signal.iter, signal.kind,
                                        cancellable, error);
          break;
        }
      break;
    case BIGDFT_SIGNAL_OPTLOOP_READY:
      switch (signal.kind)
        {
        case BIGDFT_OPTLOOP_ITER_HAMILTONIAN:
          if (g_signal_has_handler_pending(optloop,
                                           g_signal_lookup("iter-hamiltonian",
                                                           BIGDFT_OPTLOOP_TYPE),
                                           0, FALSE))
            ret = client_handle_optloop(socket, optloop, signal.kind, energs, cancellable, error);
          break;
        case BIGDFT_OPTLOOP_ITER_SUBSPACE:
          if (g_signal_has_handler_pending(optloop,
                                           g_signal_lookup("iter-subspace",
                                                           BIGDFT_OPTLOOP_TYPE),
                                           0, FALSE))
            ret = client_handle_optloop(socket, optloop, signal.kind, energs, cancellable, error);
          break;
        case BIGDFT_OPTLOOP_ITER_WAVEFUNCTIONS:
          if (g_signal_has_handler_pending(optloop,
                                           g_signal_lookup("iter-wavefunctions",
                                                           BIGDFT_OPTLOOP_TYPE),
                                           0, FALSE))
            ret = client_handle_optloop(socket, optloop, signal.kind, energs, cancellable, error);
          break;
        case BIGDFT_OPTLOOP_DONE_HAMILTONIAN:
          if (g_signal_has_handler_pending(optloop,
                                           g_signal_lookup("done-hamiltonian",
                                                           BIGDFT_OPTLOOP_TYPE),
                                           0, FALSE))
            ret = client_handle_optloop(socket, optloop, signal.kind, energs, cancellable, error);
          break;
        case BIGDFT_OPTLOOP_DONE_SUBSPACE:
          if (g_signal_has_handler_pending(optloop,
                                           g_signal_lookup("done-subspace",
                                                           BIGDFT_OPTLOOP_TYPE),
                                           0, FALSE))
            ret = client_handle_optloop(socket, optloop, signal.kind, energs, cancellable, error);
          break;
        case BIGDFT_OPTLOOP_DONE_WAVEFUNCTIONS:
          if (g_signal_has_handler_pending(optloop,
                                           g_signal_lookup("done-wavefunctions",
                                                           BIGDFT_OPTLOOP_TYPE),
                                           0, FALSE))
            ret = client_handle_optloop(socket, optloop, signal.kind, energs, cancellable, error);
          break;
        }
      break;
    default:
      break;
    }
  if (ret)
    {
      /* g_print("Client: send signal done %d.\n", signal.id); */
      answer.id = BIGDFT_SIGNAL_ANSWER_DONE;
      size = g_socket_send(socket, (const gchar*)(&answer),
                           sizeof(BigDFT_SignalReply), cancellable, error);
      if (size != sizeof(BigDFT_SignalReply))
        return FALSE;
    }
  return ret;
}

static gboolean onClientTransfer(GSocket *socket, GIOCondition condition,
                                 gpointer user_data)
{
  BigDFT_Main *bmain = (BigDFT_Main*)user_data;
  GError *error;
  
  error = (GError*)0;

  if (bmain->cancellable && g_cancellable_is_cancelled(bmain->cancellable))
    {
      g_source_destroy(g_main_current_source());
      if (bmain->destroy)
        bmain->destroy(bmain->destroyData);
      return FALSE;
    }

  if ((condition & G_IO_IN) > 0)
    {
      g_signal_handler_unblock(G_OBJECT(bmain->optloop), bmain->optloop_sync);
      if (!bigdft_signals_client_handle(socket, bmain->energs, bmain->wf, bmain->denspot,
                                        bmain->optloop, NULL, &error))
        {
          g_signal_handler_block(G_OBJECT(bmain->optloop), bmain->optloop_sync);
          if (error->code != G_IO_ERROR_CLOSED)
            g_warning("Client: %s", error->message);
          else
            {
              g_source_destroy(g_main_current_source());
              if (bmain->destroy)
                bmain->destroy(bmain->destroyData);
            }
          return (error->code != G_IO_ERROR_CLOSED);
        }
      g_signal_handler_block(G_OBJECT(bmain->optloop), bmain->optloop_sync);
    }

  if ((condition & G_IO_HUP) > 0)
    return FALSE;

  return TRUE;
}

GSource* bigdft_signals_client_create_source(GSocket *socket, BigDFT_Energs *energs,
                                             BigDFT_Wf *wf, BigDFT_LocalFields *denspot,
                                             BigDFT_OptLoop *optloop,
                                             GCancellable *cancellable,
                                             GDestroyNotify destroy, gpointer data)
{
  GSource *source;
  BigDFT_Main *bmain;

  bmain = g_malloc0(sizeof(BigDFT_Main));
  bmain->wf = wf;
  if (wf)
    g_object_ref(wf);
  bmain->energs = energs;
  if (energs)
    g_object_ref(energs);
  bmain->denspot = denspot;
  if (denspot)
    g_object_ref(denspot);
  bmain->optloop = optloop;
  if (optloop)
    {
      g_object_ref(optloop);
      bmain->optloop_sync = g_signal_connect(G_OBJECT(bmain->optloop), "sync-fortran",
                                             G_CALLBACK(onOptLoopSyncInet), (gpointer)socket);
      g_signal_handler_block(G_OBJECT(bmain->optloop), bmain->optloop_sync);
    }
  bmain->cancellable = cancellable;
  bmain->destroy = destroy;
  bmain->destroyData = data;
  
  source = g_socket_create_source(socket, G_IO_IN | G_IO_HUP, cancellable);
  g_source_set_callback(source, (GSourceFunc)onClientTransfer,
                        bmain, bigdft_signals_free_main);
  
  return source;
}
#endif
