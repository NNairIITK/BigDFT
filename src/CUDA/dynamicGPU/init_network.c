#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/types.h>    
#include <sys/socket.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <stdio.h>
#include <signal.h>
#include <netdb.h>
#include <errno.h>




#define SIZE_LINE  40
#define SIZE_NAME 100
#define SIZE_DATA 1024
int init_network(const int NUM_PARTICIPANTS)
{
  char *nomfic = "yo.yo";

  int participants[NUM_PARTICIPANTS];
  FILE * es;
  int lu = 0;

  int n,lectOK,currNum,voisin;


  int     sock;

  char buffer[SIZE_DATA];
  socklen_t tosize;

  char curr_sock_path[SIZE_NAME];
  struct  sockaddr_un servaddr; /* server adress (current address) */

  char to_sock_path[SIZE_NAME];
  struct  sockaddr_un toaddr; /* neigbour adress */

  bzero(&servaddr, sizeof(servaddr));
  bzero(&toaddr, sizeof(servaddr));
 
  


  if ((sock = socket(AF_UNIX, SOCK_DGRAM, 0)) < 0) 
    {
      perror("server: socket");
      exit(1);
    }

 
  snprintf(curr_sock_path,SIZE_NAME,"/tmp/sock%i",getpid());
  unlink(curr_sock_path);

  servaddr.sun_family = AF_UNIX;
  strncpy(servaddr.sun_path, curr_sock_path,SIZE_NAME);
  
  if (bind(sock, (struct sockaddr*)&servaddr, sizeof(servaddr)) < 0) 
    {
      close(sock);
      perror("server: bind");
      return 2;
    } 






  es = fopen(nomfic, "a+");
  fprintf(es, "%i\n", getpid());
  //  printf( "%i\n", getpid());
  while(lu < NUM_PARTICIPANTS)
    {
      lu = 0;
      
      sleep(1);
      rewind(es);
      
      do
	{
	  lectOK = fscanf(es, "%i", &n);
	  if (lectOK == 1) 
	    {
	      participants[lu] = n;
	      if(n == getpid())
		{
		  currNum = lu;
		}
	      ++lu;
	      // printf("%i-- %i\n", getpid(), n);
	    }
	}
      while (lectOK == 1 && fgetc(es) != EOF);
    }
   fclose(es);
   remove(nomfic);


   printf("pid chef %d-- %d, currnum %d \n", getpid(), participants[0],currNum);

   voisin = (currNum + 1)%NUM_PARTICIPANTS;

   printf(" %d-- voisin : %d, pid %d\n",getpid(),voisin,participants[voisin]);


   //ok we have the network topology !
  
   snprintf(to_sock_path,SIZE_NAME,"/tmp/sock%i",participants[voisin]);
    
   toaddr.sun_family = AF_UNIX;
   strncpy(toaddr.sun_path, to_sock_path,SIZE_NAME);
   
   snprintf(buffer,SIZE_DATA,"de %i vers %i\n",currNum,voisin);


   if(sendto(sock, buffer, strlen(buffer), 0, (struct sockaddr*)&toaddr, sizeof(toaddr)) < 0)
     {
       perror("sendto()");
       exit(errno);
     }

   bzero(buffer, sizeof(char)*SIZE_DATA);
   if((n = recvfrom(sock, buffer, SIZE_DATA - 1, 0, NULL, &tosize)) < 0)
     {
       perror("recvfrom()");
       exit(errno);
     }


   printf("%i recoit : %s\n",getpid(),buffer);
   close(sock);
     unlink(curr_sock_path);
  return 0;
}


int send_next()
{
  if(sendto(sock, buffer, strlen(buffer), 0, (struct sockaddr*)&toaddr, sizeof(toaddr)) < 0)
    {
      perror("sendto()");
      return errno;
    }
}

int recv_prev()
{
  if((n = recvfrom(sock, buffer, SIZE_DATA - 1, 0, NULL, &tosize)) < 0)
    {
      perror("recvfrom()");
      return errno;
    }
}
