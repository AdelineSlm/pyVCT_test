// file: write.c
#include <fcntl.h>
#include <sys/mman.h>
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
inline void write_increment(int increment)
{
	FILE *ofp;

	ofp = fopen("./output/increment_number.out","w");
	fprintf(ofp,"%d",increment);
	//fflush(ofp); 
  fclose(ofp);
}


inline void mmap_write_increment(int increment)
{
	//FILE *ofp;

  int fd = open("./output/increment_number.out", O_RDWR, 0777);

  char str[32] = {0};
  sprintf(str, "%d", increment);

  void* ptr = mmap(NULL, sizeof(str), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_SYNC, fd, 0);
  memcpy(ptr, str, strlen(str));

  munmap(ptr);
  close(fd);


	//ofp = fopen("./output/increment_number.out","w");
	//fprintf(ofp,"%d",increment);
	//fflush(ofp); 
  //fclose(ofp);
}


////////////////////////////////////////////////////////////////////////////////
inline void write_cells(VOX* pv, int increment, int NVX, int NVY)
{
	int v;
  int vx,vy;
  char filename[40];
  char astring[20];
  FILE *ofp;

  myitostr(increment, astring);
	strcpy(filename, "./output/ctags");
  strcat(filename, astring);
  strcat(filename, ".out");

	ofp = fopen(filename,"w");

  for(vx=0; vx<NVX; vx++) {
      for (vy=0; vy<NVY; vy++) {
        v = vx + vy * NVX;
   		   fprintf(ofp ,"%d ", pv[v].ctag);
      }
      fprintf(ofp, "\n");
  }
  //fflush(ofp);  
  fclose(ofp);
}

////////////////////////////////////////////////////////////////////////////////
inline void write_types(int* types, int NRc)
{
  int v;
  FILE *ofp;

  ofp = fopen("./output/types.out","w");

  for(v=0; v<NRc; v++) {
      fprintf(ofp ,"%d ", types[1+v]);
  }
  //fflush(ofp);  
  fclose(ofp);
}

////////////////////////////////////////////////////////////////////////////////
inline void write_contacts(VOX* pv, int increment, int NVX, int NVY)
{
  int v;
    int vx,vy;
    char filename[40];
    char astring[20];
    FILE *ofp;

    myitostr(increment, astring);
    strcpy(filename, "./output/contactM");
    strcat(filename, astring);
    strcat(filename, ".out");

    ofp = fopen(filename,"w");
    for(vx=0; vx<NVX; vx++) {
        for (vy=0; vy<NVY; vy++) {
            v = vx + vy * NVX;
        fprintf(ofp ,"%d ", pv[v].contact);
        }
        fprintf(ofp, "\n");
    }
    //fflush(ofp);  
    fclose(ofp);
}

////////////////////////////////////////////////////////////////////////////////
inline void write_fibers(FIBERS* pf, int NVX, int NVY)
{
  int v;
  int vx,vy;
  FILE *ofp;

  ofp = fopen("./output/fib.out","w");

  for(vx=0; vx<NVX; vx++) {
      for (vy=0; vy<NVY; vy++) {
        v = vx + vy * NVX;
        fprintf(ofp ,"%d ", pf[v].Q);
      }
      fprintf(ofp, "\n");
  }
  //fflush(ofp);  
  fclose(ofp);
}

////////////////////////////////////////////////////////////////////////////////
inline void read_cells(VOX* pv, int* types, int NRc, char filename_ctag[40], char filename_cont[40], char filename_types[40])
{
  int v,k,vx,vy;
  FILE *ofp;

  ofp = fopen(filename_ctag,"r");
    for(v=0; v<NV; v++){
        vx = v/NVY; vy = v%NVY;
        k = fscanf(ofp ,"%d", &pv[vx + vy * NVX].ctag);
    }
    //fflush(ofp);  
    fclose(ofp);

    ofp = fopen(filename_cont,"r");
    for(v=0; v<NV; v++){
        vx = v/NVY; vy = v%NVY;
        k = fscanf(ofp ,"%d", &pv[vx + vy * NVX].contact);
    }
    //fflush(ofp);  
    fclose(ofp);

    ofp = fopen(filename_types,"r");
    for(v=0; v<NRc; v++){
        k = fscanf(ofp ,"%d", &types[v+1]);
    }
    //fflush(ofp);  
    fclose(ofp);
}
