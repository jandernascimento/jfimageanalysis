#include <stdlib.h>
#include <stdio.h>
#include "Util.h"


// HEADER
void printimage(gray* image,int cols, int rows, int maxval);
gray* readimage(char* filepath);
int rows=0,cols=0,maxval=0;
#define DPC 3 //data per cell information
// BODY

void printimage(gray* image,int rows,int cols,int maxval){

    printf("P3\n");
    printf("%d %d \n", cols, rows);
    printf("%d\n",maxval);

    for(int i=0; i < rows; i++)
      for(int j=0; j < cols ; j++){
	printf("%u ",image[DPC*i * cols + DPC*j+0] );
	printf("%u ",image[DPC*i * cols + DPC*j+1] );
	printf("%u \n",image[DPC*i * cols + DPC*j+2] );
      }

}

gray* readimage(char* filepath)
    {
    FILE* ifp;
    gray* imagemap;
    int ich1, ich2;
    int validppm;
    //int j,i ;
    
    ifp = fopen(filepath,"r");
    if (ifp == NULL) {
      fprintf(stderr,"Error openning the file, check if you specified -i. Path: %s\n", filepath);
      exit(1);
    }

    /* Lecture du Magic number */
    ich1 = getc( ifp );
    if ( ich1 == EOF )
        pm_erreur( "EOF / erreur de lecture / nombre magique" );
    ich2 = getc( ifp );
    if ( ich2 == EOF )
        pm_erreur( "EOF / erreur de lecture / nombre magique" );
    //if(ich2 != '2' && ich2 != '5')
    //    pm_erreur(" mauvais type de fichier ");
    //else
      if(ich2 != '3'){
	fprintf(stderr,"Invalid file type");
        exit(1);
      }

    /* Lecture des dimensions */
    cols = pm_getint( ifp );
    rows = pm_getint( ifp );
    maxval = pm_getint( ifp );

    /* Allocation memoire  */
    imagemap = (gray *) malloc(3 * cols * rows * sizeof(gray));

    /* Lecture */
    for(int i=0; i < rows; i++){
      for(int j=0; j < cols ; j++){
	imagemap[DPC* i * cols + DPC*j + 0] = pm_getint(ifp);
	imagemap[DPC* i * cols + DPC*j + 1] = pm_getint(ifp);
	imagemap[DPC* i * cols + DPC*j + 2] = pm_getint(ifp);
	}
     }

      /* fermeture */
      fclose(ifp);
      return imagemap;
}


int main(int argc, char* argv[]){

	gray *graymap=readimage("image/clownplaintext.ppm");

	printimage(graymap,rows,cols,maxval);
}
