#include <stdlib.h>
#include <stdio.h>
#include "Util.h"


// HEADER
void printimage(double* image,int cols, int rows, int maxval);
double* readimage(char* filepath);
int rows=0,cols=0,maxval=0;
// BODY

void printimage(double* image,int rows,int cols,int maxval){

    printf("P3\n");
    printf("%d %d \n", cols, rows);
    printf("%d\n",maxval);

    for(int i=0; i < rows; i++)
      for(int j=0; j < cols ; j++){
	printf("%i\n",(int)image[i * cols + j+0] );
	printf("%i\n",(int)image[i * cols + j+1] );
	printf("%i\n",(int)image[i * cols + j+2] );
      }

}

double* readimage(char* filepath)
    {
    FILE* ifp;
    double* imagemap;
    int ich1, ich2;
    int validppm;
    int j,i ;
    
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
      if(ich2 == '3')
	validppm = 1;
      else validppm = 0;

    /* Lecture des dimensions */
    cols = pm_getint( ifp );
    rows = pm_getint( ifp );
    maxval = pm_getint( ifp );

    /* Allocation memoire  */
    imagemap = (double *) malloc(3 * cols * rows * sizeof(double));

    
    if(validppm){
    /* Lecture */
	    for(i=0; i < rows; i++){
	      for(j=0; j < cols ; j++){
		imagemap[i * cols + j + 0] = pm_getint(ifp);
		imagemap[i * cols + j + 1] = pm_getint(ifp);
		imagemap[i * cols + j + 2] = pm_getint(ifp);
		}
	    }

    }else { 
        fprintf(stderr,"Invalid file type");
        exit(1);
    }
      /* fermeture */
      fclose(ifp);
      return imagemap;
}


int main(int argc, char* argv[]){

	double *im=readimage("image/clownplaintext.ppm");

	printimage(im,rows,cols,maxval);

}

/*
int main(int argc, char* argv[]){
    FILE* ifp;
    gray* graymap;
    int ich1, ich2, rows, cols, maxval, pgmraw,j,i ;

    if ( argc != 2 ){
      printf("\nUsage : %s fichier \n\n", argv[0]);
      exit(0);
    }

    ifp = fopen(argv[1],"r");
    if (ifp == NULL) {
      printf("erreur d'ouverture du fichier %s\n", argv[1]);
      exit(1);
    }

    ich1 = getc( ifp );
    if ( ich1 == EOF )
        pm_erreur( "EOF / erreur de lecture / nombre magique" );
    ich2 = getc( ifp );
    if ( ich2 == EOF )
        pm_erreur( "EOF / erreur de lecture / nombre magique" );
    if(ich2 != '3'){// && ich2 != '5')
      pm_erreur(" mauvais type de fichier ");
      pgmraw = 0;
    }else
		pgmraw = 1;
    cols = pm_getint( ifp );
    rows = pm_getint( ifp );
    maxval = pm_getint( ifp );

    graymap = (gray *) malloc(3 * cols * rows * sizeof(gray));
    
    if(pgmraw){
	    for(i=0; i < rows; i++)
	      for(j=0; j < cols ; j++){
		graymap[i * cols + j + 0] = pm_getint(ifp);
		graymap[i * cols + j + 1] = pm_getint(ifp);
		graymap[i * cols + j + 2] = pm_getint(ifp);
		}

      	    printf("P3\n");

	    printf("%d %d \n", cols, rows);
	    printf("%d\n",maxval);
	   
	    for(i=0; i < rows; i++)
	      for(j=0; j < cols ; j++)
		printf("%i ",(graymap[i * cols + j+0]+graymap[i * cols + j+1]+graymap[i * cols + j+2])/3 );
		//printf("%u ",graymap[i * cols + j] );

    }else {
	    printf("ERROR!\n");
    }

      fclose(ifp);
      return 0;
}

*/
