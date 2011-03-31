#include <stdlib.h>
#include <stdio.h>
#include "showregion.h"
#define DPC 3 //data per cell information

void printimage(pimage_type image){

    printf("P3\n");
    printf("%d %d \n", image->cols, image->rows);
    printf("%d\n",image->maxval);

    for(int i=0; i < image->rows; i++)
      for(int j=0; j < image->cols ; j++){
	printf("%u ",image->stream[DPC*i * image->cols + DPC*j+0] );
	printf("%u ",image->stream[DPC*i * image->cols + DPC*j+1] );
	printf("%u \n",image->stream[DPC*i * image->cols + DPC*j+2] );
      }

}

pimage_type readimage(char* filepath)
    {
    FILE* ifp;
    gray* imagemap;
    int ich1, ich2;
    pimage_type image=(pimage_type)malloc(sizeof(image_type));  
 
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
    image->cols = pm_getint( ifp );
    image->rows = pm_getint( ifp );
    image->maxval = pm_getint( ifp );

    /* Allocation memoire  */
    imagemap = (gray *) malloc(3 * image->cols * image->rows * sizeof(gray));

    /* Lecture */
    for(int i=0; i < image->rows; i++){
      for(int j=0; j < image->cols ; j++){
	imagemap[DPC* i * image->cols + DPC*j + 0] = pm_getint(ifp);
	imagemap[DPC* i * image->cols + DPC*j + 1] = pm_getint(ifp);
	imagemap[DPC* i * image->cols + DPC*j + 2] = pm_getint(ifp);
	}
     }

     image->stream=imagemap;

      /* fermeture */
      fclose(ifp);
      return image; //imagemap;
}


int main(int argc, char* argv[]){

	pimage_type graymap=readimage("image/clownplaintext.ppm");

	printimage(graymap);
}
