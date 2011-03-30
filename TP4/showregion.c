#include <stdlib.h>
#include <stdio.h>
#include "Util.h"



int main(int argc, char* argv[])
    {
    FILE* ifp;
    gray* graymap;
    int ich1, ich2, rows, cols, maxval, pgmraw,j,i ;



    /* Test des arguments */
    if ( argc != 2 ){
      printf("\nUsage : %s fichier \n\n", argv[0]);
      exit(0);
    }

    /* Ouverture */
    ifp = fopen(argv[1],"r");
    if (ifp == NULL) {
      printf("erreur d'ouverture du fichier %s\n", argv[1]);
      exit(1);
    }

    /* Lecture du Magic number */
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
    //  if(ich2 == '2')
		pgmraw = 1;
      //else pgmraw = 1;

    /* Lecture des dimensions */
    cols = pm_getint( ifp );
    rows = pm_getint( ifp );
    maxval = pm_getint( ifp );

    /* Allocation memoire  */
    graymap = (gray *) malloc(3 * cols * rows * sizeof(gray));

    
    if(pgmraw){
    /* Convertendo para PGM */
	    for(i=0; i < rows; i++)
	      for(j=0; j < cols ; j++){
		graymap[i * cols + j + 0] = pm_getint(ifp);
		graymap[i * cols + j + 1] = pm_getint(ifp);
		graymap[i * cols + j + 2] = pm_getint(ifp);
		}

      	    printf("P2\n");

	    printf("%d %d \n", cols, rows);
	    printf("%d\n",maxval);
	   
	    for(i=0; i < rows; i++)
	      for(j=0; j < cols ; j++)
		printf("%i ",(graymap[i * cols + j+0]+graymap[i * cols + j+1]+graymap[i * cols + j+2])/3 );
		//printf("%u ",graymap[i * cols + j] );

    }else {
    /* Convertendo para Binario */
	    printf("ERROR!\n");
    }

    


      /* fermeture */
      fclose(ifp);
      return 0;
}
