#include <stdlib.h>
#include <stdio.h>
#include "Util.h"



int main(int argc, char* argv[])
    {
    FILE* ifp;
    bit* bitmap;
    int ich1, ich2, rows, cols ;
    int i,j ;


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

    /* Lecture et vérification du Magic number */
    ich1 = getc( ifp );
    if ( ich1 == EOF )
        pm_erreur( "EOF / erreur de lecture / nombre magique" );
    ich2 = getc( ifp );
    if ( ich2 == EOF )
        pm_erreur( "EOF / erreur de lecture / nombre magique" );
    if(ich2 != '1')
      pm_erreur(" mauvais type de fichier ");

    /* Lecture des dimensions */
    cols = pm_getint( ifp );
    rows = pm_getint( ifp );

    /* Allocation memoire  */
    bitmap = (bit *) malloc(cols * rows * sizeof(bit));

    /* Lecture */
    for(i=0; i < rows; i++)
      for(j=0; j < cols ; j++)
	bitmap[i * cols + j] = pm_getbit(ifp);

    /* Ecriture */
   printf("P2\n");
   printf("%d %d \n", cols, rows);
   printf("1\n");
   for(i=0; i < rows; i++)
      for(j=0; j < cols ; j++)
	printf("%u ",bitmap[i * cols + j] );



      /* fermeture */
      fclose(ifp);
      return 0;
}
