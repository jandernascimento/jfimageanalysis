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
    if(ich2 != '2' && ich2 != '5')
      pm_erreur(" mauvais type de fichier ");
    else
      if(ich2 == '2')
	pgmraw = 0;
      else pgmraw = 1;

    /* Lecture des dimensions */
    cols = pm_getint( ifp );
    rows = pm_getint( ifp );
    maxval = pm_getint( ifp );

    /* Allocation memoire  */
    graymap = (gray *) malloc(cols * rows * sizeof(gray));

    
    if(pgmraw){
    /* Convertendo para Plaintext */
    /* Lecture */
	    for(i=0; i < rows; i++)
	      for(j=0; j < cols ; j++)
		graymap[i * cols + j] = pm_getrawbyte(ifp);

      	    printf("P2\n");

	    printf("%d %d \n", cols, rows);
	    printf("%d\n",maxval);
	   
	    for(i=0; i < rows; i++)
	      for(j=0; j < cols ; j++)
		printf("%u ",graymap[i * cols + j] );

    }else {
    /* Convertendo para Binario */
    /* Lecture */
    for(i=0; i < rows; i++)
      for(j=0; j < cols ; j++)
	graymap[i * cols + j] = pm_getint(ifp);

     /* Ecriture */    
     printf("P5\n");

	printf("%d %d \n", cols, rows);
	    printf("%d\n",maxval);
	   
	    for(i=0; i < rows; i++)
	      for(j=0; j < cols ; j++)
		printf("%i",graymap[i * cols + j] );

    }

    


      /* fermeture */
      fclose(ifp);
      return 0;
}
