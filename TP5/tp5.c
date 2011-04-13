#include <stdlib.h>
#include <stdio.h>
#include "tp5.h"

pixel_type *get_pixel(pimage_type image,int x, int y){
  pixel_type *pixel=(pixel_type*)malloc(sizeof(pixel_type));
  
  //pixel->r = image->stream[DPC*y * image->cols + DPC*x+RED];
  //pixel->g = image->stream[DPC*y * image->cols + DPC*x+GREEN];
  pixel->g = image->stream[DPC*y * image->cols + DPC*x];
  //pixel->b = image->stream[DPC*y * image->cols + DPC*x+BLUE];
  //pixel->k = image->stream[DPC*y * image->cols + DPC*x+K];
  pixel->x = x;
  pixel->y = y;
  
  return pixel;
}


void set_pixel(pimage_type image,int x, int y, pixel_type pixel ){
  //image->stream[DPC*y * image->cols + DPC*x+RED]   = pixel.r;
  //image->stream[DPC*y * image->cols + DPC*x+GREEN] = pixel.g;
  image->stream[DPC*y * image->cols + DPC*x] = pixel.g;
  //image->stream[DPC*y * image->cols + DPC*x+BLUE]  = pixel.b;
}

void printimage(pimage_type image){

    printf("P5\n");
    printf("%d %d \n", image->cols, image->rows);
    printf("%d\n",image->maxval);

    for(int i=0; i < image->rows; i++){
      for(int j=0; j < image->cols ; j++){
		//printf("%u ",get_pixel(image,j,i)->r);
		//printf("%u ",get_pixel(image,j,i)->g);
		printf("%c",get_pixel(image,j,i)->g);
		//printf("%u \n",get_pixel(image,j,i)->b);
      }
      //printf("\n");
    }
}

pimage_type readimage(char* filepath){
    FILE* ifp;
    gray* imagemap;
    int ich1, ich2;
    pimage_type image=(pimage_type)malloc(sizeof(image_type));  
 
    ifp = fopen(filepath,"r");
    if (ifp == NULL) {
      fprintf(stderr,"Error openning the file, check if you specified -i. Path: %s\n", filepath);
      exit(1);
    }

    image->path=filepath;

    /* Lecture du Magic number */
    ich1 = getc( ifp );
    if ( ich1 == EOF )
        pm_erreur( "EOF / erreur de lecture / nombre magique" );
    ich2 = getc( ifp );
    if ( ich2 == EOF )
        pm_erreur( "EOF / erreur de lecture / nombre magique" );
    if(ich2 != '5'){
    	fprintf(stderr,"Invalid file type");
        exit(1);
    }

    /* Lecture des dimensions */
    image->cols = pm_getint( ifp );
    image->rows = pm_getint( ifp );
    image->maxval = pm_getint( ifp );

    /* Allocation memoire  */
    imagemap = (gray *) malloc(DPC * image->cols * image->rows * sizeof(gray));
    image->stream=imagemap;

    /* Lecture */
    for(int i=0; i < image->rows; i++){
      for(int j=0; j < image->cols ; j++){
		pixel_type pixel;
		//pixel.r=pm_getint(ifp);
		//pixel.g=pm_getint(ifp);
		pixel.g=pm_getrawbyte(ifp);
		//pixel.b=pm_getint(ifp);
		set_pixel(image,j,i,pixel);
	   }
     }

      /* fermeture */
      fclose(ifp);
      return image; //imagemap;
}

/** Parser method **

char* getStrParam(int argc,char* argv[],char* param,char* def){

	for(int i=0;i<argc;i++){
		if(strcmp(param,argv[i])==0){ 
			if(argv[i+1][0]=='-') return def;
			return argv[i+1];
		}
	}

	return def;

}

int getBoolParam(int argc,char* argv[],char* param){

	for(int i=0;i<argc;i++){
		if(strcmp(param,argv[i])==0){ 
			return 1;
		}
	}

	return 0;

}

int getIntParam(int argc,char* argv[],char* param,char* def){
    char* res=getStrParam(argc,argv,param,def);
    
    int val=0;
 
    for(int i=0;i<strlen(res);i++){
      	if(res[0]=='-') continue;
	val+=pow(10,(strlen(res)-i-1))*(res[i]-'0'); 
    }

    return res[0]=='-'?-1*val:val;
}
** end Parser method **/

int main(int argc, char* argv[]){

	pimage_type image=readimage("image/taz/taz001.pgm");

	printimage(image);
}
