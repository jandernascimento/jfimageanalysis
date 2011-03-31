#include <stdlib.h>
#include <stdio.h>
#include "showregion.h"

//Finds the position in between the two pixels (euclidian in 3d space)
//DO ME
ppixel_type pixel_middle(ppixel_type p1,ppixel_type p2){

  ppixel_type p;

  return p;
}

//Assign the color of the pixel as the same color as the closest key (list of k's)
//DO ME
void assign_to_group(pimage_type image, ppixel_type keys, int size){

}

//Calculates what is the middle pixel (final k)
void find_groups(pimage_type image,ppixel_type pixels,int size){

	//pixel_type pixels_k[size];

	for(int y;y<image->rows;y++){

		for(int x;x<image->cols;x++){
			
			for(int j=0;j<size;j++){

				ppixel_type k=pixels;
				ppixel_type pix=get_pixel(image,x,y);
			
				//pixels[j]=pixel_middle(k,pix);			

			}
		}

	}



}

//Computes the distance between two pixels with euclidian
//DO ME
int pixel_distance(pimage_type image,pixel_type p1,pixel_type p2){

  return 0;

}

ppixel_type get_pixel(pimage_type image,int x, int y){

  pixel_type pixel;
  pixel.r=image->stream[DPC*y * image->cols + DPC*x+0];
  pixel.g=image->stream[DPC*y * image->cols + DPC*x+1];
  pixel.b=image->stream[DPC*y * image->cols + DPC*x+2];
  pixel.x=x;
  pixel.y=y;
  return &pixel;

}


void set_pixel(pimage_type image,int x, int y, pixel_type pixel ){

  image->stream[DPC*y * image->cols + DPC*x+0]=pixel.r;
  image->stream[DPC*y * image->cols + DPC*x+1]=pixel.g;
  image->stream[DPC*y * image->cols + DPC*x+2]=pixel.b;


}

void printimage(pimage_type image){

    printf("P3\n");
    printf("%d %d \n", image->cols, image->rows);
    printf("%d\n",image->maxval);

    for(int i=0; i < image->rows; i++)
      for(int j=0; j < image->cols ; j++){
	//printf("%u ",image->stream[DPC*i * image->cols + DPC*j+0] );
	//printf("%u ",image->stream[DPC*i * image->cols + DPC*j+1] );
	//printf("%u \n",image->stream[DPC*i * image->cols + DPC*j+2] );
	printf("%u ",get_pixel(image,j,i)->r);
	printf("%u ",get_pixel(image,j,i)->g);
	printf("%u \n",get_pixel(image,j,i)->b);
        
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
    image->stream=imagemap;

    /* Lecture */
    for(int i=0; i < image->rows; i++){
      for(int j=0; j < image->cols ; j++){
	//imagemap[DPC* i * image->cols + DPC*j + 0] = pm_getint(ifp);
	//imagemap[DPC* i * image->cols + DPC*j + 1] = pm_getint(ifp);
	//imagemap[DPC* i * image->cols + DPC*j + 2] = pm_getint(ifp);
	pixel_type pixel;
	pixel.r=pm_getint(ifp);
	pixel.g=pm_getint(ifp);
	pixel.b=pm_getint(ifp);
	set_pixel(image,j,i,pixel);
	}
     }


      /* fermeture */
      fclose(ifp);
      return image; //imagemap;
}


int main(int argc, char* argv[]){

	pimage_type graymap=readimage("image/clownplaintext.ppm");

	printimage(graymap);
}
