#include <stdlib.h>
#include <stdio.h>
#include "showregion.h"

//Finds the position in between the two pixels (euclidian in 3d space)
//Reminder: the 3d space is the color
//DO ME
pixel_type *pixel_middle(pixel_type *p1,pixel_type *p2){

  pixel_type *p;

  return p;
}

//Assign the color of the pixel as the same color as the closest key (list of k's)
//DO ME
void assign_to_group(pimage_type image, pixel_type *keys, int size){

	//fprintf(stderr,"got here\n");
	for(int y=0;y<image->rows;y++){
		//fprintf(stderr,"got here\n");
		for(int x=0;x<image->cols;x++){
			fprintf(stderr,"begin:checking pixel\n");
			//pixel_type *k=keys;
			pixel_type *chosen=keys;
			pixel_type *taken=keys;
			int lowest_distance=999999999;
			fprintf(stderr,"\tbegin:checking distance\n");
			for(int j=0;j<size;j++){
				pixel_type *image_pixel=get_pixel(image,x,y);
				
				chosen=(keys+j);
				fprintf(stderr,"\tnew->chosen %i,%i\n",chosen->x,chosen->y);
				int di=pixel_distance(image,chosen,image_pixel);			
				if(di<lowest_distance) taken=(keys+j); 
			}
			fprintf(stderr,"\tend:checking distance\n");
  			
			fprintf(stderr,"\t->chosen %i,%i\n",(*chosen).x,(*chosen).y);
			pixel_type *tp=get_pixel(image,(*taken).x,(*taken).y);	
			fprintf(stderr,"\t->result setting (%i,%i,%i)\n",tp->r,tp->g,tp->b);
			fprintf(stderr,"end:checking pixel\n");
			set_pixel(image,x,y,tp);
		}

	}
}

//Calculates what is the middle pixel (final k)
//Reminder: the 3d space is the color
void find_groups(pimage_type image,pixel_type *pixels,int size){

	for(int y;y<image->rows;y++){

		for(int x;x<image->cols;x++){
			pixel_type *k=pixels;
			for(int j=0;j<size;j++){

				pixel_type *pix=get_pixel(image,x,y);
			
				pixels=pixel_middle(k,pix);			

			}
		}

	}

}
/*
void split_pixel_into_groups(pimage image, pixel_type *groups, int size){

}
*/

//Computes the distance between two pixels with euclidian
//DO ME
//Reminder: the 3d space is the color
int pixel_distance(pimage_type image,pixel_type *p1,pixel_type *p2){


  int value_r=(int)pow(p1->r-p2->r,2);
  int value_g=(int)pow(p1->g-p2->g,2);
  int value_b=(int)pow(p1->b-p2->b,2);

  int sqr=(int)sqrt(value_r+value_g+value_b);
  fprintf(stderr,"\t\tpixel p1(at %i %i color %i,%i,%i) to p2(at %i %i color %i,%i,%i) is %i\n",p1->x,p1->y,p1->r,p1->g,p1->b,p2->x,p2->y,p2->r,p2->g,p2->b,sqr);
  return sqr;

}

pixel_type *get_pixel(pimage_type image,int x, int y){

  pixel_type *pixel=(pixel_type*)malloc(sizeof(pixel_type));
  pixel->r=image->stream[DPC*y * image->cols + DPC*x+0];
  pixel->g=image->stream[DPC*y * image->cols + DPC*x+1];
  pixel->b=image->stream[DPC*y * image->cols + DPC*x+2];
  pixel->x=x;
  pixel->y=y;
  return pixel;

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

	pimage_type image=readimage("image/clownplaintext.ppm");

	pixel_type* group=(pixel_type *)malloc(2*sizeof(pixel_type));

	(group+sizeof(pixel_type)*1)->x=0;
	(group+sizeof(pixel_type)*1)->y=0;
	(group+sizeof(pixel_type)*1)->r=200;
	(group+sizeof(pixel_type)*1)->g=0;
	(group+sizeof(pixel_type)*1)->b=0;


	(group+sizeof(pixel_type)*0)->x=100;
	(group+sizeof(pixel_type)*0)->y=100;
	(group+sizeof(pixel_type)*0)->r=0;
	(group+sizeof(pixel_type)*0)->g=0;
	(group+sizeof(pixel_type)*0)->b=0;

	
	assign_to_group(image, group,2);
	
	printimage(image);
}
