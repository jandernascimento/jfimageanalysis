#include <stdlib.h>
#include <stdio.h>
#include "showregion.h"

//Finds the position in between the two pixels (euclidian in 3d space)
//Reminder: the 3d space is the color
pixel_type *pixel_random_group(pimage_type image,int size){

	pixel_type* group=(pixel_type *)malloc(size*sizeof(pixel_type));

	for(int x=0;x<size;x++){

		int xpos=(rand()%255);
		int ypos=(rand()%255);
		(group+sizeof(pixel_type)*x)->x=xpos;
		(group+sizeof(pixel_type)*x)->y=ypos;
		(group+sizeof(pixel_type)*x)->k=x;
		(group+sizeof(pixel_type)*x)->r=get_pixel(image,xpos,ypos)->r;
		(group+sizeof(pixel_type)*x)->g=get_pixel(image,xpos,ypos)->g;
		(group+sizeof(pixel_type)*x)->b=get_pixel(image,xpos,ypos)->b;
		
	//fprintf(stderr,"chosen %u,%u,%u\n",(group+sizeof(pixel_type)*x)->x,(group+sizeof(pixel_type)*x)->y,(group+sizeof(pixel_type)*x)->k);	

	}

	/*pixel_type *chosen=get_pixel(image,0,0);
	(group+sizeof(pixel_type)*0)->x=chosen->x;
	(group+sizeof(pixel_type)*0)->y=chosen->y;
	(group+sizeof(pixel_type)*0)->r=chosen->r;
	(group+sizeof(pixel_type)*0)->g=chosen->g;
	(group+sizeof(pixel_type)*0)->b=chosen->b;
	(group+sizeof(pixel_type)*0)->k=0;
	
	pixel_type *chosen2=get_pixel(image,255,255);
	(group+sizeof(pixel_type)*1)->x=chosen2->x;
	(group+sizeof(pixel_type)*1)->y=chosen2->y;
	(group+sizeof(pixel_type)*1)->r=chosen2->r;
	(group+sizeof(pixel_type)*1)->g=chosen2->g;
	(group+sizeof(pixel_type)*1)->b=chosen2->b;
	(group+sizeof(pixel_type)*1)->k=1;*/
	
	
	return group;	
}



pixel_type *pixel_middle(pixel_type *p1,pixel_type *p2){

  pixel_type *p;

  return p;
}

//Assign the color of the pixel as the same color as the closest key (list of k's)
pimage_type assign_to_group(pimage_type image, pixel_type *keys, int size){
	pimage_type image2=readimage(image->path);

	//associating each pixel to a group
	for(int y=0;y<image->rows;y++){
		for(int x=0;x<image->cols;x++){
			pixel_type *chosen;
			pixel_type *taken;
			int lowest_distance=999999999;
			pixel_type *image_pixel=get_pixel(image,x,y);			

			//Search in all groups in each one that a single pixel belongs to
			int nb_group=0;
			for(int j=0;j<size;j++){
				nb_group=j;
				chosen=(keys+sizeof(pixel_type)*j);
				int di=pixel_distance(image,chosen,image_pixel);			
				if(di<lowest_distance){ 
					taken=(keys+sizeof(pixel_type)*j); 
					lowest_distance=di; 
				}
			}
			//setting the group in the pixel
			//set_group(image2,x,y,nb_group);
			 
			//pixel_type *tp=get_pixel(image,(*taken).x,(*taken).y);	
			//set_pixel(image2,x,y,*tp);
			set_pixel(image2,x,y,*taken);
		}
	}
	//recalculating the centroids (middle of the groups)
	

	return image2;
}

int recalc_centroids(pimage_type image, int nro_groups){
	int something_changed = 0; //false
	
	int * sum_r_per_group        = (int*) malloc(sizeof(int) * nro_groups);
	int * sum_g_per_group        = (int*) malloc(sizeof(int) * nro_groups);
	int * sum_b_per_group        = (int*) malloc(sizeof(int) * nro_groups);
	int * count_pixels_per_group = (int*) malloc(sizeof(int) * nro_groups);

	initialize_array(sum_r_per_group,nro_groups);
	initialize_array(sum_g_per_group,nro_groups);
	initialize_array(sum_b_per_group,nro_groups);
	initialize_array(count_pixels_per_group,nro_groups);
	
	for(int y=0;y<image->rows;y++){
		for(int x=0;x<image->cols;x++){
			pixel_type *image_pixel=get_pixel(image,x,y);

			sum_r_per_group       [ image_pixel->k ] = sum_r_per_group[ image_pixel->k ] + image_pixel->r;
			sum_g_per_group       [ image_pixel->k ] = sum_g_per_group[ image_pixel->k ] + image_pixel->g;
			sum_b_per_group       [ image_pixel->k ] = sum_b_per_group[ image_pixel->k ] + image_pixel->b;
			count_pixels_per_group[ image_pixel->k ] = count_pixels_per_group[ image_pixel->k ] ++;
		}
	}

	for(int j=0;j<nro_groups;j++){
		int new_value_r = sum_r_per_group[j] / count_pixels_per_group[j]; 
		int new_value_g = sum_g_per_group[j] / count_pixels_per_group[j];
		int new_value_b = sum_b_per_group[j] / count_pixels_per_group[j];

		//update_group(group, j, chosen->r,chosen->g,chosen->b);
	}


	return something_changed;
}

void initialize_array(int * vec, int n){
	for(int i=0;i<n;i++)
		vec[i]=0;
}

//put the group of the pixel in the image
void set_group(pimage_type image,int x, int y,int group){
	image->stream[DPC*y * image->cols + DPC*x+K]=group;
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

  return sqr;
}

pixel_type *get_pixel(pimage_type image,int x, int y){
  pixel_type *pixel=(pixel_type*)malloc(sizeof(pixel_type));
  
  pixel->r = image->stream[DPC*y * image->cols + DPC*x+RED];
  pixel->g = image->stream[DPC*y * image->cols + DPC*x+GREEN];
  pixel->b = image->stream[DPC*y * image->cols + DPC*x+BLUE];
  pixel->k = image->stream[DPC*y * image->cols + DPC*x+K];
  pixel->x = x;
  pixel->y = y;
  
  return pixel;
}


void set_pixel(pimage_type image,int x, int y, pixel_type pixel ){
  image->stream[DPC*y * image->cols + DPC*x+RED]   = pixel.r;
  image->stream[DPC*y * image->cols + DPC*x+GREEN] = pixel.g;
  image->stream[DPC*y * image->cols + DPC*x+BLUE]  = pixel.b;
}

void printimage(pimage_type image){

    printf("P3\n");
    printf("%d %d \n", image->cols, image->rows);
    printf("%d\n",image->maxval);

    for(int i=0; i < image->rows; i++)
      for(int j=0; j < image->cols ; j++){
		printf("%u ",get_pixel(image,j,i)->r);
		printf("%u ",get_pixel(image,j,i)->g);
		printf("%u \n",get_pixel(image,j,i)->b);
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
    if(ich2 != '3'){
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

void update_group(pixel_type *group, int index, int r, int g, int b){
	(group+sizeof(pixel_type)*index)->x=0;
	(group+sizeof(pixel_type)*index)->y=0;
	(group+sizeof(pixel_type)*index)->r=r;
	(group+sizeof(pixel_type)*index)->g=g;
	(group+sizeof(pixel_type)*index)->b=b;
	(group+sizeof(pixel_type)*index)->k=index;
}


int main(int argc, char* argv[]){

	pimage_type image=readimage("image/jander.ppm");

	int nro_groups = 2;

	pixel_type* group=(pixel_type *)malloc(nro_groups*sizeof(pixel_type));
	
	//group1	
	pixel_type *chosen=get_pixel(image,79,96);
	update_group(group, 0, chosen->r,chosen->g,chosen->b);

	//group2
	pixel_type *chosen2=get_pixel(image,0,0);
	update_group(group, 1, chosen2->r,chosen2->g,chosen2->b);

	pimage_type image_grouped=assign_to_group(image, group,nro_groups);
	
	printimage(image_grouped);
}
