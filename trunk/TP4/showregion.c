#include <stdlib.h>
#include <stdio.h>
#include "showregion.h"

//generates groups randomly
pixel_type *pixel_random_group(pimage_type image,int size){

	pixel_type* group=(pixel_type *)malloc(size*sizeof(pixel_type));

	for(int m=0;m<size;m++){
		srand(m+78887777);
		int xpos=(int)(rand()%255);
		srand(m+1777);
		int ypos=(int)(rand()%255);
		
		(group+sizeof(pixel_type)*m)->x=xpos;
		(group+sizeof(pixel_type)*m)->y=ypos;
		(group+sizeof(pixel_type)*m)->k=m;
		(group+sizeof(pixel_type)*m)->r=get_pixel(image,xpos,ypos)->r;
		(group+sizeof(pixel_type)*m)->g=get_pixel(image,xpos,ypos)->g;
		(group+sizeof(pixel_type)*m)->b=get_pixel(image,xpos,ypos)->b;		
	}
	return group;	
}


//Assign the color of the pixel as the same color as the closest key (list of k's)
pimage_type assign_to_group(pimage_type image, pixel_type *keys, int size){
	pimage_type image2=readimage(image->path);

	int something_changed = 1; //true
	
	while (something_changed) {
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
					chosen=(keys+sizeof(pixel_type)*j);
					int di=pixel_distance(image,chosen,image_pixel);			
					if(di<lowest_distance){ 
						taken=(keys+sizeof(pixel_type)*j); 
						lowest_distance=di; 
						nb_group=j;
					}
				}
				//setting the group in the pixel
				set_group(image2,x,y,nb_group);
			}
		}
		//recalculating the centroids (middle of the groups)
		something_changed=recalc_centroids(image2,keys,size);
	}
	
	//setting the collor of the centroid for all the pixels of that group
	update_color_pixels(image2,keys);
	return image2;
}

void update_color_pixels(pimage_type image, pixel_type *keys){
	pixel_type *pixel_key;
	
	for(int y=0;y<image->rows;y++)
		for(int x=0;x<image->cols;x++){
			int index_group=get_group(image,x,y);
			pixel_key=(keys+sizeof(pixel_type) * index_group);
			set_pixel(image,x,y,*pixel_key);
		}
}

int recalc_centroids(pimage_type image, pixel_type *keys, int nro_groups){
	int something_changed = 0; //false

	int * sum_r_per_group        = (int*) malloc(sizeof(int) * nro_groups);
	int * sum_g_per_group        = (int*) malloc(sizeof(int) * nro_groups);
	int * sum_b_per_group        = (int*) malloc(sizeof(int) * nro_groups);
	int * count_pixels_per_group = (int*) malloc(sizeof(int) * nro_groups);

	initialize_array(sum_r_per_group,nro_groups,0);
	initialize_array(sum_g_per_group,nro_groups,0);
	initialize_array(sum_b_per_group,nro_groups,0);
	initialize_array(count_pixels_per_group,nro_groups,1);
	
	for(int y=0;y<image->rows;y++){
		for(int x=0;x<image->cols;x++){
			pixel_type *image_pixel=get_pixel(image,x,y);

			sum_r_per_group       [ image_pixel->k ] = sum_r_per_group[ image_pixel->k ] + image_pixel->r;
			sum_g_per_group       [ image_pixel->k ] = sum_g_per_group[ image_pixel->k ] + image_pixel->g;
			sum_b_per_group       [ image_pixel->k ] = sum_b_per_group[ image_pixel->k ] + image_pixel->b;
			count_pixels_per_group[ image_pixel->k ] = count_pixels_per_group[ image_pixel->k ] + 1;
		}
	}


	pixel_type *centroid;
	for(int j=0;j<nro_groups;j++){
		
		int new_value_r = sum_r_per_group[j] / count_pixels_per_group[j]; 
		int new_value_g = sum_g_per_group[j] / count_pixels_per_group[j];
		int new_value_b = sum_b_per_group[j] / count_pixels_per_group[j];

		//checking if the centroid of the group changed
		centroid=(keys+sizeof(pixel_type)*j);
		if ((centroid->r != new_value_r) || (centroid->g != new_value_g) || (centroid->b != new_value_b)){
			something_changed=1; //true
			update_group(keys, j, new_value_r, new_value_g, new_value_b);
		}
	}


	return something_changed;
}

void initialize_array(int * vec, int n, int vldef){
	for(int i=0;i<n;i++)
		vec[i]=vldef;
}

//put the group of the pixel in the image
void set_group(pimage_type image,int x, int y,int group){
	image->stream[DPC*y * image->cols + DPC*x+K]=group;
}

int get_group(pimage_type image,int x, int y){
	return image->stream[DPC*y * image->cols + DPC*x+K];
}


//Computes the distance between two pixels with euclidian
//Reminder: the 3d space is the color
int pixel_distance(pimage_type image,pixel_type *p1,pixel_type *p2){

  int value_r=(int)pow(p1->r-p2->r,2);
  int value_g=(int)pow(p1->g-p2->g,2);
  int value_b=(int)pow(p1->b-p2->b,2);

  int value_x=(int)pow(p1->x-p2->x,2);
  int value_y=(int)pow(p1->y-p2->y,2);

  int sqr=(int)sqrt(value_r+value_g+value_b+value_x+value_y);

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

/** Parser method **/

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
/** end Parser method **/

int main(int argc, char* argv[]){

	char *filepath=getStrParam(argc,argv,"-i","");
	int nro_groups=getIntParam(argc,argv,"-g","2");
	int ishelp=getBoolParam(argc,argv,"--help");

	if(ishelp){
		printf("Usage: showregion [OPTIONS]\n");	
		fprintf(stderr,"OPTIONS: \n");	
		fprintf(stderr,"-i: input file\n");	
		fprintf(stderr,"-g: number of groups, default is 2\n");
		exit(0);	
	}


	pimage_type image=readimage(filepath);

	//generating the first groups randomly
	pixel_type* group=pixel_random_group(image,nro_groups);

	//k-means algorith
	pimage_type image_grouped=assign_to_group(image, group,nro_groups);

	//printing the result
	printimage(image_grouped);
}
