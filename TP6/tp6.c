#include <stdlib.h>
#include <stdio.h>
#include "tp6.h"

/**
** read background images
**/
filelist_type readbackgrounds(char *path, int n){
	filelist_type filelist;
	filelist.size = n;
	filelist.paths=(char*)malloc(sizeof(char*)*n);
	char *i_str=(char*)malloc(sizeof(char)*4);	

	for(int i=0;i<=n;i++){ 					
		filelist.paths[i]=(char*)malloc(sizeof(char)*70);

		strcpy(filelist.paths[i], path);
		strcat(filelist.paths[i],"000\0"); //we always concatenate with these zeros
		if (i<10)
			strcat(filelist.paths[i],"00\0");
		else if (i>=10 && i<100)
			strcat(filelist.paths[i],"0\0");
		
		sprintf(i_str,"%i\0",i);//converting from int to string		
		strcat(filelist.paths[i],"\0");
		strcat(filelist.paths[i],i_str);
		strcat(filelist.paths[i],".ppm\0");	
	}

	/*/testing the structure
	FILE* ifp;
	for(int i=0;i<=n;i++){ 
		fprintf(stderr,"%s\n", filelist.paths[i]);
    	ifp = fopen(filelist.paths[i],"r");
    	if (ifp == NULL) {
    	  fprintf(stderr,"Error openning the file, check if you specified -i. Path: %s\n", filelist.paths[i]);
    	  exit(1);
    	}		
    }
    fclose(ifp);
    */

	free(i_str);
	return filelist;
}

/**
** calc mean image
**/
pimage_type calculate_mean_image(filelist_type list){

}

/**
** calc subtraction
**/
gray *matrix_subtraction(gray *matrix1,gray *matrix2,int rows,int cols){

}

/**
** get pixel from a matrix
**/
gray *get_matrix_pixel(gray *matrix,int row, int col,int dim){

	return matrix[row*dim+col];

}

/**
** calc determinant
**/
int matrix_determinant(gray *matrix,int dim){

	gray a=get_matrix_pixel(matrix,0,0,dim);
	gray b=get_matrix_pixel(matrix,0,1,dim);
	gray c=get_matrix_pixel(matrix,0,2,dim);

	gray d=get_matrix_pixel(matrix,1,0,dim);
	gray e=get_matrix_pixel(matrix,1,1,dim);
	gray f=get_matrix_pixel(matrix,1,2,dim);

	gray g=get_matrix_pixel(matrix,2,0,dim);
	gray h=get_matrix_pixel(matrix,2,1,dim);
	gray i=get_matrix_pixel(matrix,2,2,dim);

	return (a*e*i + b*f*g + c*d*h - a*f*h - b*d*i - c*e*g);

}

void matrix_determinant_test(){
	const int dim=3;
	gray *det=(gray *)malloc(sizeof(gray)*dim);

	det[dim*0+0]=1;
	det[dim*0+1]=2;
	det[dim*0+2]=1;

	det[dim*1+0]=3;
	det[dim*1+1]=4;
	det[dim*1+2]=5;

	det[dim*2+0]=5;
	det[dim*2+1]=6;
	det[dim*2+2]=7;

	int val=matrix_determinant(det,dim);

	free(det);

	printf("-->%i\n",val);
}

/**
** Get pixel information from a image type
**/
pixel_type *get_pixel(pimage_type image,int x, int y){
  pixel_type *pixel=(pixel_type*)malloc(sizeof(pixel_type));
  
  pixel->r = image->stream[DPC*y * image->cols + DPC*x+RED];
  pixel->g = image->stream[DPC*y * image->cols + DPC*x+GREEN];
  pixel->b = image->stream[DPC*y * image->cols + DPC*x+BLUE];
//  pixel->k = image->stream[DPC*y * image->cols + DPC*x+K];
  pixel->x = x;
  pixel->y = y;
  
  return pixel;
}


/**
** Set pixel information from a image type
**/
void set_pixel(pimage_type image,int x, int y, pixel_type pixel ){
  image->stream[DPC*y * image->cols + DPC*x+RED]   = pixel.r;
  image->stream[DPC*y * image->cols + DPC*x+GREEN] = pixel.g;
  image->stream[DPC*y * image->cols + DPC*x+BLUE]  = pixel.b;
}

/**
** Print image, take the type image as input
**/
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

/**
** Read image from a file 
**/
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
      return image;
}

/** Parser method **/

/**
** Get string param
**/
char* getStrParam(int argc,char* argv[],char* param,char* def){

	for(int i=0;i<argc;i++){
		if(strcmp(param,argv[i])==0){ 
			if(argv[i+1][0]=='-') return def;
			return argv[i+1];
		}
	}

	return def;

}

/**
** Get bool param
**/
int getBoolParam(int argc,char* argv[],char* param){

	for(int i=0;i<argc;i++){
		if(strcmp(param,argv[i])==0){ 
			return 1;
		}
	}

	return 0;

}

/**
** Get int param
**/
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

/**
** Main method
**/
int main(int argc, char* argv[]){

	//readbackgrounds("background_substraction/background/img_", 10);

	pimage_type image=readimage("background_substraction/img_000053.ppm");

	printimage(image);

	/*char *filepath=getStrParam(argc,argv,"-i","");
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

	printimage(image);*/

}
