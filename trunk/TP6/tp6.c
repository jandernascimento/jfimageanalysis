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
    
	return filelist;
}

/**
** calc mean image
**/
pimage_type calculate_mean_image(filelist_type list){
	fprintf(stderr,"\ncalculating mean\n");
	pimage_type image_back;
	gray* imagemap;

	pimage_type image_mean=(pimage_type)malloc(sizeof(image_type));

	for(int i=0;i<=list.size;i++){ 
		fprintf(stderr,"\n%s\n", list.paths[i]);
    	image_back=readimage(list.paths[i]);		
		fprintf(stderr,"leu\n");

		if(i==0){		
			//dimensions
		    image_mean->cols = image_back->cols;
		    image_mean->rows = image_back->rows;
		    image_mean->maxval = image_back->maxval;

		    //Allocation memoire
		    imagemap = (gray *) malloc(DPC * image_back->cols * image_back->rows * sizeof(gray));
		    image_mean->stream=imagemap;

		    //zerar image
		    for(int i=0; i < image_mean->rows; i++){
      			for(int j=0; j < image_mean->cols ; j++){
					pixel_type pixel;
					pixel.r=0;
					pixel.g=0;
					pixel.b=0;
					set_pixel(image_mean,j,i,pixel);
	   			}
     		}
		    fprintf(stderr,"zerou\n");
		}

		//Acumulation the value of the pixels of all background images
    	for(int i=0; i < image_mean->rows; i++){
      		for(int j=0; j < image_mean->cols ; j++){
				pixel_type *pixel_back=get_pixel(image_back,j,i);
				pixel_type *pixel_mean=get_pixel(image_mean,j,i);
				pixel_type newpixel;				
				newpixel.r = pixel_mean->r + pixel_back->r;
				newpixel.g = pixel_mean->g + pixel_back->g;
				newpixel.b = pixel_mean->b + pixel_back->b;
				set_pixel(image_mean,j,i,newpixel);
	   		}
     	}
    	fprintf(stderr,"acumulou\n");
    }

    //calculating the mean for each pixel
   	for(int i=0; i < image_mean->rows; i++){
   		for(int j=0; j < image_mean->cols ; j++){
			pixel_type *pixel_mean=get_pixel(image_mean,j,i);
			pixel_type newpixel;				
			newpixel.r = (int) pixel_mean->r / (list.size+1); //+1 because the first image in the folder is 000000
			newpixel.g = (int) pixel_mean->g / (list.size+1);
			newpixel.b = (int) pixel_mean->b / (list.size+1);
			set_pixel(image_mean,j,i,newpixel);
   		}
   	}
   	fprintf(stderr,"\ncalculou\n");

    printimage(image_mean);
    fprintf(stderr,"acabou\n");
}

/**
** calc the inverse of the matrix
**/
gray *matrix_inverse(gray *matrix,int rows,int cols){
	gray *inverse=(gray *)malloc(sizeof(gray)*rows*cols);
	//may generate fail duo to another dimension of matrices
	gray coefficient=1/matrix_determinant(matrix,rows);
	
	//second parameter should be the cofactor matrix	
	gray *res=matrix_multiplication_single(coefficient, matrix,rows,cols);

	return res;
}

/**
** calc the inverse of the matrix
**/
void matrix_inverse_test(){
	const int dim=3;
	gray *det=(gray *)malloc(sizeof(gray)*dim*dim);

	gray value=2;

	det[dim*0+0]=1;
	det[dim*0+1]=2;
	det[dim*0+2]=1;

	det[dim*1+0]=3;
	det[dim*1+1]=4;
	det[dim*1+2]=5;

	det[dim*2+0]=5;
	det[dim*2+1]=6;
	det[dim*2+2]=7;

	gray *res=matrix_inverse(det,dim,dim);
	matrix_print(res,dim,dim);
}

/**
** calc multiplication of a matrix by a single number test
**/
void matrix_multiplication_single_test(){

	const int dim=3;
	gray *det=(gray *)malloc(sizeof(gray)*dim*dim);

	gray value=2;

	det[dim*0+0]=1;
	det[dim*0+1]=2;
	det[dim*0+2]=1;

	det[dim*1+0]=3;
	det[dim*1+1]=4;
	det[dim*1+2]=5;

	det[dim*2+0]=5;
	det[dim*2+1]=6;
	det[dim*2+2]=7;

	gray *res=matrix_multiplication_single(value, det,dim,dim);
	matrix_print(res,3,3);

}

/**
** calc multiplication of a matrix by a single number
**/
gray *matrix_multiplication_single(gray value, gray *matrix,int rows,int cols){

	gray *result_matrix=(gray *)malloc(sizeof(gray)*rows*cols);

	for(int row=0;row<rows;row++){
		for(int col=0;col<cols;col++){
			result_matrix[row*cols+col]=value*get_matrix_pixel(matrix,row,col,cols);
		}
	}

	return result_matrix;
}

/**
** calc subtraction
**/
gray *matrix_subtraction(gray *matrix1,gray *matrix2,int rows,int cols){

	gray *result_matrix=(gray *)malloc(sizeof(gray)*rows*cols);

	for(int row=0;row<rows;row++){
		for(int col=0;col<cols;col++){
			result_matrix[row*cols+col]=get_matrix_pixel(matrix1,row,col,cols)-get_matrix_pixel(matrix2,row,col,cols);
		}
	}

	return result_matrix;
}

/**
** calc subtraction test
**/
void matrix_subtraction_test(){

	const int dim=3;

	gray *m1=(gray *)malloc(sizeof(gray)*dim*dim);
	gray *m2=(gray *)malloc(sizeof(gray)*dim*dim);

	m1[dim*0+0]=1;
	m1[dim*0+1]=2;
	m1[dim*0+2]=1;

	m1[dim*1+0]=3;
	m1[dim*1+1]=4;
	m1[dim*1+2]=5;

	m1[dim*2+0]=5;
	m1[dim*2+1]=6;
	m1[dim*2+2]=7;

	m2[dim*0+0]=1;
	m2[dim*0+1]=2;
	m2[dim*0+2]=1;

	m2[dim*1+0]=3;
	m2[dim*1+1]=4;
	m2[dim*1+2]=5;

	m2[dim*2+0]=5;
	m2[dim*2+1]=6;
	m2[dim*2+2]=7;

	gray *res=matrix_subtraction(m1,m2,dim,dim);

	matrix_print(res,3,3);

}

/**
** prints the matrix dimension
*/
void matrix_print(gray* matrix,int rows, int cols){
	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){		
			fprintf(stderr," %.3f ",matrix[i*cols+j]);	
		}
		fprintf(stderr,"\n");
	}
}

/**
** get pixel from a matrix
**/
gray get_matrix_pixel(gray *matrix,int row, int col,int dim){

	return matrix[row*dim+col];

}

/**
** calc determinant
**/
gray matrix_determinant(gray *matrix,int dim){

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
	gray *det=(gray *)malloc(sizeof(gray)*dim*dim);

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

	filelist_type list_back = readbackgrounds("background_substraction/background/img_", 2); //just work until 8

	calculate_mean_image(list_back);

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


	//jander test 
	/*
	fprintf(stderr,"Test 1\n");
	matrix_determinant_test();
	fprintf(stderr,"Test 2\n");
	matrix_subtraction_test();
	fprintf(stderr,"Test 3\n");
	matrix_multiplication_single_test(); 
	fprintf(stderr,"Test 4\n");
	matrix_inverse_test();
	// */

}
