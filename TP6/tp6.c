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
** initiate the matrix used in the mean method
**/
void initiateImageMean(pimage_type image_mean, pimage_type image_back){
	gray* imagemap;

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

/**
** calc mean image
**/
pimage_type calculate_mean_image(filelist_type list){
	fprintf(stderr,"\ncalculating mean\n");
	pimage_type image_back;

	pimage_type image_mean=(pimage_type)malloc(sizeof(image_type));

	for(int i=0;i<=list.size;i++){ 
		fprintf(stderr,"\n%s\n", list.paths[i]);
    	image_back=readimage(list.paths[i]);		
		fprintf(stderr,"leu\n");

		//printimage(image_back);

		if(i==0)
			initiateImageMean(image_mean, image_back);

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
				free(pixel_back);
				free(pixel_mean);
	   		}
     	}
    	fprintf(stderr,"acumulou\n");

    	free(image_back);
    }

    //calculating the mean for each pixel
   	for(int i=0; i < image_mean->rows; i++){
   		for(int j=0; j < image_mean->cols ; j++){
			pixel_type *pixel_mean=get_pixel(image_mean,j,i);
			pixel_type newpixel;				
			newpixel.r =  pixel_mean->r / (list.size+1); //+1 because the first image in the folder is 000000
			newpixel.g =  pixel_mean->g / (list.size+1);
			newpixel.b =  pixel_mean->b / (list.size+1);
			set_pixel(image_mean,j,i,newpixel);
			free(pixel_mean);
   		}
   	}
   	fprintf(stderr,"\ncalculou\n");

    printimage(image_mean);
    fprintf(stderr,"acabou\n");
}

/**
** test mean image method
**/
pimage_type mean_image_test(){
	filelist_type list_back = readbackgrounds("background_substraction/img_", 2);

	calculate_mean_image(list_back);
}

/**
** calc the inverse of the matrix
**/
gray *matrix_inverse(gray *matrix,int rows,int cols){
	gray *inverse=(gray *)malloc(sizeof(gray)*rows*cols);
	//may generate fail duo to another dimension of matrices
	gray coefficient=1/matrix_determinant(matrix,rows);
	
	//second parameter should be the cofactor matrix	
	gray *res=matrix_multiplication_single(coefficient, matrix_transpose(matrix_cofactor(matrix,rows,cols),rows,cols),rows,cols);

	return res;
}

/**
** calc the cofactor of the matrix
**/
gray *matrix_cofactor(gray *matrix,int rows,int cols){

	gray *result_matrix=(gray *)malloc(sizeof(gray)*rows*cols);

	gray a11=get_matrix_pixel(matrix,0,0,cols);
	gray a12=get_matrix_pixel(matrix,0,1,cols);
	gray a13=get_matrix_pixel(matrix,0,2,cols);

	gray a21=get_matrix_pixel(matrix,1,0,cols);
	gray a22=get_matrix_pixel(matrix,1,1,cols);
	gray a23=get_matrix_pixel(matrix,1,2,cols);

	gray a31=get_matrix_pixel(matrix,2,0,cols);
	gray a32=get_matrix_pixel(matrix,2,1,cols);
	gray a33=get_matrix_pixel(matrix,2,2,cols);

	gray m11=a22*a33 - a32*a23;
	gray m12=-1*(a21*a33 - a31*a23);
	gray m13=a21*a32 - a31*a22;

	gray m21=-1*(a12*a33 - a32*a13);
	gray m22=a11*a33 - a31*a13;
	gray m23=-1*(a11*a32 - a31*a12);

	gray m31=a12*a23 - a22*a13;
	gray m32=-1*(a11*a23 - a21*a13);
	gray m33=a11*a22 - a21*a12;

	result_matrix[0*cols+0]=m11;
	result_matrix[0*cols+1]=m12;
	result_matrix[0*cols+2]=m13;

	result_matrix[1*cols+0]=m21;
	result_matrix[1*cols+1]=m22;
	result_matrix[1*cols+2]=m23;

	result_matrix[2*cols+0]=m31;
	result_matrix[2*cols+1]=m32;
	result_matrix[2*cols+2]=m33;

	return result_matrix;

}

/**
** calc multiplication transpose
**/
gray *matrix_transpose(gray *matrix,int rows,int cols){

	gray *result_matrix=(gray *)malloc(sizeof(gray)*rows*cols);

	for(int row=0;row<rows;row++){
		for(int col=0;col<cols;col++){
			result_matrix[row*cols+col]=get_matrix_pixel(matrix,col,row,rows);
		}
	}

/*
	gray *result_matrix=(gray *)malloc(sizeof(gray)*rows*cols);

	gray a=get_matrix_pixel(matrix,0,0,cols);
	gray b=get_matrix_pixel(matrix,0,1,cols);
	gray c=get_matrix_pixel(matrix,0,2,cols);

	gray d=get_matrix_pixel(matrix,1,0,cols);
	gray e=get_matrix_pixel(matrix,1,1,cols);
	gray f=get_matrix_pixel(matrix,1,2,cols);

	gray g=get_matrix_pixel(matrix,2,0,cols);
	gray h=get_matrix_pixel(matrix,2,1,cols);
	gray i=get_matrix_pixel(matrix,2,2,cols);

	result_matrix[0*cols+0]=a;
	result_matrix[0*cols+1]=d;
	result_matrix[0*cols+2]=g;

	result_matrix[1*cols+0]=b;
	result_matrix[1*cols+1]=e;
	result_matrix[1*cols+2]=h;

	result_matrix[2*cols+0]=c;
	result_matrix[2*cols+1]=f;
	result_matrix[2*cols+2]=i;
*/
	return result_matrix;
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
** calc matrix multiplication
**/
gray *matrix_multiplication(gray *matrix1,int rows1,int cols1,gray *matrix2,int rows2,int cols2){

	gray *result_matrix=(gray *)malloc(sizeof(gray)*rows1*cols2);

	if(cols1!=rows2){
		fprintf(stderr,"Impossible to multiplicate those matrices\n");
		exit(0);
	}

	for(int row=0;row<rows1;row++){
		for(int col=0;col<cols2;col++){

			result_matrix[row*cols2+col]=0;

			for(int mcol=0;mcol<cols1;mcol++){
				//printf("--> %f * %f = %f \n",get_matrix_pixel(matrix1,row,mcol,cols1),get_matrix_pixel(matrix2,mcol,col,cols2),get_matrix_pixel(matrix1,row,mcol,cols1)*get_matrix_pixel(matrix2,mcol,col,cols2));
				result_matrix[row*cols2+col]=result_matrix[row*cols2+col]+get_matrix_pixel(matrix1,row,mcol,cols1)*get_matrix_pixel(matrix2,mcol,col,cols2);
			}
			//exit(0);
			//printf("--> %f\n",result_matrix[row*cols2+col]);
		}
	}

	return result_matrix;
	

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
** matrix multiplication test
**/
void matrix_multiplication_test(){


	const int dim1=3;
	const int dim2=2;

	gray *m1=(gray *)malloc(sizeof(gray)*4*3);
	gray *m2=(gray *)malloc(sizeof(gray)*3*2);

	m1[dim1*0+0]=14;
	m1[dim1*0+1]=9;
	m1[dim1*0+2]=3;

	m1[dim1*1+0]=2;
	m1[dim1*1+1]=11;
	m1[dim1*1+2]=15;

	m1[dim1*2+0]=0;
	m1[dim1*2+1]=12;
	m1[dim1*2+2]=17;

	m1[dim1*3+0]=5;
	m1[dim1*3+1]=2;
	m1[dim1*3+2]=3;

	m2[dim2*0+0]=12;
	m2[dim2*0+1]=25;

	m2[dim2*1+0]=9;
	m2[dim2*1+1]=10;

	m2[dim2*2+0]=8;
	m2[dim2*2+1]=5;

	/*m1[dim*0+0]=1;
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
	m2[dim*2+2]=7;*/

	printf("Matrix 1:\n");
	matrix_print(m1,4,3);
	printf("Matrix 2:\n");
	matrix_print(m2,3,2);

	gray *res=matrix_multiplication(m1,4,3,m2,3,2);

	printf("Result:\n");
	matrix_print(res,4,2);

}

/**
** calc determinant test
**/
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

	printf("Determinant -->%i\n",val);
}

/**
** calc covariance matrix
**/
void matrix_covariance(){

}

/**
** test covariance matrix
**/
void matrix_covariance_test(){

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

	filelist_type list_back = readbackgrounds("background_substraction/background/img_", 10); 
	calculate_mean_image(list_back);

	//mean_image_test();

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
	fprintf(stderr,"Test 5\n");
	matrix_multiplication_test();
	// */

}
