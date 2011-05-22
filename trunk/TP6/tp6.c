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
    for(int i=0; i < image_mean->rows; i++)
		for(int j=0; j < image_mean->cols ; j++){
			pixel_type pixel;
			pixel.r=0;
			pixel.g=0;
			pixel.b=0;
			set_pixel(image_mean,j,i,pixel);
		}
	
    //fprintf(stderr,"zerou\n");
}

/**
** calc mean image
**/
pimage_type calculate_mean_image(filelist_type list){
	fprintf(stderr,"\n****** CALCULATING MEAN *****\n");
	pimage_type image_back;

	pimage_type image_mean=(pimage_type)malloc(sizeof(image_type));

	for(int i=0;i<=list.size;i++){ 
		//fprintf(stderr,"\n%s\n", list.paths[i]);
    	image_back=readimage(list.paths[i]);		
		//fprintf(stderr,"leu\n");

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
    	//fprintf(stderr,"acumulou\n");

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
   	//fprintf(stderr,"\ncalculou\n");

	//printimage(image_mean);

    return image_mean;
}

/**
** test mean image method
**/
void calculate_mean_image_test(){
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
			//printf("aqui valor-->%i,%i; total:%i,%i\n",row,col,rows,cols);
			result_matrix[row*cols+col]=get_matrix_pixel(matrix,col,row,rows);
		}
	}

	//matrix_print(result_matrix,cols,rows);

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
** calc subtraction (matrix1-matrix2)
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
** get a pixel as a matrix
**/
gray *get_pixel_matrix(pimage_type image,int j, int i){
	gray *pixel_matrix=(gray *)malloc(sizeof(gray)*3);
	pixel_type *pixel=get_pixel(image,j,i);

	pixel_matrix[0]=pixel->r;
	pixel_matrix[1]=pixel->g;
	pixel_matrix[2]=pixel->b;
	return pixel_matrix;

}

/**
** calculate pb
**/
gray *calculate_pb(pimage_type image){

    printf("P3\n");
    printf("%d %d \n", image->cols, image->rows);
    printf("%d\n",image->maxval);

    filelist_type list_back = readbackgrounds("images/img_", 1);

    pimage_type imagemean=calculate_mean_image(list_back);

    for(int i=0; i < image->rows; i++){
      for(int j=0; j < image->cols ; j++){

		gray *pixel=get_pixel_matrix(image,j,i);

		gray *pixelmean=get_pixel_matrix(imagemean,j,i);

		gray *i_minus_mean_transpose = matrix_subtraction(pixel,pixelmean,1,3); //i - i_m

		gray *i_minus_mean = matrix_transpose(i_minus_mean_transpose,1,3); //(i - i_m)^t

		gray *sigma=matrix_covariance(imagemean,list_back,i,j);//sigma

		matrix_print(sigma,3,3);
		
		gray *sigma_inverse=matrix_inverse(sigma,3,3);//sigma^{-1}
		gray half=-1/2;//-1/2

		/*
		matrix_print(sigma_inverse,3,1);
		*/
		//parameters for exp function		

		gray *a1=matrix_multiplication(i_minus_mean_transpose,1,3,sigma_inverse,3,3);//a1=(i-i_m)^{t}*sigma^{-1}
		
		/*		
		gray *a2=matrix_multiplication(a1,1,3,i_minus_mean,3,1);//a2=a1*(i-i_m)
		gray *a3=matrix_multiplication_single(half,a2,1,1);//0.5*a2
				
		double exponent = exp(((double)*a3));

		//coefficient of exp function
		double b1 = pow(2*PI,2/3); //2pi^{3/2}
		double b2 = (double)matrix_determinant(sigma,3);//b2=|sigma|
		double b3 = b1*sqrt(b2);//b3=b1*sqrt(b2)
		double b4 = 1/(b3); //b3^{-1}

		double c1 = b4*exponent;
		*/
		printf("\n\n");
		exit(0);


      }  
    }


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
** calculate pb test
**/
void calculate_pb_test(){

	pimage_type image = readimage("images/img_000000.ppm");

	calculate_pb(image);

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

void matrix_transpose_test(){

	const int dim=3;
	gray *det=(gray *)malloc(sizeof(gray)*1*3);

	gray value=2;

	det[dim*0+0]=1;
	det[dim*0+1]=2;
	det[dim*0+2]=1;

	gray *res=matrix_transpose(det, 1,3);

	matrix_print(res,3,1);

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
** convert from the type <pixelType> to <gray>
**/
gray *convert_pixelTypeToGray(pixel_type *pixel_back){
	gray *result_matrix=(gray *)malloc(sizeof(gray)*3);

	result_matrix[0]=pixel_back->r;
	result_matrix[1]=pixel_back->g;
	result_matrix[2]=pixel_back->b;
	
	return result_matrix;
}

/**
** acumulate two matrices in the first one
**/
void acumulate_matrix(gray *covar_matrix,gray *mult_matrix, int cols, int rows){
   	for(int i=0; i < rows; i++)
   		for(int j=0; j < cols ; j++)
			covar_matrix[i*cols+j] += mult_matrix[i*cols+j]; 
}

/**
** calc covariance matrix
**/
gray * matrix_covariance(pimage_type image_mean,filelist_type list, int i, int j){
	pimage_type image_back;
	gray *in_matrix;
	gray *im_matrix;
	gray *diff_in_im_matrix;
	gray *transp_matrix;
	gray *mult_matrix;
	gray *covar_matrix=(gray *)malloc(sizeof(gray)*3*3);

	//each pixel of the image
//   	for(int i=0; i < image_mean->rows; i++){
//   		for(int j=0; j < image_mean->cols ; j++){

			//zerar covariance matrix
			for(int a=0; a < 3; a++)
				for(int b=0; b < 3 ; b++)
					covar_matrix[a*3+b]=0;

			//for each background image, calculate covariance matrix for each pixel and acumulate 
			for(int m=0;m<=list.size;m++){ 
				//fprintf(stderr,"\n%s\n", list.paths[m]);
		    	image_back=readimage(list.paths[m]);		
				//fprintf(stderr,"leu\n");

				//i_n
				pixel_type *pixel_back=get_pixel(image_back,j,i);
				in_matrix=convert_pixelTypeToGray(pixel_back);
				//fprintf(stderr,"\n  pixel: r:%3.0f g:%3.0f b:%3.0f \n",in_matrix[0],in_matrix[1],in_matrix[2]);
				//i_m
				pixel_type *pixel_mean=get_pixel(image_mean,j,i);
				im_matrix=convert_pixelTypeToGray(pixel_mean);
				//fprintf(stderr,"m pixel: r:%3.0f g:%3.0f b:%3.0f \n",im_matrix[0],im_matrix[1],im_matrix[2]);
				//subtraction
				diff_in_im_matrix=matrix_subtraction(in_matrix,im_matrix,3,1);
				//fprintf(stderr,"d pixel: r:%3.0f g:%3.0f b:%3.0f \n",diff_in_im_matrix[0],diff_in_im_matrix[1],diff_in_im_matrix[2]);
				//transpose
				transp_matrix=matrix_transpose(diff_in_im_matrix, 3, 1);
				//fprintf(stderr,"\nt pixel: r:%3.0f g:%3.0f b:%3.0f \n",transp_matrix[0],transp_matrix[1],transp_matrix[2]);
				//multiplication
				mult_matrix=matrix_multiplication(diff_in_im_matrix,3,1,transp_matrix,1,3);
				//fprintf(stderr,"covariance matrix\n");
				//matrix_print(mult_matrix,3,3);
				//acumulating in the covariance matrix
				acumulate_matrix(covar_matrix,mult_matrix,3,3);
				//fprintf(stderr,"covariance matrix acumulate\n");
				//matrix_print(covar_matrix,3,3);
				
				free(in_matrix);
				free(pixel_back);
				free(im_matrix);
				free(pixel_mean);
				free(diff_in_im_matrix);
				free(transp_matrix);
				free(mult_matrix);

		    	free(image_back);
	   		}

    		//calculating the mean of the covariance matrix
		   	for(int i=0; i < 3; i++)
		   		for(int j=0; j < 3 ; j++)
					covar_matrix[i*3+j] = covar_matrix[i*3+j] / (list.size+1); //+1 because the first image in the folder is 000000
			//fprintf(stderr,"mean of the covariance matrix\n");
			//matrix_print(covar_matrix,3,3);
//     	}
//    }


	//fprintf(stderr,"\ncalculou\n");
    return covar_matrix;
}

/**
** test covariance matrix
**/
void matrix_covariance_test(){
	fprintf(stderr,"\n****** CALCULATING COVARIANCE MATRIX *****\n");
	//filelist_type list_back = readbackgrounds("background_substraction/background/img_", 10); 
	//pimage_type image_mean=calculate_mean_image(list_back);

	filelist_type list_back = readbackgrounds("background_substraction/img_", 2);
	pimage_type image_mean=calculate_mean_image(list_back);

	//each pixel of the image
   	for(int i=0; i < image_mean->rows; i++)
   		for(int j=0; j < image_mean->cols ; j++)
			//covariance matrix
			matrix_covariance(image_mean, list_back, i, j);

	free(image_mean);	
}

/** Parser method **/

/**
** Get string param
**/
char* getStrParam(int argc,char* argv[],char* param,char* def){
	for(int i=0;i<argc;i++)
		if(strcmp(param,argv[i])==0){ 
			if(argv[i+1][0]=='-') 
				return def;
			return argv[i+1];
		}
	return def;
}

/**
** Get bool param
**/
int getBoolParam(int argc,char* argv[],char* param){
	for(int i=0;i<argc;i++)
		if(strcmp(param,argv[i])==0)
			return 1;
	return 0;
}

/**
** Get int param
**/
int getIntParam(int argc,char* argv[],char* param,char* def){
    char* res=getStrParam(argc,argv,param,def);
    
    int val=0;
 
    for(int i=0;i<strlen(res);i++){
      	if(res[0]=='-') 
      		continue;
		val+=pow(10,(strlen(res)-i-1))*(res[i]-'0'); 
    }

    return res[0]=='-'?-1*val:val;
}
/** end Parser method **/

/**
** Main method
**/
int main(int argc, char* argv[]){

	//matrix_covariance_test();

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
	//*	
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
	fprintf(stderr,"Test 6\n");
	matrix_transpose_test();
	fprintf(stderr,"Test 7\n");
	calculate_pb_test();
	// */

}
