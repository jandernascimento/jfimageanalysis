#include <stdlib.h>
#include <stdio.h>
#include "tp5.h"

int lcols, lrows, lmaxval;

/**
** find the max
*/
int findMax(double* array, int len){
    int max = array[0];
    for (int i=1;i<len; ++i)
        if(array[i]>max)
            max = array[i];
    return max;
}


/**
** find the minimum
*/
int findMin(double* array, int len){
    int min = array[0];
    for (int i=1;i<len; ++i)
        if (array[i]<min)
            min = array[i];
    return min;
}

/**
** equalize
*/
void minMax(double* oldArr, int oldMin, int oldMax, int newMin, int newMax, int len){
    double tmp;
    for (int i=0;i<len; ++i)
    {
        tmp = (((oldArr[i]-oldMin)/(double)(oldMax-oldMin))*(newMax-newMin))+newMin;
        oldArr[i] = (int)tmp;
    }
}

/**
** performs the convolution
*/
double* ApplyConvolution(int dim, double* kernel, double* image, int imageH, int imageW){
       int j;  // row    index of the current image
       int i;  // column index of the current image
       int jk; // row    index of the kernel;
       int ik; // column index of the kernel;
	int ii, jj;
       double newval; // new colors
       int size = imageH*imageW;
       double* tmpImage = (double*)(malloc(sizeof(double)*size));
       int kernelCenteri; // index of the central column of the kernel
       int kernelCenterj; // index of the central row of the kernel
       double kernelTotalValue;
       kernelCenteri = dim / 2;
       kernelCenterj = dim / 2;
       kernelTotalValue = 0.0;
       for (j = 0; j < dim; j++)
         for(i = 0; i < dim; i++)
           kernelTotalValue += (double)(kernel[j*dim+i]);
       if (kernelTotalValue<=0)
           kernelTotalValue=1;
       // convolution computation
       for (j = 0; j < imageH; j++) {
         for (i = 0; i < imageW; i++) {
           newval = 0;
           for (jk = 0; jk < dim; jk++) {
             for (ik = 0; ik < dim; ik++) {
               ii = i + ik - kernelCenteri;
               jj = j + jk - kernelCenterj;
               if ((jj >= 0) && (jj <imageH ) && (ii >= 0) && (ii < imageW))
               {
                   newval += image[jj*imageW+ii]*(double)(kernel[jk*dim+ik]);
               }
             }
           }
           newval = newval / kernelTotalValue;

           tmpImage[j*imageW+i] = newval;

         }
       }
       minMax(tmpImage,findMin(tmpImage,size),findMax(tmpImage,size),0,255,size);
       return tmpImage;
}

/**
** horizontal sobel filter
*/
double* sobelfilterH(int val){
	double *kernel = (double*)malloc(sizeof(double)*val*val);
	
	if(val==3){
		kernel[0*3+0]=-1;
		kernel[0*3+1]=-2;
		kernel[0*3+2]=-1;

		kernel[1*3+0]=0;
		kernel[1*3+1]=0;
		kernel[1*3+2]=0;

		kernel[2*3+0]=1;
		kernel[2*3+1]=2;
		kernel[2*3+2]=1;
		return kernel;
	}	
	exit(-1);
}

/**
** vertical sobel filter
*/
double* sobelfilterV(int val){
	double *kernel = (double*)malloc(sizeof(double)*val*val);
	
	if(val==3){
		kernel[0*3+0]=-1;
		kernel[0*3+1]=0;
		kernel[0*3+2]=1;

		kernel[1*3+0]=-2;
		kernel[1*3+1]=0;
		kernel[1*3+2]=2;

		kernel[2*3+0]=-1;
		kernel[2*3+1]=0;
		kernel[2*3+2]=1;
		return kernel;
	}	
	exit(-1);
}

/**
** apply the norm of the gradient
*/
double * normGradient(double *gx, double *gy, int lrows, int lcols){
	double * g = (double *) malloc(sizeof(double)*lrows*lcols);
	for (int x=0;x<lcols;x++)
		for (int y=0;y<lrows;y++){
			//sqrt (val1*val1 + val2*val2)
			int val1 = gx[y*lcols+x] * gx[y*lcols+x];
			int val2 = gy[y*lcols+x] * gy[y*lcols+x];
			g[y*lcols+x] = sqrt( val1 + val2);//*/

			/*/mod(val1)+mod(val2)
			int val1 = gx[y*lcols+x];
			int val2 = gy[y*lcols+x];
			g[y*lcols+x] = abs( val1 - val2);*/
		}

	return g;
}

/**
** gradient in x
*/
double *callgradx(int times, int kernelsize, double *img, int height, int width){
	double *kernel=sobelfilterV(3);
	double *resultImage;
	resultImage=ApplyConvolution(3, kernel, img, height, width);
	free(kernel);
	return resultImage;
}

/**
** gradient in y 
*/
double *callgrady(int times, int kernelsize, double *img, int height, int width){
	double *kernel=sobelfilterH(3);
	double *resultImage;
	resultImage=ApplyConvolution(3, kernel, img, height, width);
	free(kernel);
	return resultImage;
}

/**
** gradient in xy
*/
double *callgradxy(int times, int kernelsize, double *img, int height, int width){
	double * kernelV=sobelfilterV(3);
	double * kernelH=sobelfilterH(3);
	double * gx=ApplyConvolution(3, kernelV, img, height, width);
	double * gy=ApplyConvolution(3, kernelH, img, height, width);
	double *resultImage=normGradient(gx, gy, height, width);
	free(gx);
	free(gy);	
	free(kernelV);
	free(kernelH);	
	return resultImage;			
}

/**
** gradient in x^2 
*/
double * calc_xx(double *gx, int lrows, int lcols){
	double * g = (double *) malloc(sizeof(double)*lrows*lcols);
	for (int x=0;x<lcols;x++)
		for (int y=0;y<lrows;y++){
			g[y*lcols+x] = gx[y*lcols+x] * gx[y*lcols+x];
		}

	return g;
}

/**
** gradient in y^2 
*/
double * calc_yy(double *gy, int lrows, int lcols){
	double * g = (double *) malloc(sizeof(double)*lrows*lcols);
	for (int x=0;x<lcols;x++)
		for (int y=0;y<lrows;y++){
			g[y*lcols+x] = gy[y*lcols+x] * gy[y*lcols+x];
		}

	return g;
}

/**
** calculates the c term of harris 
*/
double * calc_xy(double *gx, double *gy, int lrows, int lcols){
	double * g = (double *) malloc(sizeof(double)*lrows*lcols);
	for (int x=0;x<lcols;x++)
		for (int y=0;y<lrows;y++){
			g[y*lcols+x] = gx[y*lcols+x] * gy[y*lcols+x];
		}

	return g;
}


/**
** reads the image 
*/
double* readimage(char* filepath){
    FILE* ifp;
    double* imagemap;
    int ich1, ich2;
    int pgmraw;
    int j,i ;

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
    if(ich2 != '2' && ich2 != '5')
      pm_erreur(" mauvais type de fichier ");
    else if(ich2 == '2')
		pgmraw = 1;
    else pgmraw = 0;

    /* Lecture des dimensions */
    lcols = pm_getint( ifp );
    lrows = pm_getint( ifp );
    lmaxval = pm_getint( ifp );

    /* Allocation memoire  */
    imagemap = (double *) malloc(lcols * lrows * sizeof(double));

    
    if(pgmraw){
    /* Convertendo para Plaintext */
    /* Lecture */
	    for(i=0; i < lrows; i++)
	      for(j=0; j < lcols ; j++)
		imagemap[i * lcols + j] = pm_getint(ifp);

    }else { 
        printf("Invalid file type");
        exit(1);
    }
      /* fermeture */
      fclose(ifp);
      return imagemap;
}

/**
** prints the image
*/
void printimage(double* image,int cols, int rows, int maxval){
	printf("P2\n");
    printf("%d %d \n", cols, rows);
    printf("%d\n",maxval);

    for(int i=0; i < rows; i++)
      for(int j=0; j < cols ; j++)
	printf("%i ",(int)image[i * cols + j] );
}

/**
** calculates the inverse 
*/
double* matrixinverse(double* matrix,int dim){

	fprintf(stderr,"antes maloc. dimensao:%i\n",dim);
	
	double * result;
	result = (double *)malloc(sizeof(double)*dim*dim);

	fprintf(stderr,"depois maloc");

	double diag1=matrix[dim*0+0]*matrix[dim*1+1];
	double diag2=(-1*matrix[dim*1+0]*matrix[dim*0+1]);

	double det=diag1+diag2;
	fprintf(stderr,"determinante:%f",det);

	result[dim*0+0]=(1/det)*matrix[dim*1+1];
	result[dim*0+1]=(1/det)*-1*matrix[dim*0+1];
	result[dim*1+0]=(1/det)*-1*matrix[dim*1+0];
	result[dim*1+1]=(1/det)*matrix[dim*0+0];

	return result;
}

/**
** prints the matrix
*/
void matrixprint(double* matrix,int dim){

	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){		
			printf(" %0.2f ",matrix[i*dim+j]);	
		}
		printf("\n");
	}

}

/**
** test the matrix
*/
void matrixtest(){

	double * kernel = (double *)malloc(sizeof(double)*2*2);

	kernel[0*2+0]= -91;
	kernel[0*2+1]= 200;

	kernel[1*2+0]= 3.5;
	kernel[1*2+1]= -4.7;

	matrixprint(kernel,2);	
	printf("----\n");
	matrixprint(matrixinverse(kernel,2),2);	

}

/**
** performs the crop 
*/
double *imageextract(double *image,int lcols,int lrows,int x, int y, int xn,int ym){

	double* newimage=(double*)malloc(sizeof(double)*xn*ym);

	for(int row=y;row<lrows;row++){
		for(int col=x;col<lcols;col++){
			//make sure that the dimension is right
			if(((ym<=row))||((xn<=col))) break;
			newimage[(row-y)*xn+col-x]=image[row*lcols+col];	
		}
	}

	return newimage;
}

/**
** make the diff between images 
*/
double *imagediff(double *image1,double *image2,int cols,int rows){

	double* newimage=(double*)malloc(sizeof(double)*cols*rows);

	for(int y=0;y<rows;y++){
		for(int x=0;x<cols;x++){
			newimage[y*cols+x]=image1[y*cols+x]-image2[y*cols+x];
		}
	}

	return newimage;
}


/**
** test the extract 
*/
void imageextracttest(){
	double * image=readimage("image/tazplain/taz001.pgm");
	double *resultImage=imageextract(image,lcols,lrows,0,0,200,200);
	resultImage=callgradxy(1, 3, image, lrows, lcols);
	printimage(resultImage,lrows,lcols,lmaxval);
}


/**
** initialize the array 
*/
void initialize_array(double * vec, int n, int vldef){
	for(int i=0;i<n;i++)
		vec[i]=vldef;
}

/**
** interpolate the image
*/
double * interpolate_image(double *image, int rows, int cols, double delta_i, double delta_j){
	double * image_w = (double *)malloc(sizeof(double)*rows*cols);

	int delta_e_i     = (int) delta_i;
	double delta_f_i;
	int delta_e_j     = (int) delta_j;
	double delta_f_j;

	if (delta_i < 0)	delta_f_i = (double) (-1)*delta_i + delta_e_i; 
	else				delta_f_i = (double)      delta_i - delta_e_i; 
	if (delta_j < 0)	delta_f_j = (double) (-1)*delta_j + delta_e_j; 
	else				delta_f_j = (double)      delta_j - delta_e_j; 
	fprintf(stderr,"delta_i:%f delta_e_i:%i delta_f_i:%f delta_j:%f delta_e_j:%i delta_f_j:%f,\n",delta_i,delta_e_i,delta_f_i,delta_j,delta_e_j,delta_f_j);

	//initializing the array with zeros
	initialize_array(image_w,rows*cols,0);

	//interpolating
	for(int i=0;i<rows;i++)
		for(int j=0;j<cols;j++){
			//point 1
			if ( ( ((i+delta_e_i) * rows) + (j+delta_e_j)) < (rows*cols))
				image_w[((i+delta_e_i) * rows) + (j+delta_e_j)] += ( (1.0-delta_f_i) * (1.0-delta_f_j) * (image[i*rows+j]) );

			//point 2
			if ( ( ((i+delta_e_i) * rows) + (j+delta_e_j+1)) < (rows*cols))
				image_w[((i+delta_e_i) * rows) + (j+delta_e_j+1)] += ( (delta_f_i) * (1.0-delta_f_j) * (image[i*rows+j]) );

			//point 3
			if ( ( ((i+delta_e_i+1) * rows) + (j+delta_e_j)) < (rows*cols))
				image_w[((i+delta_e_i+1) * rows) + (j+delta_e_j)] += ( (1.0-delta_f_i) * (delta_f_j) * (image[i*rows+j]) );

			//point 4
			if ( ( ((i+delta_e_i+1) * rows) + (j+delta_e_j+1)) < (rows*cols))
				image_w[((i+delta_e_i+1) * rows) + (j+delta_e_j+1)] += ( (delta_f_i) * (delta_f_j) * (image[i*rows+j]) );
		}

	return image_w;	
}


/**
** creates first term of the harris 
*/
double *harrispart(double *image,int rows, int cols){
			//image=readimage("image/tazplain/taz001.pgm");
			int nr=1;
			int bin=3;
			double *image_grad_x2=callgradx(nr, bin, image, lrows, lcols);	
			double *image_grad_y2=callgrady(nr, bin, image, lrows, lcols);
			double *image_grad_xy=calc_xy(image_grad_x2, image_grad_y2, lrows, lcols);

			image_grad_x2=calc_xx(image_grad_x2, lrows, lcols);
			//image_grad_x2=callbinomial(1, 3, image_grad_x2, lrows, lcols);	//smothing
			image_grad_y2=calc_yy(image_grad_y2, lrows, lcols);
			//image_grad_y2=callbinomial(1, 3, image_grad_y2, lrows, lcols);  //smothing
			//image_grad_xy=callbinomial(1, 3, image_grad_xy, lrows, lcols);  //smothing		

			double a=0;
			double b=0;
			double c=0;

			for(int y=0;y<lrows;y++){
				for(int x=0;x<lcols;x++){
					int pos;
					double part=0;
					pos=y*lcols+x;
					a+=image_grad_x2[pos];
					b+=image_grad_y2[pos];
					c+=image_grad_xy[pos];
				}
			}


			double *harris=(double *)malloc(sizeof(double)*4);
			harris[0*2+0]=a;
			harris[0*2+1]=c;
			harris[1*2+0]=c;
			harris[1*2+1]=b;

			//free(image_grad_x2);
			//free(image_grad_y2);
			//free(image_grad_xy);
			
			return harris;
}


/**
** multiply matrices 
*/
double * mult_pixels_matrices(double *mat1, double *mat2, int rows, int cols){
	double * result = (double *)malloc(sizeof(double)*rows*cols);

	for(int y=0;y<rows;y++)
		for(int x=0;x<cols;x++)	
			result[y*rows+x] = mat1[y*rows+x] * mat2[y*rows+x];

	return result;
}

/**
** multiply mathematically matrices 
*/
double * mult_matrices(double *mat1, double *mat2){
	double * result = (double *)malloc(sizeof(double)*2*1);

	result[0] = (mat1[0] * mat2[0]) + (mat1[1] * mat2[1]);
	result[1] = (mat1[2] * mat2[0]) + (mat1[3] * mat2[1]);

	return result;
}

/**
** sum the values of a matrix 
*/
double sum_matrices_values(double *matrix,int dim){
	double total=0;
        for(int y=0;y<dim;y++)
	  for(int x=0;x<dim;x++)	
		total+= matrix[y*dim+x];

	return total;


}

/**
** run the ex3
*/
void exercise3(int iteration){

	double *image1,*t;

	image1=readimage("image/tazplain/taz001.pgm");
		
	image1=imageextract(image1,lcols,lrows,130, 50, 100, 84);

	t=readimage("image/tazplain/taz.pgm");

	//initialize
	int deltai=0,deltaj=0;
	double* tw=t;

	double *displacement;
	for(int x=0;x<iteration;x++){
		
		displacement=calc_displacements(t,image1,lrows,lcols);

		deltai-=displacement[0];
		deltaj-=displacement[1];
		
		tw=interpolate_image(t, lrows, lcols, deltaj, deltai);

	}
	printimage(tw,lcols,lrows,lmaxval);

	//free(image1);
	//free(t);
	//free(tw);
	//free(displacement);
}

/**
** estimating the displacement alog i and j (equation 1 of the TP5)
** the input are the images Io(i,j) and T(i,j)
** returns a matrix 2x1
*/
double * calc_displacements(double * T, double * Io, int rows, int cols){
	//part1 of the equation
	double * mat1=harrispart(Io, rows, cols);	
	mat1 = matrixinverse(mat1,2);


	//part2  of the equation
	double * mat_diff=imagediff(T,Io,cols,rows);
		//(0,0)
	double * gradientI = callgrady(1, 3, Io, rows, cols); //grandient y = gradient i
	double * mat_temp = mult_pixels_matrices (gradientI, mat_diff, rows, cols); 
	double value = sum_matrices_values(mat_temp,cols);
	double * mat2=(double *)malloc(sizeof(double)*2*1);
	mat2[0]= value;
		//(1,0)
	gradientI = callgradx(1, 3, Io, rows, cols); //grandient x = gradient j
	//mat_temp = mult_pixels_matrices (gradientI, mat_diff, rows, cols); 
	//value = sum_matrices_values(mat_temp,cols);
	mat2[1]=value;


	//part1 * part2
	double * displ=(double *)malloc(sizeof(double)*2*1);//=mult_matrices(mat1, mat2);

	//free(mat1);

	return displ;	


}


/**
** MAIN 
*/
int main(int argc, char* argv[]){
	//double * image=readimage("taz_ascii.pgm");

	exercise3(1); //<param> is the number of iterations

	//free(image);
}
