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
** gradient_x^2 
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
** gradient_y^2 
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

    // Lecture du Magic number 
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

    // Lecture des dimensions 
    lcols = pm_getint( ifp );
    lrows = pm_getint( ifp );
    lmaxval = pm_getint( ifp );

    // Allocation memoire  
    imagemap = (double *) malloc(lcols * lrows * sizeof(double));

    
    if(pgmraw){
	    // Convertendo para Plaintext 
	    // Lecture 
	    for(i=0; i < lrows; i++)
	      for(j=0; j < lcols ; j++)
		imagemap[i * lcols + j] = pm_getint(ifp);

    }else { 
        printf("Invalid file type");
        exit(1);
    }
      // fermeture 
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
	double * result;
	result = (double *)malloc(sizeof(double)*dim*dim);

	double diag1=matrix[dim*0+0]*matrix[dim*1+1];
	double diag2=(-1*matrix[dim*1+0]*matrix[dim*0+1]);

	double det=diag1+diag2;

	result[dim*0+0]=(1/det)*matrix[dim*1+1];
	result[dim*0+1]=(1/det)*-1*matrix[dim*0+1];
	result[dim*1+0]=(1/det)*-1*matrix[dim*1+0];
	result[dim*1+1]=(1/det)*matrix[dim*0+0];

	return result;
}

/**
** prints the matrix different dimension
*/
void matrixprint(double* matrix,int cols, int rows){
	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){		
			fprintf(stderr," %f ",matrix[i*cols+j]);	
		}
		fprintf(stderr,"\n");
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

	matrixprint(kernel,2,2);	
	printf("----\n");
	matrixprint(matrixinverse(kernel,2),2,2);	
}

/**
** performs the crop 
*/
double *imageextract(double *image,int lcols,int lrows,int x, int y, int xn,int ym){
	double* newimage=(double*)malloc(sizeof(double)*xn*ym);
	int new_x=0;
	int new_y=0;
	for(int row=y;row<lrows;row++){
		for(int col=x,new_x=0;col<lcols;col++){
			//make sure that the dimension is right
			//if ( ((ym+y)<=row) || ((xn+x)<=col) ) break;
			//if(((ym<=row))||((xn<=col))) break;
			//newimage[(row-y)*xn+col-x]=image[row*lcols+col];	
			if (new_x>=xn || new_y>=ym) break;
			newimage[new_y*xn+new_x]=image[row*lcols+col];	
			new_x++;
		}
		new_y++;
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
	//fprintf(stderr,"delta_i:%f delta_e_i:%i delta_f_i:%f delta_j:%f delta_e_j:%i delta_f_j:%f,\n",delta_i,delta_e_i,delta_f_i,delta_j,delta_e_j,delta_f_j);

	//initializing the array with zeros
	initialize_array(image_w,rows*cols,0);

	//interpolating
	for(int i=0;i<rows;i++)
		for(int j=0;j<cols;j++){
			//point 1
			if ( ( ((i+delta_e_i) * cols) + (j+delta_e_j)) < (rows*cols))
				image_w[((i+delta_e_i) * cols) + (j+delta_e_j)] += ( (1.0-delta_f_i) * (1.0-delta_f_j) * (image[i*cols+j]) );

			//point 2
			if ( ( ((i+delta_e_i) * cols) + (j+delta_e_j+1)) < (rows*cols))
				image_w[((i+delta_e_i) * cols) + (j+delta_e_j+1)] += ( (delta_f_i) * (1.0-delta_f_j) * (image[i*cols+j]) );

			//point 3
			if ( ( ((i+delta_e_i+1) * cols) + (j+delta_e_j)) < (rows*cols))
				image_w[((i+delta_e_i+1) * cols) + (j+delta_e_j)] += ( (1.0-delta_f_i) * (delta_f_j) * (image[i*cols+j]) );

			//point 4
			if ( ( ((i+delta_e_i+1) * cols) + (j+delta_e_j+1)) < (rows*cols))
				image_w[((i+delta_e_i+1) * cols) + (j+delta_e_j+1)] += ( (delta_f_i) * (delta_f_j) * (image[i*cols+j]) );
		}
	return image_w;	
}


/**
** creates first term of the harris 
*/
double *harrispart(double *image,int rows, int cols){
	int nr=1;
	int bin=3;
	double *image_grad_x2=callgradx(nr, bin, image, lrows, lcols);	
	double *image_grad_y2=callgrady(nr, bin, image, lrows, lcols);
	double *image_grad_xy=calc_xy(image_grad_x2, image_grad_y2, lrows, lcols);

	image_grad_x2=calc_xx(image_grad_x2, lrows, lcols);

	image_grad_y2=calc_yy(image_grad_y2, lrows, lcols);

	double a=0;
	double b=0;
	double c=0;

	for(int y=0;y<lrows;y++){
		for(int x=0;x<lcols;x++){
			int pos;
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
	
	free(image_grad_x2);
	free(image_grad_y2);
	free(image_grad_xy);

	return harris;
}


/**
** multiply matrices 
*/
double * mult_pixels_matrices(double *mat1, double *mat2, int rows, int cols){
	double * m = (double *)malloc(sizeof(double)*rows*cols);

	for(int y=0;y<rows;y++)
		for(int x=0;x<cols;x++)	
			m[y*cols+x] = mat1[y*cols+x] * mat2[y*cols+x];

	return m;
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
double sum_matrices_values(double *matrix, int rows, int cols){
	double total=0;
        for(int y=0;y<rows;y++)
	  for(int x=0;x<cols;x++)	
		total+= matrix[y*cols+x];

	return total;
}

double *findDeltas(double* image, int cols, int rows, pixel_type start, pixel_type end ,int iteration, int should_print){
	double *image1;
	double *t;

	//image1=readimage("image/tazplain/taz001.pgm");		
	image1=imageextract(image,lcols,lrows, start.x, start.y, end.x,end.y);

	t=readimage("image/tazplain/taz.pgm");

	//initialize
	double deltai=0,deltaj=0;
	double* tw=readimage("image/tazplain/taz.pgm");
	double *displacement;
	
	for(int x=0;x<iteration;x++){	
		fprintf(stderr,"1\n");	
		displacement=calc_displacements(tw,image1,lrows,lcols);
		fprintf(stderr,"2\n");	

		deltai-=displacement[0];
		deltaj-=displacement[1];
		fprintf(stderr," Iteration %i delta_i:%f delta_j:%f\n",x+1,deltai,deltaj);
		
		fprintf(stderr,"3\n");		
		tw=interpolate_image(t, lrows, lcols, deltai, deltaj);
		fprintf(stderr,"4\n");	
	}
	
	if(should_print)
		printimage(tw,lcols,lrows,lmaxval);
		
	displacement[0]=deltai;
	displacement[1]=deltaj;

	free(t);
	free(tw);
	//free(displacement);
	return displacement;
	free(image1);
}

/**
** run the ex3
*/
void exercise3(int iteration){
	double *image1,*t;

	image1=readimage("image/tazplain/taz001.pgm");	

	pixel_type start;
	start.x=50;
	start.y=130;

	pixel_type end;
	end.x=84;
	end.y=100;

	findDeltas(image1, lcols, lrows, start, end, iteration, 1);

	
}


/**
*** Save one image to the disk
**/
void saveimage(char *path,double* image,int cols,int rows,int maxval){
	FILE *fp;

	if((fp=fopen(path, "wb"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}

    fprintf(fp, "P2\n");
    fprintf(fp, "%d %d \n", cols, rows);
    fprintf(fp, "%d\n",maxval);

	for(int i=0; i < rows; i++)
    	for(int j=0; j < cols ; j++)
			fprintf(fp, "%i ",(int) image[i * cols + j] );

	fclose(fp);
}

/**
*** Create the boundbox in an image
**/
double *boundbox(double *image, int cols, int rows, int x1, int y1, int x2,int y2){
	int colorvalue=180;

	if(x1<0 || y1<0 || x2<0 || y2<0 ){
		fprintf(stderr,"function boundbox:Invalid boundary, image(%i,%i) bound(%i,%i,%i,%i)",cols,rows,x1,y1,x2,y2);
		exit(1);
	}
	if(x1 > cols || y1 > rows || x2 > cols || y2 > rows){
		fprintf(stderr,"function boundbox:Invalid boundary, image(%i,%i) bound(%i,%i,%i,%i)",cols,rows,x1,y1,x2,y2);
		exit(1);
	}

	for(int y=y1;y<rows;y++){
		for(int x=x1;x<cols;x++){
			int pos=y*cols+x;
			//left
			if(y1==y && x<x2) image[pos]=colorvalue;
			//right
			if(x==x2 && y>=y1 && y<=y2 ) image[pos]=colorvalue;
			//upper
			if(y>y1 && x==x1 && y<y2) image[pos]=colorvalue;				       //lower
			if(y==y2 && x<x2) image[pos]=colorvalue;	
		}	
	}

	return image;
}


/**
*** Test boundbox
**/
void boundboxtest(){

	double *image=readimage("image/tazplain/taz001.pgm");
	
	boundbox(image,lcols,lrows,100,100,150,120);

	printimage(image,lcols,lrows,lmaxval);
}

/**
** estimating the displacement alog i and j (equation 1 of the TP5)
** the input are the images Io(i,j) and T(i,j)
** returns a matrix 2x1
*/
double * calc_displacements(double * T, double * Io, int rows, int cols){
	//part1 of the equation
	fprintf(stderr,"antes harris\n");
	double * mat1=harrispart(Io, rows, cols);	
	fprintf(stderr,"depois harris\n");
	mat1 = matrixinverse(mat1,2);

	//part2  of the equation
	double * mat2=(double *)malloc(sizeof(double)*2*1);
	double * mat_diff=imagediff(T,Io,cols,rows);
		//(0,0)
	double * gradientI = callgrady(1, 3, Io, rows, cols); //grandient y = gradient i
	double * mat_temp = mult_pixels_matrices (gradientI, mat_diff, rows, cols); 
	double value = sum_matrices_values(mat_temp, rows,cols);
	mat2[0]= value;

	free(gradientI);
	free(mat_temp);

		//(1,0)
	gradientI = callgradx(1, 3, Io, rows, cols); //grandient x = gradient j
	mat_temp = mult_pixels_matrices (gradientI, mat_diff, rows, cols); 

	value = sum_matrices_values(mat_temp,rows,cols);
	mat2[1]=value;

	//part1 * part2
	double * displ=mult_matrices(mat1, mat2);

	free(gradientI);
	free(mat_temp);
	free(mat1);
	free(mat2);
	free(mat_diff);

	return displ;	
}

/**
** MAIN 
*/
int main(int argc, char* argv[]){
	exercise3(10); //<param> is the number of iterations
	//boundboxtest();

	//double *image1;
	//image1=readimage("taz.pgm");		
	//saveimage("out.pgm", image1, lcols, lrows, lmaxval);
	//free(image1);
}
