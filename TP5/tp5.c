#include <stdlib.h>
#include <stdio.h>
#include "tp5.h"

int lcols, lrows, lmaxval;

/** Parser method **

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
** end Parser method **/


/*************************** GRADIENT *************************************************************/
int findMax(double* array, int len){
    int max = array[0];
    for (int i=1;i<len; ++i)
        if(array[i]>max)
            max = array[i];
    return max;
}

int findMin(double* array, int len){
    int min = array[0];
    for (int i=1;i<len; ++i)
        if (array[i]<min)
            min = array[i];
    return min;
}

void minMax(double* oldArr, int oldMin, int oldMax, int newMin, int newMax, int len){
    double tmp;
    for (int i=0;i<len; ++i)
    {
        tmp = (((oldArr[i]-oldMin)/(double)(oldMax-oldMin))*(newMax-newMin))+newMin;
        oldArr[i] = (int)tmp;
    }
}

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
			g[y*lcols+x] = abs( val1 - val2);//*/
		}

	return g;
}
double *callgradx(int times, int kernelsize, double *img, int height, int width){
	double *kernel=sobelfilterV(3);
	double *resultImage;
	resultImage=ApplyConvolution(3, kernel, img, height, width);
	free(kernel);
	return resultImage;
}

double *callgrady(int times, int kernelsize, double *img, int height, int width){
	double *kernel=sobelfilterH(3);
	double *resultImage;
	resultImage=ApplyConvolution(3, kernel, img, height, width);
	free(kernel);
	return resultImage;
}

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

double * calc_xx(double *gx, int lrows, int lcols){
	double * g = (double *) malloc(sizeof(double)*lrows*lcols);
	for (int x=0;x<lcols;x++)
		for (int y=0;y<lrows;y++){
			g[y*lcols+x] = gx[y*lcols+x] * gx[y*lcols+x];
		}

	return g;
}

double * calc_yy(double *gy, int lrows, int lcols){
	double * g = (double *) malloc(sizeof(double)*lrows*lcols);
	for (int x=0;x<lcols;x++)
		for (int y=0;y<lrows;y++){
			g[y*lcols+x] = gy[y*lcols+x] * gy[y*lcols+x];
		}

	return g;
}

double * calc_xy(double *gx, double *gy, int lrows, int lcols){
	double * g = (double *) malloc(sizeof(double)*lrows*lcols);
	for (int x=0;x<lcols;x++)
		for (int y=0;y<lrows;y++){
			g[y*lcols+x] = gx[y*lcols+x] * gy[y*lcols+x];
		}

	return g;
}

/*************************************************************************************************/

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

void printimage(double* image,int cols, int rows, int maxval){
	printf("P2\n");
    printf("%d %d \n", cols, rows);
    printf("%d\n",maxval);

    for(int i=0; i < rows; i++)
      for(int j=0; j < cols ; j++)
	printf("%i ",(int)image[i * cols + j] );
}

//Crop an image
double * crop_image(double *image, int x1, int y1, int xn, int xm){
}

double* matrixinverse(double* matrix,int dim){

	double * res = (double *)malloc(sizeof(double)*dim*dim);

	double diag1=matrix[dim*0+0]*matrix[dim*1+1];
	double diag2=(-1*matrix[dim*1+0]*matrix[dim*0+1]);

	double det=diag1+diag2;

	res[dim*0+0]=(1/det)*matrix[dim*1+1];
	res[dim*0+1]=(1/det)*-1*matrix[dim*0+1];
	res[dim*1+0]=(1/det)*-1*matrix[dim*1+0];
	res[dim*1+1]=(1/det)*matrix[dim*0+0];

	return res;
}

void matrixprint(double* matrix,int dim){

	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){		
			printf(" %0.2f ",matrix[i*dim+j]);	
		}
		printf("\n");
	}

}

void matrixtest(){

	double * kernel = (double *)malloc(sizeof(double)*2*2);

	kernel[0*2+0]=1;
	kernel[0*2+1]=2;

	kernel[1*2+0]=3;
	kernel[1*2+1]=4;

	matrixprint(kernel,2);	
	printf("----\n");
	matrixprint(matrixinverse(kernel,2),2);	

}

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

//the values here may end up negative, due to the subtraction without norm
double *imagediff(double *image1,double *image2,int cols,int rows){

	double* newimage=(double*)malloc(sizeof(double)*cols*rows);

	for(int y=0;y<rows;y++){
		for(int x=0;x<cols;x++){
			newimage[y*cols+x]=image1[y*cols+x]-image2[y*cols+x];
		}
	}
}


void imageextracttest(){
	double * image=readimage("image/tazplain/taz001.pgm");
	double *resultImage=imageextract(image,lcols,lrows,0,0,200,200);
	resultImage=callgradxy(1, 3, image, lrows, lcols);
	printimage(resultImage,lrows,lcols,lmaxval);
}


//********************************************* INTERPOLATION ********************************************/
void initialize_array(double * vec, int n, int vldef){
	for(int i=0;i<n;i++)
		vec[i]=vldef;
}

double * interpolate_image(double *image, int rows, int cols, double delta_x, double delta_y){
	double * image_w = (double *)malloc(sizeof(double)*rows*cols);

	int delta_e_x     = (int) delta_x;
	double delta_f_x;
	int delta_e_y     = (int) delta_y;
	double delta_f_y;

	if (delta_x < 0)	delta_f_x = (double) (-1)*delta_x + delta_e_x; 
	else				delta_f_x = (double)      delta_x - delta_e_x; 
	if (delta_y < 0)	delta_f_y = (double) (-1)*delta_y + delta_e_y; 
	else				delta_f_y = (double)      delta_y - delta_e_y; 
	fprintf(stderr,"delta_x:%f delta_e_x:%i delta_f_x:%f delta_y:%f delta_e_y:%i delta_f_y:%f,\n",delta_x,delta_e_x,delta_f_x,delta_y,delta_e_y,delta_f_y);

	//initializing the array with zeros
	initialize_array(image_w,rows*cols,0);

	//interpolating
	for(int y=0;y<rows;y++)
		for(int x=0;x<cols;x++){
			//point 1
			if ( ( ((y+delta_e_y) * rows) + (x+delta_e_x)) < (rows*cols))
				image_w[((y+delta_e_y) * rows) + (x+delta_e_x)] += ( (1.0-delta_f_x) * (1.0-delta_f_y) * (image[y*rows+x]) );

			//point 2
			if ( ( ((y+delta_e_y) * rows) + (x+delta_e_x+1)) < (rows*cols))
				image_w[((y+delta_e_y) * rows) + (x+delta_e_x+1)] += ( (delta_f_x) * (1.0-delta_f_y) * (image[y*rows+x]) );

			//point 3
			if ( ( ((y+delta_e_y+1) * rows) + (x+delta_e_x)) < (rows*cols))
				image_w[((y+delta_e_y+1) * rows) + (x+delta_e_x)] += ( (1.0-delta_f_x) * (delta_f_y) * (image[y*rows+x]) );

			//point 4
			if ( ( ((y+delta_e_y+1) * rows) + (x+delta_e_x+1)) < (rows*cols))
				image_w[((y+delta_e_y+1) * rows) + (x+delta_e_x+1)] += ( (delta_f_x) * (delta_f_y) * (image[y*rows+x]) );
		}

	return image_w;	
}
//********************************************************************************************************/


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


			double *result=(double *)malloc(sizeof(double)*4);
			result[0*4+0]=a;
			result[0*4+1]=c;
			result[1*4+0]=c;
			result[1*4+1]=b;
			
			return result;
}


//multiplies two matrices, taking each pixel and multiply each other
double * mult_pixels_matrices(double *mat1, double *mat2, int rows, int cols){
	double * result = (double *)malloc(sizeof(double)*rows*cols);

	for(int y=0;y<rows;y++)
		for(int x=0;x<cols;x++)	
			result[y*rows+x] = mat1[y*rows+x] * mat2[y*rows+x];

	return result;
}

int main(int argc, char* argv[]){
	double * image=readimage("taz_ascii.pgm");
	//interpolation
	double *resultImage=interpolate_image(image, lrows, lcols, 5, 5);
	printimage(resultImage,lcols,lrows,lmaxval);

	//matrixtest();
	//imageextracttest();	
	/*double * image=readimage("taz_ascii.pgm");
	double * kernel = (double *)malloc(sizeof(double)*3*3);
	
	kernel[0*3+0]=-1;
	kernel[0*3+1]=0;
	kernel[0*3+2]=1;

	kernel[1*3+0]=-2;
	kernel[1*3+1]=0;
	kernel[1*3+2]=2;

	kernel[2*3+0]=-1;
	kernel[2*3+1]=0;
	kernel[2*3+2]=1;
	matrixprint(kernel,3);	
	printf("------------\n");
	double *resultImage=ApplyConvolution(3, kernel, image, lrows, lcols);

	printimage(resultImage,lcols,lrows,lmaxval);
	*/
}
