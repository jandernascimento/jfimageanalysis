#include <stdlib.h>
#include <stdio.h>
#include "tp5.h"

int lcols, lrows, lmaxval;

pixel_type *get_pixel(pimage_type image,int x, int y){
  pixel_type *pixel=(pixel_type*)malloc(sizeof(pixel_type));
  
  //pixel->r = image->stream[DPC*y * image->cols + DPC*x+RED];
  //pixel->g = image->stream[DPC*y * image->cols + DPC*x+GREEN];
  pixel->g = image->stream[DPC*y * image->cols + DPC*x];
  //pixel->b = image->stream[DPC*y * image->cols + DPC*x+BLUE];
  //pixel->k = image->stream[DPC*y * image->cols + DPC*x+K];
  pixel->x = x;
  pixel->y = y;
  
  return pixel;
}


void set_pixel(pimage_type image,int x, int y, pixel_type pixel ){
  //image->stream[DPC*y * image->cols + DPC*x+RED]   = pixel.r;
  //image->stream[DPC*y * image->cols + DPC*x+GREEN] = pixel.g;
  image->stream[DPC*y * image->cols + DPC*x] = pixel.g;
  //image->stream[DPC*y * image->cols + DPC*x+BLUE]  = pixel.b;
}

/* prints the image in the structure format
void printimage(pimage_type image){

    printf("P5\n");
    printf("%d %d \n", image->cols, image->rows);
    printf("%d\n",image->maxval);

    for(int i=0; i < image->rows; i++){
      for(int j=0; j < image->cols ; j++){
		//printf("%u ",get_pixel(image,j,i)->r);
		//printf("%u ",get_pixel(image,j,i)->g);
		printf("%c",get_pixel(image,j,i)->g);
		//printf("%u \n",get_pixel(image,j,i)->b);
      }
      //printf("\n");
    }
}
*/

/* reads the image in the structure format
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

    // Lecture du Magic number 
    ich1 = getc( ifp );
    if ( ich1 == EOF )
        pm_erreur( "EOF / erreur de lecture / nombre magique" );
    ich2 = getc( ifp );
    if ( ich2 == EOF )
        pm_erreur( "EOF / erreur de lecture / nombre magique" );
    if(ich2 != '5'){
    	fprintf(stderr,"Invalid file type");
        exit(1);
    }

    // Lecture des dimensions 
    image->cols = pm_getint( ifp );
    image->rows = pm_getint( ifp );
    image->maxval = pm_getint( ifp );

    // Allocation memoire  
    imagemap = (gray *) malloc(DPC * image->cols * image->rows * sizeof(gray));
    image->stream=imagemap;

    // Lecture 
    for(int i=0; i < image->rows; i++){
      for(int j=0; j < image->cols ; j++){
		pixel_type pixel;
		//pixel.r=pm_getint(ifp);
		//pixel.g=pm_getint(ifp);
		pixel.g=pm_getrawbyte(ifp);
		//pixel.b=pm_getint(ifp);
		set_pixel(image,j,i,pixel);
	   }
     }

      // fermeture 
      fclose(ifp);
      return image; //imagemap;
}
*/

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
	//return kernel;

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




int main(int argc, char* argv[]){	
	double * image=readimage("taz_ascii.pgm");

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
		
	double *resultImage=ApplyConvolution(3, kernel, image, lrows, lcols);

	printimage(resultImage,lcols,lrows,lmaxval);
}
