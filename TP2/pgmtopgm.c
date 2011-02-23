#include <stdlib.h>
#include <stdio.h>
#include "Util.h"
#include "math.h"


/** Code for histogram:raquel **/
int rows, cols, maxval;

void displayImageFile(gray * graymap_strectched){
	printf("P2\n");

	printf("%d %d \n", cols, rows);
	printf("%d\n",maxval);
	   
	for(int i=0; i < rows; i++)
		for(int j=0; j < cols ; j++)
			printf("%u ",graymap_strectched[i * cols + j] );
}

void printArray(char * title, int* valueshisto){
	printf("%s",title);

	for(int i=0; i < 256; i++)
		printf("value[%d]:%d\n",i,valueshisto[i]);
}

int getMin(gray* graymap){
	int val=2147483647;
	for(int i=0; i < rows; i++)
		for(int j=0; j < cols ; j++){
			if (graymap[i * cols + j] < val )
				val=graymap[i * cols + j];
		}
	return val;
}

int getMax(gray* graymap){
	int val=0;
	for(int i=0; i < rows; i++)
		for(int j=0; j < cols ; j++){
			if (graymap[i * cols + j] > val )
				val=graymap[i * cols + j];
		}
	return val;
}

int getcdfMin(int * valueshisto){
	for(int i=0; i < 256; i++)
		if (valueshisto[i] > 0 )
			return valueshisto[i];
	return 0;
}

int * generateHistogramValues(char * title, gray* graymap){
	//printf("\n\ndimensions: rows:%d columns:%d\n",rows,cols);

	int * valueshisto=(int *) malloc(256*sizeof(int));
	int value;

	//initializing the array
	for(int i=0; i < 256; i++)
		valueshisto[i] = 0;

	//storing the number of pixels of each intensity
	for(int i=0; i < rows; i++)
		for(int j=0; j < cols ; j++){
			value=graymap[i * cols + j];
			valueshisto[value] = valueshisto[value]++;
		}

	//printing the values of the Histogram
	//printArray(title, valueshisto);

	return valueshisto;
}

void stretchingHistogram(char * title,double newRangeMin, double newRangeMax, gray* graymap){
	double actualRangeMin=(double) getMin(graymap);
	double actualRangeMax=(double) getMax(graymap);
	//printf("\nactualmax:%f actualmin:%f newRangeMax:%f newRangeMin:%f\n",actualRangeMax,actualRangeMin,newRangeMax,newRangeMin);

	double constant = ( (newRangeMax-newRangeMin) / (actualRangeMax-actualRangeMin) );

	gray* graymap_strechted;
	graymap_strechted = (gray *) malloc(cols * rows * sizeof(gray));

	//performing the stretching
	//formula: y[n] = newRangeMin + (newRangeMax-newRangeMin / actualRangeMax-actualRangeMin) * (x[n]-actuanRangeMin)
	for(int i=0; i < rows; i++)
		for(int j=0; j < cols ; j++){
			graymap_strechted[i * cols + j] = newRangeMin + (constant*(graymap[i * cols + j]-actualRangeMin) );
		}

	//generate histogram of the new image
	generateHistogramValues(title, graymap_strechted);
	
	//write the image on the screen
	displayImageFile(graymap_strechted);

	free(graymap_strechted);
}

void equalizingHistogram(char * title, gray* graymap, int* valueshisto){
	gray* graymap_equalized;
	graymap_equalized = (gray *) malloc(cols * rows * sizeof(gray));

	//acumulating the histogram
	for(int i=1;i<256;i++)
		valueshisto[i]=valueshisto[i]+valueshisto[i-1];
	//printArray("\n\nCumulative distribution function\n\n",valueshisto);


	double cdfMin=(double) getcdfMin(valueshisto); //cdf=cumulative distribution function
	double constant = ( 255.00 / ((rows*cols)-cdfMin) );

	//performing the equalization
	//formula: h(v) = round ( (cfd(v)-cdfmin / (MxN)-cdfmin) x (L-1) )
	for(int i=0; i < rows; i++)
		for(int j=0; j < cols ; j++){
			graymap_equalized[i * cols + j] = round( (valueshisto[graymap[i * cols + j]] - cdfMin) * constant);
		}


	//generate histogram of the new image
	generateHistogramValues(title, graymap_equalized);


	//write the image on the screen
	displayImageFile(graymap_equalized);
	
	free(graymap_equalized);
}

gray * readimage_int(int argc, char* argv[]){
    FILE* ifp;
    gray* imagemap;
    int ich1, ich2;
    int pgmraw;
    int j,i ;

    /* Test des arguments */
    if ( argc != 2 ){
      printf("\nUsage : %s fichier \n\n", argv[0]);
      exit(0);
    }

    /* Ouverture */
    ifp = fopen(argv[1],"r");
    if (ifp == NULL) {
      printf("erreur d'ouverture du fichier %s\n", argv[1]);
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
    else
      if(ich2 == '2')
	pgmraw = 1;
      else pgmraw = 0;

    /* Lecture des dimensions */
    cols = pm_getint( ifp );
    rows = pm_getint( ifp );
    maxval = pm_getint( ifp );

    /* Allocation memoire  */
    imagemap = (gray *) malloc(cols * rows * sizeof(gray));

    
    if(pgmraw){
    /* Convertendo para Plaintext */
    /* Lecture */
	    for(i=0; i < rows; i++)
	      for(j=0; j < cols ; j++)
		imagemap[i * cols + j] = pm_getint(ifp);

    }else { 
        printf("Invalid file type");
        exit(1);
    }
      /* fermeture */
      fclose(ifp);
      return imagemap;
}

/** Code for histogram:end **/

/** Code for filtering:jander **/

int lrows, lcols, lmaxval;

void printkernel(int dim, double *kernel){
 for(int i=0;i<3;i++){
  	for(int j=0;j<3;j++){
		printf("%f ",kernel[dim*i+j]);
	}
	printf("\n");
  }
}

void setPixel(double *image,int width, int height, int x,int y, int value){

  if(x>width || y>height){
  	printf("Setpixel: Out of range x:%i y:%i",x,y);
	exit(1);
  }

  image[width*y+x]=value;

}

int getPixel(char *image,int width, int height, int x,int y){

  if(x>width || y>height){
  	printf("Getpixel: Out of range x:%i y:%i",x,y);
	exit(1);
  }

  return image[width*y+x];

}

int findMax(double* array, int len)
{
    int max = array[0];
    for (int i=1;i<len; ++i)
        if(array[i]>max)
            max = array[i];
    return max;
}
int findMin(double* array, int len)
{
    int min = array[0];
    for (int i=1;i<len; ++i)
        if (array[i]<min)
            min = array[i];
    return min;
}
void minMax(double* oldArr, int oldMin, int oldMax, int newMin, int newMax, int len)
{
    double tmp;
    for (int i=0;i<len; ++i)
    {
        tmp = (((oldArr[i]-oldMin)/(double)(oldMax-oldMin))*(newMax-newMin))+newMin;
        oldArr[i] = (int)tmp;
    }
}

void printimage(double* image,int cols, int rows, int maxval){

printf("P2\n");
    printf("%d %d \n", cols, rows);
    printf("%d\n",maxval);

    for(int i=0; i < rows; i++)
      for(int j=0; j < cols ; j++)
	printf("%i ",(int)image[i * cols + j] );

}

double* readimage(int argc, char* argv[])
    {
    FILE* ifp;
    double* imagemap;
    int ich1, ich2;
    int pgmraw;
    int j,i ;

    /* Test des arguments */
    if ( argc != 2 ){
      printf("\nUsage : %s fichier \n\n", argv[0]);
      exit(0);
    }

    /* Ouverture */
    ifp = fopen(argv[1],"r");
    if (ifp == NULL) {
      printf("erreur d'ouverture du fichier %s\n", argv[1]);
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
    else
      if(ich2 == '2')
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

double* ApplyConvolution(int dim, double* kernel, double* image, int imageH, int imageW){
        //minMaxDouble(kernel,-1,1,0,1,9);
    /*
        for (int i=0;i<dim;++i)
            for (int j=0;j<dim;++j)
                qDebug("%f KERNEL", (double)(kernel[i*dim+j]));
     */
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

double* binomialfilter()
{
	double *kernel = (double*)malloc(sizeof(double)*9);

	kernel[0*3+0]=1;
	kernel[0*3+1]=2;
	kernel[0*3+2]=1;

	kernel[1*3+0]=2;
	kernel[1*3+1]=4;
	kernel[1*3+2]=2;

	kernel[2*3+0]=1;
	kernel[2*3+1]=2;
	kernel[2*3+2]=1;
	return kernel;
}

/** Code for filtering:end **/

int main(int argc, char* argv[]){

/*	double* image;
	double *kernel;


	kernel=binomialfilter();
	image=readimage(argc,argv);
	image=ApplyConvolution(3, kernel, image, lrows, lcols);
	printimage(image,lcols,lrows,lmaxval);
	free(image);
*/ 
//****** RAQUEL *************//
int * valueshisto;
gray * image_int;
image_int=readimage_int(argc,argv);

//write the image on the screen
//displayImageFile(image_int);
	
//generate histogram
valueshisto=generateHistogramValues("\nHistogram's values\n\n",image_int);

//histogram stretching
//stretchingHistogram("\nStretched Histogram's values\n\n", 0,255,image_int);

//histogram equalization
equalizingHistogram("\nEqualized Histogram's values\n\n",image_int,valueshisto);

free(image_int);
return 0;

}
