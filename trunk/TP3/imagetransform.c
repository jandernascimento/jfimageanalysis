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

int * generateHistogramValues(char * title, gray* graymap, int print){
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
	if (print)
		printArray(title, valueshisto);

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
	generateHistogramValues(title, graymap_strechted,0);
	
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
	generateHistogramValues(title, graymap_equalized,0);


	//write the image on the screen
	displayImageFile(graymap_equalized);
	
	free(graymap_equalized);
}

gray * readimage_int(char *filepath){
    FILE* ifp;
    gray* imagemap;
    int ich1, ich2;
    int pgmraw;
    int j,i ;

    ifp = fopen(filepath,"r");
    if (ifp == NULL) {
      printf("erreur d'ouverture du fichier %s\n", filepath);
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

double* readimage(char* filepath)
    {
    FILE* ifp;
    double* imagemap;
    int ich1, ich2;
    int pgmraw;
    int j,i ;
/*
    if ( argc != 2 ){
      printf("\nUsage : %s fichier \n\n", argv[0]);
      exit(0);
    }
*/
    ifp = fopen(filepath/*argv[1]*/,"r");
    if (ifp == NULL) {
      printf("erreur d'ouverture du fichier %s\n", filepath);
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

void windowSort(double *window,int dim){
  for (int j = 0; j < (dim*dim); ++j){
    //   Find position of minimum element
    int min = j;
    for (int l = j + 1; l < (dim*dim); ++l)
      if (window[l] < window[min])
        min = l;
      //Put found minimum element in its place
      double temp = window[j];
      window[j] = window[min];
      window[min] = temp;
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

double * normGradient(double *gx, double *gy, int dim){
	double * g = (double *) malloc(sizeof(double)*dim*dim);
	for (int i=0;i<dim;i++)
		for (int j=0;j<dim;j++){
			int val1 = pow(gx[j*dim+i], 2);
			int val2 = pow(gy[j*dim+i], 2);
			g[j*dim+i] = sqrt( val1 + val2);
		}

	return g;
}

double* binomialfilter(int val)
{
	double *kernel = (double*)malloc(sizeof(double)*val*val);
	if(val==3){
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
	}else{
		kernel[0*3+0]=1;
		kernel[0*3+1]=4;
		kernel[0*3+2]=6;
		kernel[0*3+3]=4;
		kernel[0*3+4]=1;

		kernel[1*3+0]=4;
		kernel[1*3+1]=16;
		kernel[1*3+2]=24;
		kernel[1*3+3]=16;
		kernel[1*3+4]=4;

		kernel[2*3+0]=6;
		kernel[2*3+1]=24;
		kernel[2*3+2]=36;
		kernel[2*3+3]=24;
		kernel[2*3+4]=6;

		kernel[3*3+0]=4;
		kernel[3*3+1]=16;
		kernel[3*3+2]=24;
		kernel[3*3+3]=16;
		kernel[3*3+4]=4;

		kernel[4*3+0]=1;
		kernel[4*3+1]=4;
		kernel[4*3+2]=6;
		kernel[4*3+3]=4;
		kernel[4*3+4]=1;

		return kernel;
	}
	exit(-1);
}

/** Code for filtering:end **/


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

void medianFilter(double * image, double * result, int N, int M){
	//Move window through all elements of the image
	for (int m = 1; m < M - 1; ++m)
		for (int n = 1; n < N - 1; ++n){
			//Pick up window elements
			int k = 0;
			double window[9]; //9 is the size of the kernel
			for (int j = m - 1; j <= m + 1; ++j)
				for (int i = n - 1; i <= n + 1; ++i)
					window[k++] = image[j * N + i];

			//Order elements
			windowSort(window,3); //3 is one dimension of the kernel

			//Get result - the middle element
			result[m*N+n] = window[4];
		}
}

int main(int argc, char* argv[]){
	double* image;
	double *kernel,*kernel2;
	
	int isfilter=getBoolParam(argc,argv,"-f");
	int ishist=getBoolParam(argc,argv,"-h");
	int isstretch=getBoolParam(argc,argv,"-s");
	int isequa=getBoolParam(argc,argv,"-e");
	
	char *filepath=getStrParam(argc,argv,"-i","");

	int * valueshisto;
	gray * image_int;

	if(ishist|isstretch|isequa) image_int=readimage_int(filepath);

	if(isfilter){
		int bin=getIntParam(argc,argv,"-s","0");
		int nr=getIntParam(argc,argv,"-n","1");
		char *method=getStrParam(argc,argv,"-f","binomial");
		image=readimage(filepath);
		if(strcmp(method,"binomial")==0){
			kernel=binomialfilter(bin);

			//applying the filter n times
			for(int x=0;x<nr;x++){
				image=ApplyConvolution(bin, kernel, image, lrows, lcols);
			}
			printimage(image,lcols,lrows,lmaxval);
				
		}else if(strcmp(method,"median")==0){
			double *nn;
			nn=(double *)malloc(sizeof(double)*lcols*lrows);

			//applying the filter n times
			for(int x=0;x<nr;x++){			

				//making a copy of the image
				for(int x=0;x<lcols;x++)
				for(int y=0;y<lrows;y++)
					nn[y*lcols+x]=image[y*lcols+x];

				medianFilter(image,nn, lcols,lrows);

				//replacing the image
				for(int x=0;x<lcols;x++)
				for(int y=0;y<lrows;y++)
					image[y*lcols+x]=nn[y*lcols+x];

			}

			printimage(nn,lcols,lrows,lmaxval);
		}else if(strcmp(method,"gradienty")==0){
			kernel=sobelfilterH(3);
			double *resultImage;
			//resultImage=(double *)malloc(sizeof(double)*lcols*lrows);
			//for(int x=0;x<nr;x++){
				resultImage=ApplyConvolution(3, kernel, image, lrows, lcols);
			//}
			printimage(resultImage,lcols,lrows,lmaxval);				
		}else if(strcmp(method,"gradientx")==0){
			kernel=sobelfilterV(3);
			double *resultImage;
			//resultImage=(double *)malloc(sizeof(double)*lcols*lrows);
			//for(int x=0;x<nr;x++){
				resultImage=ApplyConvolution(3, kernel, image, lrows, lcols);
			//}
			printimage(resultImage,lcols,lrows,lmaxval);				
		}else if(strcmp(method,"gradientxy")==0){
			kernel=sobelfilterV(3);
			kernel2=sobelfilterH(3);
			double *gx, *gy, *resultImage;
			gx=ApplyConvolution(3, kernel, image, lrows, lcols);
			gy=ApplyConvolution(3, kernel2, image, lrows, lcols);

			//temp code
			/*double* tmpImage = (double*)(malloc(sizeof(double)*9));
			kernel[0*3+0]=1;
		kernel[0*3+1]=2;
		kernel[0*3+2]=1;

		kernel[1*3+0]=2;
		kernel[1*3+1]=4;
		kernel[1*3+2]=2;

		kernel[2*3+0]=1;
		kernel[2*3+1]=2;
		kernel[2*3+2]=1;*/

			//
			resultImage=normGradient(gx, gy, 3);
			printimage(resultImage,lcols,lrows,lmaxval);
			free(gx);
			free(gy);				
		}
	free(image);
	free(kernel);
	}else if(ishist){
		valueshisto=generateHistogramValues("\nHistogram's values\n\n",image_int,1);
	}else if(isstretch){
		stretchingHistogram("\nStretched Histogram's values\n\n", 0,255,image_int);
	}else if(isequa){
		valueshisto=generateHistogramValues("\nHistogram's values\n\n",image_int,0);
		equalizingHistogram("\nEqualized Histogram's values\n\n",image_int,valueshisto);
	}else{
		printf("Usage: cmd -i -[fsn|h|e|s]\n");	
		printf("OPTIONS: \n");	
		printf("-i S: input file\n");	
		printf("-f [median|binomial|gradientx|gradienty|gradientxy]: chooser the filter\n");	
		printf("-s N: size of the kernel\n");	
		printf("-n N: number of times that the filter will be applyed\n");	
		printf("-h: Histogram\n");	
		printf("-e: Equalization\n");	
		printf("-s: Stretching\n");
		exit(0);	
	}


        if(ishist|isstretch|isequa) free(image_int);
	return 0;
}
