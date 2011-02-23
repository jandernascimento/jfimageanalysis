#include <stdlib.h>
#include <stdio.h>
#include "Util.h"


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

void generateHistogramValues(char * title, gray* graymap){
	//printf("\ndimensions: rows:%d columns:%d\n",rows,cols);

	int valueshisto[256];
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
	//printArray(title, rows, cols, valueshisto);

}

gray * stretchingHistogram(char * title,double newRangeMin, double newRangeMax, gray* graymap){
	double actualRangeMin=(double) getMin(graymap);
	double actualRangeMax=(double) getMax(graymap);

	//printf("\nactualmax:%f actualmin:%f newRangeMax:%f newRangeMin:%f",actualRangeMax,actualRangeMin,newRangeMax,newRangeMin);

	double constant = ( (newRangeMax-newRangeMin) / (actualRangeMax-actualRangeMin) );

	//printf("\nconstant:%3.4f\n",constant);

	gray* graymap_strechted;
	graymap_strechted = (gray *) malloc(cols * rows * sizeof(gray));

	for(int i=0; i < rows; i++)
		for(int j=0; j < cols ; j++){
			graymap_strechted[i * cols + j] = newRangeMin + (constant*(graymap[i * cols + j]-actualRangeMin) );
		}

	//generate histogram
	generateHistogramValues(title, graymap_strechted);
	
	return graymap_strechted; 
	//free(graymap_strechted);
}

/** Code for histogram:end **/

/** Code for filtering:jander **/

void printkernel(int dim, double *kernel){
 for(int i=0;i<3;i++){
  	for(int j=0;j<3;j++){
		printf("%f ",kernel[dim*i+j]);
	}
	printf("\n");
  }
}

void setPixel(char *image,int width, int height, int x,int y, int value){

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

int findMax(int* array, int len)
{
    int max = array[0];
    for (int i=1;i<len; ++i)
        if(array[i]>max)
            max = array[i];
    return max;
}
int findMin(int* array, int len)
{
    int min = array[0];
    for (int i=1;i<len; ++i)
        if (array[i]<min)
            min = array[i];
    return min;
}
void minMax(int* oldArr, int oldMin, int oldMax, int newMin, int newMax, int len)
{
    double tmp;
    for (int i=0;i<len; ++i)
    {
        tmp = (((oldArr[i]-oldMin)/(double)(oldMax-oldMin))*(newMax-newMin))+newMin;
        oldArr[i] = (int)tmp;
    }
}

int ApplyConvolution(int dim, double* kernel, int width, int height, gray *image){
       int j;  // row    index of the current image
       int i;  // column index of the current image
       int jk; // row    index of the kernel;
       int ik; // column index of the kernel;
       int newval[3]; // new colors
       int size = height*width;
       //int* tmpImageR = (int*)(malloc(sizeof(int)*size)); // separated by colors
       int* tmpImageG = (int*)(malloc(sizeof(int)*size)); // separated by colors
       //int* tmpImageB = (int*)(malloc(sizeof(int)*size)); // separated by colors
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
       for (j = 0; j < height; j++) {
         for (i = 0; i < width; i++) {
           //newval[0] = 0;
           newval[1] = 0;
           //newval[2] = 0;
           for (jk = 0; jk < dim; jk++) {
             for (ik = 0; ik < dim; ik++) {
               int ii = i + ik - kernelCenteri;
               int jj = j + jk - kernelCenterj;
               if ((jj >= 0) && (jj <height ) && (ii >= 0) && (ii < width))
               {
                   // DO ME
 		   //newval[0] += getPixelColorIntensity(ImageAbstraction::red,jj,ii) * (double)(kernel[jk*dim+ik]);
                   //newval[1] += getPixelColorIntensity(ImageAbstraction::green,jj,ii) * (double)(kernel[jk*dim+ik]);
                   newval[1] += ((double)getPixel(image,width,height,ii,jj)) * 
					((double)kernel[jk*dim+ik]);
                   //newval[2] += getPixelColorIntensity(ImageAbstraction::blue,jj,ii) * (double)(kernel[jk*dim+ik]);
               }
             }
           }
           //newval[0] = newval[0] / kernelTotalValue;
           newval[1] = newval[1] / kernelTotalValue;
           //newval[2] = newval[2] / kernelTotalValue;

           //tmpImageR[i*height+j] = newval[0];
           tmpImageG[j*width+i] = newval[1];
           //tmpImageB[i*height+j] = newval[2];
         }
       }

       //minMax(tmpImageR,findMin(tmpImageR,size),findMax(tmpImageR,size),0,255,size);
       //minMax(tmpImageG,findMin(tmpImageG,size),findMax(tmpImageG,size),0,255,size);
       //minMax(tmpImageB,findMin(tmpImageB,size),findMax(tmpImageB,size),0,255,size);
       for (int i=0; i<height;++i)
            for (int j=0;j<width;++j){
                //DO ME
		//setPixel(i,j,tmpImageR[j*this->height()+i],tmpImageG[j*this->height()+i],tmpImageB[j*this->height()+i]);
		//setPixel(image,width,height,j,i,getPixel(tmpImageG,width,height,j,i));		
		setPixel(image,width,height,j,i,tmpImageG[i*width+j]);
		}
       return 1;
}

void printimage(gray* image,int cols, int rows, int maxval){

printf("P2\n");
    printf("%d %d \n", cols, rows);
    printf("%d\n",maxval);

    for(int i=0; i < rows; i++)
      for(int j=0; j < cols ; j++)
	printf("%i ",image[i * cols + j] );

}

gray* readimage(int argc, char* argv[])
    {
    FILE* ifp;
    gray* graymap;
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
    graymap = (gray *) malloc(cols * rows * sizeof(gray));

    
    if(pgmraw){
    /* Convertendo para Plaintext */
    /* Lecture */
	    for(i=0; i < rows; i++)
	      for(j=0; j < cols ; j++)
		graymap[i * cols + j] = pm_getint(ifp);

    }else { 
        printf("Invalid file type");
        exit(1);
    }
      /* fermeture */
      fclose(ifp);
      return graymap;
}


/** Code for filtering:end **/

int main(int argc, char* argv[]){
gray* image;
gray* image_strectched;

double *kernel;

kernel=(double *)malloc(sizeof(double)*9);

kernel[0*3+0]=1;
kernel[0*3+1]=2;
kernel[0*3+2]=1;

kernel[1*3+0]=2;
kernel[1*3+1]=4;
kernel[1*3+2]=2;

kernel[2*3+0]=1;
kernel[2*3+1]=2;
kernel[2*3+2]=1;

image=readimage(argc,argv);

//ApplyConvolution(3, kernel, cols, rows, image);

//printimage(image,cols,rows,maxval);
 
//****** RAQUEL *************//
//generate histogram
generateHistogramValues("\nHistogram's values\n\n",image);

//write the image on the screen
//displayImageFile(image);

//histogram stretching
image_strectched=stretchingHistogram("\nStretched Histogram's values\n\n", 0,255,image);

//write the image on the screen
displayImageFile(image_strectched);

free(image);
free(image_strectched);
return 0;

}

