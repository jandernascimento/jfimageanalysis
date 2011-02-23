#include <stdlib.h>
#include <stdio.h>
#include "Util.h"


/** Code for histogram:raquel **/

/** Code for histogram:end **/

/** Code for filtering:jander **/
void setPixel(char *image,int width, int height, int x,int y, int value){

  if(x>width || y>height){
  	printf("Out of range x:%i y:%i",x,y);
	exit(1);
  }

  image[width*y+x]=value;

}

int getPixel(char *image,int width, int height, int x,int y){

  if(x>width || y>height){
  	printf("Out of range x:%i y:%i",x,y);
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

int ApplyConvolution(int dim, double* kernel, int width, int height, char *image){
       int j;  // row    index of the current image
       int i;  // column index of the current image
       int jk; // row    index of the kernel;
       int ik; // column index of the kernel;
       int newval[3]; // new colors
       int size = height*width;
       int* tmpImageR = (int*)(malloc(sizeof(int)*size)); // separated by colors
       int* tmpImageG = (int*)(malloc(sizeof(int)*size)); // separated by colors
       int* tmpImageB = (int*)(malloc(sizeof(int)*size)); // separated by colors
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
           newval[0] = 0;
           newval[1] = 0;
           newval[2] = 0;
           for (jk = 0; jk < dim; jk++) {
             for (ik = 0; ik < dim; ik++) {
               int ii = i + ik - kernelCenteri;
               int jj = j + jk - kernelCenterj;
               if ((jj >= 0) && (jj <height ) && (ii >= 0) && (ii < width))
               {
                   // DO ME
 		   //newval[0] += getPixelColorIntensity(ImageAbstraction::red,jj,ii) * (double)(kernel[jk*dim+ik]);
                   //newval[1] += getPixelColorIntensity(ImageAbstraction::green,jj,ii) * (double)(kernel[jk*dim+ik]);
                   //newval[2] += getPixelColorIntensity(ImageAbstraction::blue,jj,ii) * (double)(kernel[jk*dim+ik]);
               }
             }
           }
           newval[0] = newval[0] / kernelTotalValue;
           newval[1] = newval[1] / kernelTotalValue;
           newval[2] = newval[2] / kernelTotalValue;

           tmpImageR[i*height+j] = newval[0];
           tmpImageG[i*height+j] = newval[1];
           tmpImageB[i*height+j] = newval[2];
         }
       }

       minMax(tmpImageR,findMin(tmpImageR,size),findMax(tmpImageR,size),0,255,size);
       minMax(tmpImageG,findMin(tmpImageG,size),findMax(tmpImageG,size),0,255,size);
       minMax(tmpImageB,findMin(tmpImageB,size),findMax(tmpImageB,size),0,255,size);
       for (int i=0; i<height;++i)
            for (int j=0;j<width;++j)
                //DO ME
		//setPixel(i,j,tmpImageR[j*this->height()+i],tmpImageG[j*this->height()+i],tmpImageB[j*this->height()+i]);

       return 1;
}


/** Code for filtering:end **/

int main(int argc, char* argv[]){

}
