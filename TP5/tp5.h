#include "Util.h"
#define DPC 1   //data per cell information
#define RED 0   //index of the red in DPC
#define GREEN 1 //index of the green in DPC
#define BLUE 2  //index of the blue in DPC
#define K 3     //index of the group k in DPC

// STRUCTURE

typedef gray istream, *pistream;

typedef struct tpixel {

  int r,g,b;
  int x,y;
  int k; //number of the group that the pixel belongs to
} pixel_type; //, *ppixel_type;


typedef struct timage {

  int cols,rows;
  int maxval;
  gray* stream;
  char* path;
} image_type, *pimage_type;

// METHODS
double* readimage(char* filepath);
void printimage(double* image,int cols, int rows, int maxval);
void set_pixel(pimage_type image,int x, int y, pixel_type pixel );

double * crop_image(double *image, int x1, int y1, int xn, int xm);

pixel_type *get_pixel(pimage_type image,int x, int y);

double *callgradxy(int times, int kernelsize, double *img, int height, int width);
double *callgrady(int times, int kernelsize, double *img, int height, int width);
double *callgradx(int times, int kernelsize, double *img, int height, int width);
double * normGradient(double *gx, double *gy, int lrows, int lcols);
double* sobelfilterV(int val);
double* sobelfilterH(int val);
double* ApplyConvolution(int dim, double* kernel, double* image, int imageH, int imageW);
int findMax(double* array, int len);
int findMin(double* array, int len);
void minMax(double* oldArr, int oldMin, int oldMax, int newMin, int newMax, int len);

//interpolation
double * interpolate_image(double *image, int rows, int cols, double delta_x, double delta_y);
//
