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

double* ApplyConvolution(int dim, double* kernel, double* image, int imageH, int imageW);
double * calc_xx(double *gx, int lrows, int lcols);
double * calc_xy(double *gx, double *gy, int lrows, int lcols);
double * calc_yy(double *gy, int lrows, int lcols);
double *callgradx(int times, int kernelsize, double *img, int height, int width);
double *callgradxy(int times, int kernelsize, double *img, int height, int width);
double *callgrady(int times, int kernelsize, double *img, int height, int width);
double *harrispart(double *image,int rows, int cols);
double *imagediff(double *image1,double *image2,int cols,int rows);
double *imageextract(double *image,int lcols,int lrows,int x, int y, int xn,int ym);
double * interpolate_image(double *image, int rows, int cols, double delta_i, double delta_j);
double* matrixinverse(double* matrix,int dim);
double * mult_matrices(double *mat1, double *mat2);
double * mult_pixels_matrices(double *mat1, double *mat2, int rows, int cols);
double * normGradient(double *gx, double *gy, int lrows, int lcols);
double* readimage(char* filepath);
double* sobelfilterH(int val);
double* sobelfilterV(int val);
double sum_matrices_values(double *matrix,int dim);
int findMax(double* array, int len);
int findMin(double* array, int len);
int lcols, lrows, lmaxval;
int main(int argc, char* argv[]);
void exercise3();
void imageextracttest();
void initialize_array(double * vec, int n, int vldef);
//void matrixprint(double* matrix,int dim);
void matrixprint(double* matrix,int cols,int rows);
void matrixtest();
void minMax(double* oldArr, int oldMin, int oldMax, int newMin, int newMax, int len);
void printimage(double* image,int cols, int rows, int maxval);
double * calc_displacements(double * T, double * Io, int rows, int cols);
