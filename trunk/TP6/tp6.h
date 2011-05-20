#include "Util.h"
#define DPC 3   //data per cell information
#define RED 0   //index of the red in DPC
#define GREEN 1 //index of the green in DPC
#define BLUE 2  //index of the blue in DPC

// STRUCTURE
typedef struct tpixel {

  int r,g,b;
  int x,y;
} pixel_type; 


typedef struct timage {

  int cols,rows;
  int maxval;
  gray* stream;
  char* path;
} image_type, *pimage_type;

typedef struct tfilelist {

  int size;
  char **paths;

} filelist_type, *pfilelist_type;

// METHODS
//read image
pimage_type readimage(char* filepath);
//print image
void printimage(pimage_type image);
//get pixel
pixel_type *get_pixel(pimage_type image,int x, int y);
//get matrix position
gray get_matrix_pixel(gray *matrix,int row, int col,int dim);
//iterate over image
filelist_type readbackgrounds(char *path, int n);
//calculate the mean image
pimage_type calculate_mean_image(filelist_type list);
//prints the matrix
void matrix_print(gray* matrix,int rows, int cols);
//calculate the multiplication of a matrix for a single number
gray *matrix_multiplication_single(gray value,gray *matrix,int rows,int cols);
//calculate the inverse of the matrix
gray *matrix_inverse(gray *matrix,int rows,int cols);
//calculate the inverse of the matrix test
void matrix_inverse_test();
//calculate the multiplication of a matrix for a single number test
void matrix_multiplication_single_test();
//calculate difference matrix
gray *matrix_subtraction(gray *matrix1,gray *matrix2,int rows,int cols);
//calculate difference matrix test
void matrix_subtraction_test();
//calculate determinant
gray matrix_determinant(gray *matrix,int dim);
//calculate determinant (test)
void matrix_determinant_test();
