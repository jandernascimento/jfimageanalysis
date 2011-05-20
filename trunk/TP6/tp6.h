#include "Util.h"
#define DPC 4   //data per cell information
#define RED 0   //index of the red in DPC
#define GREEN 1 //index of the green in DPC
#define BLUE 2  //index of the blue in DPC
#define K 3     //index of the group k in DPC

// STRUCTURE
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
gray *get_matrix_pixel(gray *matrix,int row, int col,int dim);
//iterate over image
filelist_type readbackgrounds(char *path);
//calculate the mean image
pimage_type calculate_mean_image(filelist_type list);
//calculate difference matrix
gray *matrix_subtraction(gray *matrix1,gray *matrix2,int rows,int cols);
//calculate determinant
int matrix_determinant(gray *matrix,int dim);
//calculate determinant (test)
void matrix_determinant_test();
