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

// METHODS
pimage_type readimage(char* filepath);
void printimage(pimage_type image);
pixel_type *get_pixel(pimage_type image,int x, int y);
int pixel_distance(pimage_type image,pixel_type *p1,pixel_type *p2);
pixel_type *pixel_middle(pixel_type *p1,pixel_type *p2);
pimage_type assign_to_group(pimage_type image, pixel_type *keys, int size);

