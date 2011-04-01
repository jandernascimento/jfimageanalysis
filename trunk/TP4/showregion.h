#include "Util.h"
#define DPC 3 //data per cell information
// STRUCTURE

typedef struct tpixel {

  int r,g,b;
  int x,y;

} pixel_type; //, *ppixel_type;

typedef struct timage {

  int cols,rows;
  int maxval;
  gray* stream;

} image_type, *pimage_type;

// METHODS
pimage_type readimage(char* filepath);
void printimage(pimage_type image);
pixel_type *get_pixel(pimage_type image,int x, int y);
int pixel_distance(pimage_type image,pixel_type *p1,pixel_type *p2);
pixel_type *pixel_middle(pixel_type *p1,pixel_type *p2);
void assign_to_group(pimage_type image, pixel_type *keys, int size);
