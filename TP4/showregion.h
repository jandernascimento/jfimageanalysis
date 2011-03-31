#include "Util.h"
// STRUCTURE
struct pixel {

  int x,y;

};

typedef struct timage {

  int cols,rows;
  int maxval;
  gray* stream;

} image_type, *pimage_type;

// METHODS
pimage_type readimage(char* filepath);
void printimage(pimage_type image);
