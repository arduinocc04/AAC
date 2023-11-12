#include <cstdlib>
#include <iostream>

#include <png.h>
struct Pixel {
  png_byte r, g, b, a;
};

int getGrayScaledImage(FILE* inputImage, ) {
  if (!inputImage) {
    std::cout << "Open Image Failed." << std::endl;
    return -2;
  }
  unsigned char* header = (unsigned char*)malloc(sizeof(unsigned char) * 8);
  fread(header, 1, 8, inputImage);
  bool isPng = !png_sig_cmp(header, 0, 8);
  free(header);
  if (!isPng) {
    std::cout << "Input File isn't png" << std::endl;
    return -3;
  }

  png_structp pngPtr =
      png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (pngPtr != NULL) {
    png_infop infoPtr = png_create_info_struct(pngPtr);
    if (infoPtr != NULL) {
      if (setjmp(png_jmpbuf(pngPtr)) == 0) {
        png_uint_32 width, height;
        int bitDepth, colorType, interlaceMethod, compressionMethod,
            filterMethod;

        png_init_io(pngPtr, inputImage);
        png_set_sig_bytes(pngPtr, 8);
        png_read_info(pngPtr, infoPtr);

        if (png_get_IHDR(pngPtr, infoPtr, &width, &height, &bitDepth,
                         &colorType, &interlaceMethod, &compressionMethod,
                         &filterMethod)) {
          struct Pixel *rowPointers[height];
          for(int i = 0; i < height; i++) {
            rowPointers[i] = (Pixel*)malloc(width*sizeof(Pixel));
          }
          png_read_image(pngPtr, (png_bytepp)rowPointers);
          for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                png_byte gray = (rowPointers[i][j].r + rowPointers[i][j].g + rowPointers[i][j].b)/3;
                rowPointers[i][j].r = gray;
                rowPointers[i][j].g = gray;
                rowPointers[i][j].b = gray;
            }
          }
        }
        else
          png_error(pngPtr, "pngpixel: png_get_IHDR failed");
      }
    }
  }
  return 0;
}
