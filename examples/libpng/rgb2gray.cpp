#include <cstdlib>
#include <iostream>

#include <png.h>

struct Pixel {
  png_byte r, g, b, a;
};

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "You didn't give input file name." << std::endl;
    return -1;
  }
  std::string inputFile = argv[1];
  std::string outputFile = "gray-" + inputFile;
  FILE* inputImage = fopen(inputFile.c_str(), "rb");

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
          std::cout << "Width: " << width << " Height: " << height << std::endl;
          std::cout << "Color type: " << ((colorType == PNG_COLOR_TYPE_RGB_ALPHA)?"RGBA":"IDK") << std::endl;

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

          FILE *outputImage = fopen(outputFile.c_str(), "wb");
          if (!outputImage) {
            std::cout << "Could not open output file" << std::endl;
            return -4;
          }

          png_structp outPtr;
          outPtr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
          if(!outPtr) {
            std::cout << "Could not initialize output png struct" << std::endl;
            return -5;
          }
          
          png_infop outInfoPtr;
          outInfoPtr = png_create_info_struct(outPtr);
          if(!outInfoPtr) {
            png_destroy_write_struct(&outPtr, (png_infopp)NULL);
            std::cout << "Could not initialize output png info struct" << std::endl;
            return -6;
          }

          png_init_io(outPtr, outputImage);
          png_set_IHDR(outPtr, outInfoPtr, width, height, bitDepth, colorType, interlaceMethod, compressionMethod, filterMethod);
          png_write_info(outPtr, outInfoPtr);
          png_set_rgb_to_gray(outPtr, 1, -1, -1);
          png_write_image(outPtr, (png_byte **)rowPointers);
          png_write_end(outPtr, NULL);
          png_destroy_write_struct(&outPtr, &outInfoPtr);

          fclose(outputImage);
          outputImage = NULL;
        }
        else
          png_error(pngPtr, "pngpixel: png_get_IHDR failed");
      }
    }
  }
  return 0;
}
