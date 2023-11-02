#include <cstdlib>
#include <iostream>

#include <png.h>

#include <ncurses.h>

struct Pixel {
  png_byte r, g, b, a;
};

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "You didn't give input file name." << std::endl;
    return -1;
  }
  std::string inputFile = argv[1];
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
          png_bytepp gray = (png_bytepp)malloc(sizeof(png_bytep)*height);
          for(int i = 0; i < height; i++) {
            gray[i] = (png_bytep)malloc(sizeof(png_byte)*width);
          }
          for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                gray[i][j] = (rowPointers[i][j].r + rowPointers[i][j].g + rowPointers[i][j].b)/3;
            }
          }
          initscr();
          for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                if(gray[i][j] > 200) {
                    mvprintw(i, j, "#");
                }
                else if(gray[i][j] > 150) {
                    mvprintw(i, j, "o");
                }
                else if(gray[i][j] > 100) {
                    mvprintw(i, j, "^");
                }
                else if(gray[i][j] > 50) {
                    mvprintw(i, j, ".");
                }
            }
            std::cout << std::endl;
          }
          refresh();
          getch();
          endwin();
        }
        else
          png_error(pngPtr, "pngpixel: png_get_IHDR failed");
      }
    }
  }
  return 0;
}
