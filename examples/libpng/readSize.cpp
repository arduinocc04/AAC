#include <iostream>
#include <string>
#include <cstdlib>
#include <cstring>

#include <png.h>

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cout << "You didn't give input file name." << std::endl;
		return -1;
	}
	std::string inputFile = argv[1];
	FILE *inputImage = fopen(inputFile.c_str(), "rb");

	if (!inputImage) {
		std::cout << "Open Image Failed." << std::endl;
		return -2;
	}
	unsigned char *header = (unsigned char*)malloc(sizeof( unsigned char)*8);
	fread(header, 1, 8, inputImage);
	bool isPng = !png_sig_cmp(header, 0, 8);
	free(header);
	if(!isPng) {
		std::cout << "Input File isn't png" << std::endl;
		return -3;
	}
	png_image image;
	memset(&image, 0, sizeof(image));
	image.version = PNG_IMAGE_VERSION;
	image.opaque = NULL;
	if (png_image_begin_read_from_file(&image, argv[1])) {
		png_bytep buffer;
		image.format = PNG_FORMAT_RGBA;
		buffer = (png_bytep)malloc(PNG_IMAGE_SIZE(image));

		if (buffer != NULL) {
			if (png_image_finish_read(&image, NULL, buffer, 0, NULL)) {
				std::cout << "Width: " << image.width << std::endl;
				std::cout << "Height: " << image.height << std::endl;
			}
		}
	}
}
