#include <stdio.h>
#include <stdlib.h>

#define PNG_SETJMP_NOT_SUPPORTED
#include <png.h>

#define WIDTH 256
#define HEIGHT 32
#define COLOR_DEPTH 8

struct Pixel {
	png_byte r, g, b, a;
};

int main(int argc, char *argv[]) {
	srand(time(NULL));

	/* open PNG file for writing */
	FILE *f = fopen("out.png", "wb");
	if (!f) {
		fprintf(stderr, "could not open out.png\n");
		return 1;
	}

	/* initialize png data structures */
	png_structp png_ptr;
	png_infop info_ptr;

	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr) {
		fprintf(stderr, "could not initialize png struct\n");
		return 1;
	}

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
		fclose(f);
		return 1;
	}

	/* begin writing PNG File */
	png_init_io(png_ptr, f);
	png_set_IHDR(png_ptr, info_ptr, WIDTH, HEIGHT, COLOR_DEPTH,
	             PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
	             PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_ptr, info_ptr);

	/* allocate image data */
	struct Pixel *row_pointers[HEIGHT];
	for (int row = 0; row < HEIGHT; row++) {
		row_pointers[row] = (Pixel *)calloc(WIDTH*2, sizeof(struct Pixel));
	}

	/* draw a bunch of vertical lines */
	for (int col = 0; col < WIDTH; col++) {
		int bar_height = rand() % HEIGHT;
		for (int row = 0; row < HEIGHT; row++) {
			if (HEIGHT - row <= bar_height) {
				row_pointers[row][col].r = 255; // red
				row_pointers[row][col].g = 0; // green
				row_pointers[row][col].b = 0; // blue
				row_pointers[row][col].a = 255; // alpha (opacity)
			} else {
				row_pointers[row][col].r = 0; // red
				row_pointers[row][col].g = 0; // green
				row_pointers[row][col].b = 0; // blue
				row_pointers[row][col].a = 0; // alpha (opacity)
			}
		}
	}

	/* write image data to disk */
	png_write_image(png_ptr, (png_byte **)row_pointers);

	/* finish writing PNG file */
	png_write_end(png_ptr, NULL);

	/* clean up PNG-related data structures */
	png_destroy_write_struct(&png_ptr, &info_ptr);

	/* close the file */
	fclose(f);
	f = NULL;

	return 0;
}
