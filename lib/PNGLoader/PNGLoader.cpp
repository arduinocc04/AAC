#include "PNGLoader.hpp"

#define PNG_HEAD_BYTE 8

bool eg::PNG::getMetadata() {
    return png_get_IHDR(pngStructp, pngInfop, &info.width,
                     &info.height, &info.bitDepth,
                     &info.colorType,
                     &info.interlaceMethod,
                     &info.compressionMethod,
                     &info.filterMethod);
}

bool eg::PNG::isPNG() {
    unsigned char * header =
        new unsigned char[PNG_HEAD_BYTE];

    fread(header, sizeof(unsigned char),
          PNG_HEAD_BYTE, fimage);
    return !png_sig_cmp(header, 0, PNG_HEAD_BYTE);
}

void eg::PNG::openImage(std::string _inputPath) {
    fimage = fopen(inputPath.c_str(), "rb");
    inputPath = _inputPath;

    if(!fimage) throw exceptions::FileNotFound;

    if(!isPNG()) throw exceptions::InvalidFormat;

    pngStructp = png_create_read_struct(
                               PNG_LIBPNG_VER_STRING,
                               NULL, NULL, NULL
                 );
    if(!pngStructp)
        throw exceptions::StructpGenFailed;

    pngInfop = png_create_info_struct(pngStructp);
    if(!pngInfop) throw exceptions::InfopGenFailed;

    if(setjmp(png_jmpbuf(pngStructp)) != 0)
        throw exceptions::SetjmpFailed;

    png_init_io(pngStructp, fimage);
    png_set_sig_bytes(pngStructp, PNG_HEAD_BYTE);

    if(!getMetadata())
        throw exceptions::GetMetadataFailed;

    image = new Pixel * [info.height];
    for(int i = 0; i < info.height; i++) {
        image[i] = new Pixel[info.width];
    }
    png_read_image(pngStructp, (png_bytepp)image);

    if(!image) throw exceptions::ReadImageFailed;
}

eg::Image * eg::PNG::getImage() {
    return &image;
}

eg::Image * eg::PNG::copy() {
    Image res = new Pixel * [info.height];
    for(int i = 0; i < info.height; i++) {
        res[i] = new Pixel[info.width];
    }

    for(int i = 0; i < info.height; i++) {
        for(int j = 0; j < info.width; j++) {
            res[i][j].r = image[i][j].r;
            res[i][j].g = image[i][j].g;
            res[i][j].b = image[i][j].b;
            res[i][j].a = image[i][j].a;
        }
    }

    return &res;
}
