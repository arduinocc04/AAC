#include "PNGLoader.hpp"
#define PNG_HEAD_BYTE 4

eg::PNG::PNG() {
    info.initialized = false;
    pngStructp = nullptr;
    pngInfop = nullptr;
};

eg::PNG::~PNG() {
    if(image) {
        if(info.initialized) {
            for(int i = 0; i < info.height; i++) {
                delete image[i];
            }
        }
        delete image;
    }

    if(fimage) fclose(fimage);

    if(pngStructp)
        png_destroy_read_struct(&pngStructp,
                                NULL, NULL
        );
};

bool eg::PNG::getMetadata() {
    bool success = png_get_IHDR(pngStructp, pngInfop,
                     &info.width,
                     &info.height, &info.bitDepth,
                     &info.colorType,
                     &info.interlaceMethod,
                     &info.compressionMethod,
                     &info.filterMethod);
    info.initialized = success;
    return success;
}

bool eg::PNG::isPNG() {
    unsigned char * header =
        new unsigned char[PNG_HEAD_BYTE];

    fread(header, sizeof(unsigned char),
          PNG_HEAD_BYTE, fimage);

    bool ans = !png_sig_cmp(header, 0, PNG_HEAD_BYTE);
    delete header;

    return ans;
}

void eg::PNG::openImage(std::string _inputPath) {
    inputPath = _inputPath;
    fimage = fopen(inputPath.c_str(), "rb");

    if(!fimage) throw exceptions::FileNotFound();

    if(!isPNG()) throw exceptions::InvalidFormat();

    pngStructp = png_create_read_struct(
                               PNG_LIBPNG_VER_STRING,
                               nullptr, nullptr,
                               nullptr
                 );
    if(!pngStructp)
        throw exceptions::StructpGenFailed();

    pngInfop = png_create_info_struct(pngStructp);
    if(!pngInfop) throw exceptions::InfopGenFailed();

    if(setjmp(png_jmpbuf(pngStructp)))
        throw exceptions::SetjmpFailed();

    png_init_io(pngStructp, fimage);
    png_set_sig_bytes(pngStructp, PNG_HEAD_BYTE);
    png_read_info(pngStructp, pngInfop);

    if(!getMetadata())
        throw exceptions::GetMetadataFailed();

    image = new Pixel * [info.height];
    for(int i = 0; i < info.height; i++) {
        image[i] = new Pixel[info.width];
    }
    png_read_image(pngStructp, (png_bytepp)image);

    if(!image) throw exceptions::ReadImageFailed();
}

void eg::PNG::cvtGrayMean() {
    if(!info.initialized || !image)
        throw exceptions::ImageNotOpened();
    for(int i = 0; i < info.height; i++) {
        for(int j = 0; j < info.width; j++) {
            Pixel * p = &image[i][j];
            png_byte mean = (p->r + p->g + p->b)/3; // May overflow? ceil? floor? round?
            p->r = p->g = p->b = mean;
        }
    }
}

void eg::PNG::cvtGray(int method) {
    if(!info.initialized)
        throw exceptions::GetMetadataFailed();
    if(info.colorType != PNG_COLOR_TYPE_RGB_ALPHA)
        throw exceptions::InvalidFormat();
    switch(method) {
        case grayCvtMethod::mean:
            cvtGrayMean();
            break;
        default:
            throw exceptions::InvalidParameter();
            break;
    }
}

void eg::PNG::saveImage(std::string _outputPath) {
    outputPath = _outputPath;
    FILE * foutImage = fopen(outputPath.c_str(), "wb");

    if(!foutImage) exceptions::FileNotFound();

    png_structp opngStructp =
        png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                NULL, NULL, NULL);
    if(!opngStructp)
        throw exceptions::StructpGenFailed();

    png_infop oInfop =
        png_create_info_struct(opngStructp);
    if(!oInfop) throw exceptions::InfopGenFailed();

    png_init_io(opngStructp, foutImage);

    if(!info.initialized)
        throw exceptions::GetMetadataFailed();
    png_set_IHDR(opngStructp, oInfop, info.width,
                 info.height, info.bitDepth,
                 info.colorType, info.interlaceMethod,
                 info.compressionMethod,
                 info.filterMethod);
    png_write_info(opngStructp, oInfop);
    png_write_image(opngStructp, (png_byte **)image);
    png_write_end(opngStructp, NULL);
    png_destroy_write_struct(&opngStructp, &oInfop);

    fclose(foutImage);
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
