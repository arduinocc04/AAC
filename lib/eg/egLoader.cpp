/**
 * @file PNGLoader.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include "egLoader.hpp"
#define PNG_HEAD_BYTE 8
#define ALPHA_THRESHOLD 245

eg::PNG::PNG() {
    info.initialized = false;
    pngStructp = nullptr;
    buffer = nullptr;
    pngInfop = nullptr;
};

eg::PNG::~PNG() {
    if(pngStructp) {
        if (pngInfop) {
            png_destroy_read_struct(&pngStructp,
                                    &pngInfop, NULL
            );
        } else {
            png_destroy_read_struct(&pngStructp,
                                    NULL, NULL
            );
        }
    }
};

bool eg::PNG::getMetadata() {
    bool success = png_get_IHDR(pngStructp, pngInfop,
                     &info.width,
                     &info.height, &info.bitDepth,
                     &info.colorType,
                     &info.interlaceMethod,
                     &info.compressionMethod,
                     &info.filterMethod);
    return info.initialized = success;
}

bool eg::PNG::isPNG() {
    unsigned char * header =
        new unsigned char[PNG_HEAD_BYTE];

    fread(header, sizeof(unsigned char),
          PNG_HEAD_BYTE, fimage);

    bool ans = !png_sig_cmp(header, 0, PNG_HEAD_BYTE);
    delete[] header;

    return ans;
}

void eg::PNG::setImage(eg::Image & a) {
    int h = image.dimensions()[0];
    int w = image.dimensions()[1];
    if(h != a.dimensions()[0] || w != a.dimensions()[1])
        throw eg::exceptions::InvalidParameter();

    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            for(int k = 0; k < 4; k++)
                image(i, j, k) = a(i, j, k);
}

void eg::PNG::allocBuffer() {
    if(buffer) throw exceptions::BufferNotNull();
    buffer = new Pixel * [info.height];
    for(int i = 0; i < info.height; i++) {
        buffer[i] = new Pixel[info.width];
    }
}

void eg::PNG::freeBuffer() {
    for(int i = 0; i < info.height; i++) delete[] buffer[i];
    delete[] buffer;
    buffer = nullptr;
}

void eg::PNG::allocImage() {
    image = Eigen::Tensor<png_byte, 3>(info.height, info.width, 4);
}

void eg::PNG::copyBufferToImage() {
    for(int i = 0; i < info.height; i++) {
        for(int j = 0; j < info.width; j++) {
            image(i, j, 0) = buffer[i][j].r;
            image(i, j, 1) = buffer[i][j].g;
            image(i, j, 2) = buffer[i][j].b;
            image(i, j, 3) = buffer[i][j].a;
        }
    }
}

void eg::PNG::copyImageToBuffer() {
    for(int i = 0; i < info.height; i++) {
        for(int j = 0; j < info.width; j++) {
            buffer[i][j].r = image(i, j, 0);
            buffer[i][j].g = image(i, j, 1);
            buffer[i][j].b = image(i, j, 2);
            buffer[i][j].a = image(i, j, 3);
        }
    }
}

void eg::PNG::readImageBuffer(std::string _inputPath) {
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

    allocBuffer();
    png_read_image(pngStructp, (png_bytepp)buffer);

    if(info.colorType == PNG_COLOR_TYPE_RGB) {
        info.colorType = PNG_COLOR_TYPE_RGB_ALPHA;
        for(int i = 0; i < info.height; i++) {
            for(int j = 0; j < info.width; j++) {
                buffer[i][j].a = 255;
            }
        }
    }

    for(int i = 0; i < info.height; i++) {
        for(int j = 0; j < info.width; j++) {
            if(buffer[i][j].a < ALPHA_THRESHOLD) {
                buffer[i][j].r = buffer[i][j].g = buffer[i][j].b = 0;
            }
        }
    }

    fclose(fimage);
}

void eg::PNG::openImage(std::string _inputPath) {
    readImageBuffer(_inputPath);
    allocImage();
    copyBufferToImage();
    freeBuffer();
    if(!image.data()) throw exceptions::ReadImageFailed();
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

    allocBuffer();
    copyImageToBuffer();
    png_write_image(opngStructp, (png_bytepp)buffer);
    freeBuffer();
    png_write_end(opngStructp, NULL);
    png_destroy_write_struct(&opngStructp, &oInfop);

    fclose(foutImage);
}

eg::Image * eg::PNG::getImage() {
    return &image;
}

eg::Image eg::PNG::copy() {
    Image res = Eigen::Tensor<png_byte, 3>(info.height, info.width, 4);
    for(int i = 0; i < info.height; i++)
        for(int j = 0; j < info.width; j++)
            for(int k = 0; k < 4; k++)
                res(i, j, k) = image(i, j, k);

    return res;
}

void eg::PNG::divideImageByLength(int _gridHeight, int _gridWidth) {
    gridWidth = _gridWidth, gridHeight = _gridHeight;
    gridRowCnt = (info.height/gridHeight) + ((info.height % gridHeight) > 0);
    gridColCnt = (info.width/gridWidth) + ((info.width % gridWidth) > 0);
}

eg::Image eg::PNG::getImageAtGrid(int r, int c) {
    if(r < 0 || r >= gridRowCnt)
        throw exceptions::InvalidParameter();
    if(c < 0 || c >= gridColCnt)
        throw exceptions::InvalidParameter();

    Eigen::array<Eigen::Index, 3> offsets = {r*gridHeight, c*gridWidth, 0};
    Eigen::array<Eigen::Index, 3> extents = {std::min((int)info.height - r*gridHeight, gridHeight),
                                             std::min((int)info.width - c*gridWidth, gridWidth),
                                             4};

    return image.slice(offsets, extents);
}
