#include "PNGLoader.hpp"
#define PNG_HEAD_BYTE 4

eg::PNG::PNG() {
    info.initialized = false;
    pngStructp = nullptr;
    pngInfop = nullptr;
};

eg::PNG::~PNG() {
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
    return info.initialized = success;
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

void eg::PNG::allocPlayground() {
    playground = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(info.height, info.width);
}

void eg::PNG::allocBuffer() {
    if(buffer) throw exceptions::BufferNotNull();
    buffer = new Pixel * [info.height];
    for(int i = 0; i < info.height; i++) {
        buffer[i] = new Pixel[info.width];
    }
}

void eg::PNG::freeBuffer() {
    for(int i = 0; i < info.height; i++) delete buffer[i];
    delete buffer;
    buffer = nullptr;
}

void eg::PNG::allocImage() {
    image = Eigen::Tensor<png_byte, 3>(info.height, info.width, 4);
}

void eg::PNG::copyImageToPlayground() {
    for(int i = 0; i < info.height; i++)
        for(int j = 0; j < info.width; j++)
            // playground(i, j) = (image(i, j, 0) + image(i, j, 1) + image(i, j, 2))/3;
            playground(i, j) = image(i, j, 0); // Because we will grayscale image.
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
}

void eg::PNG::openImage(std::string _inputPath) {
    readImageBuffer(_inputPath);
    allocImage();
    copyBufferToImage();
    freeBuffer();
    if(!image.data()) throw exceptions::ReadImageFailed();
}

void eg::PNG::cvtGrayMean() {
    if(!info.initialized || !image.data())
        throw exceptions::ImageNotOpened();

    for(int i = 0; i < info.height; i++) {
        for(int j = 0; j < info.width; j++) {
            int mean = 0;
            for(int k = 0; k < 3; k++)
                mean += image(i, j, k);
            mean /= 3;
            image(i, j, 0) = image(i, j, 1) = image(i, j, 2) = mean;
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

eg::Image * eg::PNG::copy() {
    Image res = Eigen::Tensor<png_byte, 3>(info.height, info.width, 4);
    for(int i = 0; i < info.height; i++)
        for(int j = 0; j < info.width; j++)
            for(int k = 0; k < 4; k++)
                res(i, j, k) = image(i, j, k);

    return &res;
}
