/**
 * @file PNGLoader.h
 *
 * @brief Load png image easy. Easy pnG
 * @author Daniel Cho
 * @date 2023.11.12
 * @version 0.0.1
 */
#include <png.h>
#include <string>

#include "unsupported/Eigen/CXX11/Tensor"

#ifndef __PNGEXCEPTIONS_H
#include "PNGExceptions.hpp"
#endif
#define __PNGEXCEPTIONS_H

namespace eg {

struct Pixel {
    png_byte r, g, b, a;
};

struct PNGInfo {
    png_uint_32 width, height;
    int bitDepth, colorType, interlaceMethod,
        compressionMethod, filterMethod;
    bool initialized;
};

/**
 * @enum eg::grayCvtMethod
 */
enum grayCvtMethod {
    mean
};

/**
 * @enum eg::edgeDetectMethod
 */
enum edgeDetectMethod {
    gradient
};

/**
 * @enum eg::blurMethod
 */
enum blurMethod {
    gaussian
};

using Image = Eigen::Tensor<png_byte, 3>;
/**
 * @class PNG
 * @author Daniel Cho
 * @date 2023.11.12
 * @version 0.0.1
 */
class PNG {
private:
    Image image;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> playground;
    std::string inputPath;
    std::string outputPath;
    FILE * fimage;
    png_structp pngStructp;
    png_infop pngInfop;
    struct Pixel ** buffer;

    bool isPNG();
    bool getMetadata();

    void allocBuffer();
    void freeBuffer();
    void allocPlayground();
    void allocImage();
    void copyBufferToImage();
    void copyImageToBuffer();
    void copyImageToPlayground();
    void readImageBuffer(std::string _inputPath);
    /**
     * @brief convert rgba image to grayscale by average pixel values
     */
    void cvtGrayMean();
    void getEdgeGrad();
public:
    PNGInfo info;

    PNG();
    ~PNG();
    /**
     * @brief open png image at given path.
     *
     * @param _inputPath an string
     *
     * @throws FileNotFound
     * @throws InvalidFormat
     * @throws StructpGenFailed
     * @throws InfopGenFailed
     * @throws SetjmpFailed
     * @throws GetMetadataFailed
     * @throws ReadImageFailed
     *
     * @see PNGExceptions.hpp
     */
    void openImage(std::string _inputPath);

    /**
     * @brief get pointer of image
     *
     * @attention this doesn't copy image.
     *
     * @return Image pointer
     */
    Image * getImage();

    /**
     * @brief convert opened image to grayscale.
     * @attention you need open image before call this function. Call this method will change opened image.
     * @param method an integer
     * @see eg::PNG::cvtGrayMean
     * @see eg::grayCvtMethod
     * @throws GetMetadataFailed
     * @throws InvalidFormat
     * @throws InvalidParameter
     */
    void cvtGray(int method);

    /**
     * @brief save opened image.
     * @param _outputPath a string
     * @throws FileNotFound
     * @throws StructpGenFailed
     * @throws InfopGenFailed
     * @throws GetMetadataFailed
     */
    void saveImage(std::string _outputPath);

    /**
     * @brief get pointer of copied image
     *
     * @attention Do not use this function as far as you can. This function can be deprecated because it's hard to free Image type.
     *
     * @see getImage
     *
     * @return Image pointer
     */
    Image * copy();

    /**
     * @brief get Edge of image.
     * @attention This will change opened image.
     */
    void getEdge(int method);
};

}

