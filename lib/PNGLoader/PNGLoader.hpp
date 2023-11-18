/**
 * @file PNGLoader.hpp
 *
 * @brief Load png image easy. Easy pnG
 * @author Daniel Cho
 * @date 2023.11.12
 * @version 0.0.1
 */
#include <string>
#include <cmath>
#include <algorithm>

#include <png.h>

#ifndef __TENSOR
#include "unsupported/Eigen/CXX11/Tensor"
#endif
#define __TENSOR

#ifndef __PNGEXCEPTIONS_H
#include "PNGExceptions.hpp"
#endif
#define __PNGEXCEPTIONS_H

#ifndef __PNGMATH_H
#include "PNGMath.hpp"
#endif
#define __PNGMATH_H

#include "PNGMethods.hpp"

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

using Image = Eigen::Tensor<png_byte, 3>;
/**
 * @class PNG
 * @author Daniel Cho
 * @date 2023.11.12
 * @version 0.0.1
 */
class PNG {
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
     * @attention you need open image before call this function. Call this method will change playground and opened image.
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
     * @brief get Edge of playground.
     * @attention This will change playground and opened image.
     */
    void getEdge(int method);

    /**
     * @brief blur playground.
     * @attention This will change playground and opened image
     */
    void blur(int method);
    void binary(double threshold);
    Eigen::Tensor<double, 2> * getPlayground() {
        return &playground;
    }

    void dividePlaygroundByLength(int _gridHeight, int _gridWidth);

    /**
     * @todo implement correctly
     */
    void dividePlaygroundByCnt(int _gridRowCnt, int _gridColCnt);
    int getGridColCnt() { return gridColCnt; }
    int getGridRowCnt() { return gridRowCnt; }
    int getGridWidth() { return gridWidth; }
    int getGridHeight() { return gridHeight; }
    Eigen::Tensor<double, 2> getPlaygroundAtGrid(int r, int c);

    std::string getInputPath() { return inputPath; }
private:
    Image image;
    Eigen::Tensor<double, 2> playground;
    Eigen::Tensor<double, 2> ** playgroundGrid;
    int gridWidth, gridHeight, gridColCnt, gridRowCnt;
    std::string inputPath;
    std::string outputPath;
    FILE * fimage;
    png_structp pngStructp;
    png_infop pngInfop;
    struct Pixel ** buffer;

    bool isPNG();
    bool getMetadata();

    /**
     * @attention to call this function, buffer must nullptr
     * @see eg::PNG::freeBuffer
     */
    void allocBuffer();
    void freeBuffer();

    /**
     * @attention calling this function will overwrite data
     */
    void allocPlayground();

    /**
     * @attention calling this function will overwrite data
     */
    void allocImage();

    /**
     * @attention to call this function, you must allocate image
     * @see eg::PNG::allocImage
     */
    void copyBufferToImage();

    /**
     * @attention to call this function, you must allocate buffer
     * @see eg::PNG::allocBuffer
     */
    void copyImageToBuffer();

    /**
     * @attention to call this function, you must allocate playground
     * @see eg::PNG::allocPlayground
     */
    void copyImageToPlayground();

    /**
     * @attention to call this function, you must allocate image
     * @see eg::PNG::allocImage
     */
    void copyPlaygroundToImage();

    /**
     * @attention this will overwrite inputPath, fimage, pngStructp, pngInfop, info, buffer
     * @attention After call this, please copy buffer to image
     * @see eg::PNG::copyBufferToImage
     */
    void readImageBuffer(std::string _inputPath);

    /**
     * @brief convert rgba image to grayscale by average pixel values
     */
    void cvtGrayMean();
    void getEdgeGrad();
    void blurGaussian();
};

}

