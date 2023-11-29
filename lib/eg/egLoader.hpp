/**
 * @file egLoader.hpp
 *
 * @brief Load png image easy. Easy pnG
 * @author Daniel Cho
 * @date 2023.11.12
 * @version 0.0.1
 */
#include <string>

#include <png.h>

#ifndef __TENSOR
#include "unsupported/Eigen/CXX11/Tensor"
#endif
#define __TENSOR

#ifndef __EGEXCEPTIONS_H
#include "egExceptions.hpp"
#endif
#define __EGEXCEPTIONS_H

#ifndef __EGTYPES_H
#include "egTypes.hpp"
#endif
#define __EGTYPES_H

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

    void setImage(Image & a);

    /**
     * @brief get pointer of image
     *
     * @attention this doesn't copy image.
     *
     * @return Image pointer
     */
    Image * getImage();

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
     * @brief get copied image
     *
     * @see getImage
     *
     * @return Image
     */
    Image copy();

    void divideImageByLength(int _gridHeight, int _gridWidth);

    int getGridColCnt() { return gridColCnt; }
    int getGridRowCnt() { return gridRowCnt; }
    int getGridWidth() { return gridWidth; }
    int getGridHeight() { return gridHeight; }
    Image getImageAtGrid(int r, int c);

    std::string getInputPath() { return inputPath; }
private:
    Image image;
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
     * @attention this will overwrite inputPath, fimage, pngStructp, pngInfop, info, buffer
     * @attention After call this, please copy buffer to image
     * @see eg::PNG::copyBufferToImage
     */
    void readImageBuffer(std::string _inputPath);
};

}

