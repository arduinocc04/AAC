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
#include "PNGExceptions.hpp"

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

using Image = Pixel **;

/**
 * @class PNG
 * @author Daniel Cho
 * @date 2023.11.12
 * @version 0.0.1
 */
class PNG {
private:
    Image image;
    std::string inputPath;
    std::string outputPath;
    FILE * fimage;
    png_structp pngStructp;
    png_infop pngInfop;

    bool isPNG();
    bool getMetadata();
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
     * @todo Implement.
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
};

}

