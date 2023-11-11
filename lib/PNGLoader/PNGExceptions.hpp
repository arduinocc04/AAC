/**
 * @file PNGExceptions.hpp
 *
 * @brief Exceptions
 * @author Daniel Cho
 * @date 2023.11.12
 * @version 0.0.1
 */

#include <exception>

namespace eg::exceptions { // You need C++17
class FileNotFound: public std::exception {
} FileNotFound;

class InvalidFormat: public std::exception {
} InvalidFormat;

/**
 * @brief failed to create png_structp
 */
class StructpGenFailed: public std::exception {
} StructpGenFailed;

/**
 * @brief failed to create png_infop
 */
class InfopGenFailed: public std::exception {
} InfopGenFailed;

class SetjmpFailed: public std::exception {
} SetjmpFailed;

class GetMetadataFailed: public std::exception {
} GetMetadataFailed;

class ReadImageFailed: public std::exception {
} ReadImageFailed;
}
