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
class FileNotFound: public std::exception {};

class InvalidFormat: public std::exception {};

/**
 * @brief failed to create png_structp
 */
class StructpGenFailed: public std::exception {};

/**
 * @brief failed to create png_infop
 */
class InfopGenFailed: public std::exception {};

class SetjmpFailed: public std::exception {};

class GetMetadataFailed: public std::exception {};

class ReadImageFailed: public std::exception {};
}
