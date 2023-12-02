/**
 * @file egMethods.hpp
 * @author Daniel Cho
 * @version 0.0.1
 */
namespace eg {

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

enum contourMethod {
    suzuki
};

enum distanceMethod {
    taxi
};

enum centerlineMethod {
    grassfire
};

enum paddingMethod {
    zero
};

enum histCmpMethod {
    bhattacharyya
};

enum matCmpMethod {
    rmse, shape, logpolar
};

enum pathDecomMethod {
    greedy, potrace
};

enum resizeMethod {
    vector
};

}
