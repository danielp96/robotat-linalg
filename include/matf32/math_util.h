/**
 * @file math_util.h
 *
 * Linear algebra routines that do not depend on the matrix datatype
 *
 */

#ifndef ROBOTAT_MATH_UTIL_H_
#define ROBOTAT_MATH_UTIL_H_


#include <string.h>                     // For memcpy, memset etc.
#include <stdio.h>                      // For printf.
#include <stdlib.h>                     // Standard library.
#include <stdint.h>                     // For uint8_t, uint16_t and uint16_t.
#include <math.h>                       // For sqrtf.
#include <float.h>                      // Required for FLT_EPSILON.
#include <stdbool.h>                    // For bool datatype.
#include <time.h>                       // For srand, clock.


/**
 * @brief   Find the dot product of two vectors, pointed by p_srca and p_srcb, of the same size.
 *
 * @param[in]   p_srca  Points to vector 1.
 * @param[in]   p_srcb  Points to vector 2.
 * @param[in]   length  Length of vectors.
 *
 * @return Dot product between the vectors.
 */
float
dot(float* p_srca, float* p_srcb, uint16_t length);


/**
 * @brief   Create an identity matrix array with size row x column.
 *
 * @param[in, out]  p_dst   Points to array to allocate the identity matrix.
 * @param[in]       row     Number of required rows.
 * @param[in]       column  Number of required columns.
 *
 * @return None.
 */
void
eye(float* p_dst, uint16_t row, uint16_t column);


/**
 * @brief   Creates a diagonal matrix array pointed by p_dst with the size row x column,
 * from a vector pointed by p_src. Notice that the row of vector x need to be the same
 * length as the column of the matrix.
 *
 * @param[in]       p_src       Points to vector with diagonal entries.
 * @param[in,out]   p_dst       Points to array to allocate the diagonal matrix.
 * @param[in]       row_d       Number of required rows.
 * @param[in]       column_d    Number of required columns.
 *
 * @return None.
 */
void
diag(float* p_src, float* p_dst, int row_d, int column_d);


/**
 * @brief   Turn all elements of the matrix array pointed by p_dst, size row x column, into 0.
 *
 * @param[in,out]   p_dst   Points to zero matrix array.
 * @param[in]       row     Number or rows.
 * @param[in]       column  Number of columns.
 *
 * @return None.
 */
void
zeros(float* p_dst, int row, int column);


/**
 * @brief   Turn all elements of the matrix array pointed by p_dst, size row x column, into 1.
 *
 * @param[in,out]   p_dst   Points to ones matrix array.
 * @param[in]       row     Number or rows.
 * @param[in]       column  Number of columns.
 *
 * @return None.
 */
void
ones(float* p_dst, int row, int column);


/**
 * @brief   Creates a random array with values sampled from a normal (Gaussian) distribution.
 *
 * @param[in, out]  p_dst   Points to random vector to create.
 * @param[in]       length  Vector length.
 * @param[in]       mu      Mean.
 * @param[in]       sigma   Standar deviation.
 *
 * @return  None.
 */
void
randn(float* p_dst, uint16_t length, float mu, float sigma);


// ====================================================================================================
// Miscellaneous
// ====================================================================================================
/**
 * @brief   Size-aware matrix copy.
 *
 * @param[in]       p_src   Points to matrix to copy from.
 * @param[in, out]  p_dst   Points to matrix to copy to.
 *
 * @return  Execution status
 *              MATH_SUCCESS :          Operation successful.
 *              MATH_SIZE_MISMATCH :    Matrix size check failed.
 */
void
copy(float* p_src, float* p_dst, int row, int column);


/**
 * @brief   Prints array to console (formatted).
 *
 * @param[in]   p_src   Points to array to print.
 * @param[in]   row     Number of rows.
 * @param[in]   column  Number of columns.
 *
 * @return  None.
 */
void
print(float* p_src, uint16_t row, uint16_t column);


/**
 * @brief   Saturates the input, given upper and lower limits.
 * 
 * @param[in]   input           Input value.
 * @param[in]   lower_limit     Lower saturation threshold.
 * @param[in]   upper_limit     Upper saturation threshold.
 * 
 * @return  Saturated output.
 */
float
saturation(float input, float lower_limit, float upper_limit);


/**
 * @brief   Gets the input's sign.
 *
 * @param[in]   number      Input value.
 * 
 * @return  Sign of input.
 */
float
sign(float number);


/**
 * @brief   Gets the mean of a given array.
 *
 * @param[in]   p_src   Points to input array.
 * @param[in]   length  Array length.
 *
 * @return  Mean of array.
 */
float
mean(float* p_src, uint16_t length);


/**
 * @brief   Gets the standard deviation of a given array.
 *
 * @param[in]   p_src   Points to input array.
 * @param[in]   length  Array length.
 *
 * @return  Standard deviation of array.
 */
float
std(float* p_src, uint16_t length);


/**
 * @brief   Compares two arrays.
 *
 * @param[in]   p_a   Points to first array.
 * @param[in]   p_b   Points to second array.
 *
 * @return  If arrays values are equal.
 */
bool
is_equal(float* p_a, float* p_b, uint16_t length);


#endif // ROBOTAT_MATH_UTIL_H_