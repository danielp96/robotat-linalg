/**
 * @file constants.h
 *
 * Library configuration macros.
 *
 */


// ====================================================================================================
// Constant macro definitions
// ====================================================================================================
#define MAX_ITERATION_COUNT_SVD (30)    /**< Maximum number of iterations for svd_jacobi_one_sided.c */
#define MAX_VEC_SIZE            (100)   /**< Maximum number of elements allowed for a single row vector. */
#define MAX_MAT_SIZE            (MAX_VEC_SIZE*MAX_VEC_SIZE)     /**< Maximum number of elements allowed for a matrix. */
#define MATH_MATRIX_CHECK               /**< Comment this to disable matrix size checking. */
