/**
 * @file matf32_check.c
 */

#include <stdint.h>
#include <stdbool.h>

#include "matf32/matf32_def.h"
#include "matf32/matf32_check.h"
#include "matf32/math_util.h"


// ====================================================================================================
// Matrix datatype-based checks
// ====================================================================================================


bool
matf32_check_triangular_upper(const matf32_t* const p_mat)
{
#ifdef MATH_MATRIX_CHECK
    if (p_mat->num_rows != p_mat->num_cols)
    {
        return false;
    }
#endif

    float* p_data_src = p_mat->p_data;

    for (uint16_t i = 1; i < p_mat->num_rows; ++i)
    {
        for (uint16_t j = 0; j < i; ++j)
        {
            if (0 != p_data_src[i*p_mat->num_rows + j])
            {
                return false;
            }
        }
    }

    return true;
}

bool
matf32_check_triangular_lower(const matf32_t* const p_mat)
{
#ifdef MATH_MATRIX_CHECK
    if (p_mat->num_rows != p_mat->num_cols)
    {
        return false;
    }
#endif

    float* p_data_src = p_mat->p_data;

    for (uint16_t i = 0; i < p_mat->num_rows-1; ++i)
    {
        for (uint16_t j = p_mat->num_cols-1; j > i; --j)
        {
            if (0 != p_data_src[i*p_mat->num_rows + j])
            {
                return false;
            }
        }
    }

    return true;
}

bool
matf32_is_equal(const matf32_t* const p_mat_a, const matf32_t* const p_mat_b)
{
#ifdef MATH_MATRIX_CHECK
    if ((p_mat_a->num_rows != p_mat_b->num_rows) || (p_mat_a->num_cols != p_mat_b->num_cols))
    {
        return false;
    }
#endif
    return is_equal(p_mat_a->p_data, p_mat_b->p_data, p_mat_a->num_rows * p_mat_a->num_cols);
}
