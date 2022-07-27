/**
 * @file linsolve.c
 */

#include "linsolve.h"



err_status_t
matf32_forward_substitution(const matf32_t* const p_l, float* const p_b, float* p_x)
{
#ifdef MATH_MATRIX_CHECK
    if (!matf32_check_triangular_lower(p_l))
    {
        return MATH_ARGUMENT_ERROR;
    }
#endif

    float* p_data_src = p_l->p_data;

    float lx = 0; // sum accumulator

    for (uint16_t i = 0; i < p_l->num_rows; ++i) 
    {
        // reset lx
        lx = 0;

        // calculate sum x_i * l_(i,j)
        for (uint16_t j = 0; j < i; ++j)
        {
            lx += p_x[j]*p_data_src[i*p_l->num_rows + j];
        }

        // calculate x_i
        p_x[i] = (p_b[i] - lx)/p_data_src[i*p_l->num_rows + i];
    }

    return MATH_SUCCESS;
}


err_status_t
matf32_backward_substitution(const matf32_t* const p_u, float* const p_b, float* p_x)
{
#ifdef MATH_MATRIX_CHECK
    if (!matf32_check_triangular_upper(p_u))
    {
        return MATH_ARGUMENT_ERROR;
    }
#endif

    float* p_data_src = p_u->p_data;

    float ux = 0; // sum accumulator

    for (int16_t i = p_u->num_rows-1; i >= 0; --i) 
    {
        // reset ux
        ux = 0;
        // calculate sum x_i * u_(i,j)
        for (uint16_t j = p_u->num_cols-1; j>i; --j)
        {
            ux += p_x[j]*p_data_src[i*p_u->num_rows + j];
        }

        // calculate x_i
        p_x[i] = (p_b[i] - ux)/p_data_src[i*p_u->num_rows + i];
    }
    return MATH_SUCCESS;
}
