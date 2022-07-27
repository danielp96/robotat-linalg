/**
 * @file robotat_linalg.c
 */
#include "matf32/matf32_def.h"


void
matf32_init(matf32_t* const instance, uint16_t num_rows, uint16_t num_cols, float* p_data)
{
    instance->num_rows = num_rows;
    instance->num_cols = num_cols;
    instance->p_data = p_data;
}


void
matf32_print(const matf32_t* p_src)
{
    float* p_data_src = p_src->p_data;

    for (uint16_t i = 0; i < p_src->num_rows; i++)
    {
        for (uint16_t j = 0; j < p_src->num_cols; j++)
        {
            printf("%0.18f\t", *(p_data_src++));
        }
        printf("\n");
    }
    printf("\n");
}


// ====================================================================================================
// Special matrix initializations
// ====================================================================================================


void
matf32_eye(matf32_t* const p_dst)
{
    eye(p_dst->p_data, p_dst->num_rows, p_dst->num_cols);
}


void
matf32_diag(float* p_src, matf32_t* const p_dst)
{
    diag(p_src, p_dst->p_data, p_dst->num_rows, p_dst->num_cols);
}


void
matf32_zeros(matf32_t* const p_dst)
{
    zeros(p_dst->p_data, p_dst->num_rows, p_dst->num_cols);
}


void
matf32_ones(matf32_t* const p_dst)
{
    ones(p_dst->p_data, p_dst->num_rows, p_dst->num_cols);
}


void
matf32_randn(matf32_t* const p_dst, float mu, float sigma)
{
    randn(p_dst->p_data, p_dst->num_rows * p_dst->num_cols, mu, sigma);
}
