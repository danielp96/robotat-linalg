
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "robotat_linalg.h"

float A_data[] = {0, 1, 2,
                  3, 4, 5,
                  6, 7, 8};

float B_data[] = {0, 1, 2,
                  3, 4, 5,
                  6, 7, 8};

float C_data[] = {0, 0, 0,
                  0, 0, 0,
                  0, 0, 0};

float Upper_data[] = {1, 2, 3,
                      0, 4, 5,
                      0, 0, 6};

float Lower_data[] = {1, 0, 0,
                      2, 3, 0,
                      4, 5, 6};

float b1_vector[] = {1, 8, 32};
float b2_vector[] = {10, 13, 6};
float x_vector[] = {0, 0, 0};

void
main(void)
{
    matf32_t A, B, C, Upper, Lower;

    matf32_init(&A, 3, 3, A_data);
    matf32_init(&B, 3, 3, B_data);
    matf32_init(&C, 3, 3, C_data);

    printf("Matrix A: \n");
    matf32_print(&A);
    printf("Matrix B: \n");
    matf32_print(&B);

    printf("Testing add: \n");
    matf32_add(&A, &B, &C);
    matf32_print(&C);


    printf("Testing sub: \n");
    matf32_sub(&A, &B, &C);
    matf32_print(&C);

    printf("Testing scale: \n");
    matf32_scale(&A, 10, &C);
    matf32_print(&C);

    printf("Testing transpose: \n");
    matf32_trans(&A, &C);
    matf32_print(&C);

    printf("Testing equal: \n");

    if (matf32_is_equal(&A, &B))
    {
        printf("A and B are equal.\n\n");
    }
    else
    {
        printf("A and B are NOT equal.\n\n");
    }


    matf32_init(&Upper, 3, 3, Upper_data);
    matf32_init(&Lower, 3, 3, Lower_data);

    // add test for not square matrix
    printf("Testing upper triangular: \n");
    matf32_print(&Upper);
    if (matf32_check_triangular_upper(&Upper))
    {
        printf("Matrix is upper triangular.\n\n");
    }
    else
    {
        printf("Matrix is NOT upper triangular.\n\n");
    }

    matf32_print(&Lower);
    if (matf32_check_triangular_upper(&Lower))
    {
        printf("Matrix is upper triangular.\n\n");
    }
    else
    {
        printf("Matrix is NOT upper triangular.\n\n");
    }



    printf("Testing lower triangular: \n");
    matf32_print(&Upper);
    if (matf32_check_triangular_lower(&Upper))
    {
        printf("Matrix is lower triangular.\n\n");
    }
    else
    {
        printf("Matrix is NOT lower triangular.\n\n");
    }

    matf32_print(&Lower);
    if (matf32_check_triangular_lower(&Lower))
    {
        printf("Matrix is lower triangular.\n\n");
    }
    else
    {
        printf("Matrix is NOT lower triangular.\n\n");
    }

    printf("Testing forward substitution: \n");
    matf32_print(&Lower);

    printf("vector b: \n");
    for (int i=0; i<3; ++i)
    {
        printf("%f\n", b1_vector[i]);
    }
    printf("\n");

    matf32_forward_substitution(&Lower, b1_vector, x_vector);

    // expected x = {1, 2, 3}
    printf("vector x: \n");
    for (int i=0; i<3; ++i)
    {
        printf("%f\n", x_vector[i]);
    }


    printf("\n\n");
    printf("Testing backward substitution: \n");
    x_vector[0] = 0;
    x_vector[1] = 0;
    x_vector[2] = 0;
    matf32_print(&Upper);

    printf("vector b: \n");
    for (int i=0; i<3; ++i)
    {
        printf("%f\n", b2_vector[i]);
    }
    printf("\n");

    matf32_backward_substitution(&Upper, b2_vector, x_vector);

    // expected x = {3, 2, 1}
    printf("vector x: \n");
    for (int i=0; i<3; ++i)
    {
        printf("%f\n", x_vector[i]);
    }

}

err_status_t
matf32_equal(const matf32_t* p_srca, const matf32_t* p_srcb)
{

}