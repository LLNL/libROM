/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This source file is a test runner that uses the Google Test
// Framework to run unit tests on the CAROM::Vector class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "linalg/Vector.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

/**
 *  Test methods that do not require assigning data to the Vector. These
 *  methods are:
 *
 *  TODO(oxberry1@llnl.gov): Do more exhaustive testing to test all branches
 *  of each method. For now, simple tests are better than nothing.
 *
 *  * Vector(int, bool)
 *
 *  * bool distributed() const
 *
 *  * int dim() const
 *
 *  * void setSize(int)
 *
 */

TEST(VectorSerialTest, Test_distributed)
{
    CAROM::Vector v(2, false);
    EXPECT_FALSE(v.distributed());

    CAROM::Vector w(2, true);
    EXPECT_TRUE(w.distributed());
}

TEST(VectorSerialTest, Test_dim)
{
    CAROM::Vector v(2, false);
    EXPECT_EQ(v.dim(), 2);
}

TEST(VectorSerialTest, Test_setSize)
{
    CAROM::Vector v(2, false);
    EXPECT_EQ(v.dim(), 2);
    v.setSize(3);
    EXPECT_EQ(v.dim(), 3);
}

/** Test methods that require assigning data to the Vector
 *
 *  TODO(oxberry1@llnl.gov): Do more exhaustive testing to test all branches
 *  of each method. For now, simple tests are better than nothing.
 *
 *  * Vector(double*, int, bool, bool)
 *
 *  * const double& operator() (int, int) const
 *
 *  * double& operator() (int, int)
 *
 *  * const double& item(int, int) const
 *
 *  * double& item(int, int)
 *
 *  * Vector(const Vector&)
 *
 *  * Vector& operator=(const Vector&) for copy-assignment
 *
 *  * Vector& operator=(const Vector&) for assignment
 */

TEST(VectorSerialTest, Test_const_call_operator)
{
    double v_data[2] = {1, 2};
    const CAROM::Vector v(v_data, 2, false, true);

    EXPECT_DOUBLE_EQ(v(0), 1);
    EXPECT_DOUBLE_EQ(v(1), 2);
}

TEST(VectorSerialTest, Test_nonconst_call_operator)
{
    CAROM::Vector v(2, false);
    v(0) = 1;
    v(1) = 2;

    EXPECT_DOUBLE_EQ(v(0), 1);
    EXPECT_DOUBLE_EQ(v(1), 2);
}

TEST(VectorSerialTest, Test_const_item)
{
    double v_data[2] = {1, 2};
    const CAROM::Vector v(v_data, 2, false, true);

    EXPECT_DOUBLE_EQ(v.item(0), 1);
    EXPECT_DOUBLE_EQ(v.item(1), 2);
}

TEST(VectorSerialTest, Test_nonconst_item)
{
    CAROM::Vector v(2, false);
    v.item(0) = 1;
    v.item(1) = 2;

    EXPECT_DOUBLE_EQ(v.item(0), 1);
    EXPECT_DOUBLE_EQ(v.item(1), 2);
}

TEST(VectorSerialTest, Test_copy_constructor)
{
    double v_data[2] = {1, 2};
    const CAROM::Vector v(v_data, 2, false, true);
    CAROM::Vector w(v);

    EXPECT_FALSE(w.distributed());
    EXPECT_EQ(w.dim(), 2);
    EXPECT_DOUBLE_EQ(w(0), 1);
    EXPECT_DOUBLE_EQ(w(1), 2);
}

TEST(VectorSerialTest, Test_copy_assignment_operator)
{
    double v_data[2] = {1, 2};
    const CAROM::Vector v(v_data, 2, false, true);
    CAROM::Vector w = v;

    EXPECT_FALSE(w.distributed());
    EXPECT_EQ(w.dim(), 2);
    EXPECT_DOUBLE_EQ(w(0), 1);
    EXPECT_DOUBLE_EQ(w(1), 2);
}

TEST(VectorSerialTest, Test_assignment_operator)
{
    double v_data[2] = {1, 2};
    const CAROM::Vector v(v_data, 2, false, true);
    CAROM::Vector w(2, false);
    w = v;

    EXPECT_FALSE(w.distributed());
    EXPECT_EQ(w.dim(), 2);
    EXPECT_DOUBLE_EQ(w(0), 1);
    EXPECT_DOUBLE_EQ(w(1), 2);
}

/** Test methods that operate on Vector objects
 *
 *  * double norm() const
 *
 *  * double normalize()
 *
 *  * Vector* plus(const Vector&) const
 *
 *  * Vector* plus(const Vector*) const
 *
 *  * void plus(const Vector&, Vector*&) const
 *
 *  * void plus(const Vector&, Vector&) const
 *
 *  * Vector* plusAx(double, const Vector&)
 *
 *  * Vector* plusAx(double, const Vector*)
 *
 *  * void plusAx(double, const Vector&, Vector*&) const
 *
 *  * void plusAx(double, const Vector&, Vector&) const
 *
 *  * void plusEqAx(double, const Vector&)
 *
 *  * void plusEqAx(double, const Vector*)
 *
 *  * Vector* minus(const Vector&) const
 *
 *  * Vector* minus(const Vector*) const
 *
 *  * void minus(const Vector&, Vector*&) const
 *
 *  * void minus(const Vector&, Vector&) const
 *
 *  * Vector* mult(double) const
 *
 *  * void mult(double, Vector*&) const
 *
 *  * void mult(double, Vector&) const
 */

TEST(VectorSerialTest, Test_norm)
{
    CAROM::Vector v(2, false);

    /* 2-norm of (1, 1) = sqrt(2) */
    v(0) = 1;
    v(1) = 1;
    EXPECT_DOUBLE_EQ(v.norm(), M_SQRT2);

    /* 2-norm of (-1, 1) = sqrt(2) */
    v(0) = -1;
    v(1) = 1;
    EXPECT_DOUBLE_EQ(v.norm(), M_SQRT2);

    /* 2-norm of (3, 4) = 5 */
    v(0) = 3;
    v(1) = 4;
    EXPECT_DOUBLE_EQ(v.norm(), 5);

    /* 2-norm of (5, 12) = 13 */
    v(0) = 5;
    v(1) = 12;
    EXPECT_DOUBLE_EQ(v.norm(), 13);

    /** TODO(oxberry1@llnl.gov): Add more thorough tests **/
}

TEST(VectorSerialTest, Test_normalize)
{
    CAROM::Vector v(2, false);

    /* (3, 4) normalized is (.6, .8) */
    v(0) = 3;
    v(1) = 4;
    EXPECT_DOUBLE_EQ(v.normalize(), 5);
    EXPECT_DOUBLE_EQ(v(0), 0.6);
    EXPECT_DOUBLE_EQ(v(1), 0.8);
    EXPECT_DOUBLE_EQ(v.norm(), 1);
}

TEST(VectorSerialTest, Test_inner_product_const_reference)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    /* [ 1,  1]^{T} [ 1,  1] =   2 */
    EXPECT_DOUBLE_EQ(v.inner_product(v),   2);

    /* [ 1,  1]^{T} [-1,  1] =   0 */
    EXPECT_DOUBLE_EQ(v.inner_product(w),   0);

    /* [ 1,  1]^{T} [ 3,  4] =   7 */
    EXPECT_DOUBLE_EQ(v.inner_product(x),   7);

    /* [ 1,  1]^{T} [ 5, 12] =  17 */
    EXPECT_DOUBLE_EQ(v.inner_product(y),  17);

    /* [-1,  1]^{T} [ 1,  1] =   0 */
    EXPECT_DOUBLE_EQ(w.inner_product(v),   0);

    /* [-1,  1]^{T} [-1,  1] =   2 */
    EXPECT_DOUBLE_EQ(w.inner_product(w),   2);

    /* [-1,  1]^{T} [ 3,  4] =   1 */
    EXPECT_DOUBLE_EQ(w.inner_product(x),   1);

    /* [-1,  1]^{T} [ 5, 12] =   7 */
    EXPECT_DOUBLE_EQ(w.inner_product(y),   7);

    /* [ 3,  4]^{T} [ 1,  1] =   7 */
    EXPECT_DOUBLE_EQ(x.inner_product(v),   7);

    /* [ 3,  4]^{T} [-1,  1] =   1 */
    EXPECT_DOUBLE_EQ(x.inner_product(w),   1);

    /* [ 3,  4]^{T} [ 3,  4] =  25 */
    EXPECT_DOUBLE_EQ(x.inner_product(x),  25);

    /* [ 3,  4]^{T} [ 5, 12] =  63 */
    EXPECT_DOUBLE_EQ(x.inner_product(y),  63);

    /* [ 5, 12]^{T} [ 1,  1] =  17 */
    EXPECT_DOUBLE_EQ(y.inner_product(v),  17);

    /* [ 5, 12]^{T} [-1,  1] =   7 */
    EXPECT_DOUBLE_EQ(y.inner_product(w),   7);

    /* [ 5, 12]^{T} [ 3,  4] =  63 */
    EXPECT_DOUBLE_EQ(y.inner_product(x),  63);

    /* [ 5, 12]^{T} [ 5, 12] = 169 */
    EXPECT_DOUBLE_EQ(y.inner_product(y), 169);
}

TEST(VectorSerialTest, Test_inner_product_const_pointer)
{
    CAROM::Vector *v = new CAROM::Vector(2, false);
    (*v)(0) =  1;
    (*v)(1) =  1;
    CAROM::Vector *w = new CAROM::Vector(2, false);
    (*w)(0) = -1;
    (*w)(1) =  1;
    CAROM::Vector *x = new CAROM::Vector(2, false);
    (*x)(0) =  3;
    (*x)(1) =  4;
    CAROM::Vector *y = new CAROM::Vector(2, false);
    (*y)(0) =  5;
    (*y)(1) = 12;

    /* [ 1,  1]^{T} [ 1,  1] =   2 */
    EXPECT_DOUBLE_EQ(v->inner_product(v),   2);

    /* [ 1,  1]^{T} [-1,  1] =   0 */
    EXPECT_DOUBLE_EQ(v->inner_product(w),   0);

    /* [ 1,  1]^{T} [ 3,  4] =   7 */
    EXPECT_DOUBLE_EQ(v->inner_product(x),   7);

    /* [ 1,  1]^{T} [ 5, 12] =  17 */
    EXPECT_DOUBLE_EQ(v->inner_product(y),  17);

    /* [-1,  1]^{T} [ 1,  1] =   0 */
    EXPECT_DOUBLE_EQ(w->inner_product(v),   0);

    /* [-1,  1]^{T} [-1,  1] =   2 */
    EXPECT_DOUBLE_EQ(w->inner_product(w),   2);

    /* [-1,  1]^{T} [ 3,  4] =   1 */
    EXPECT_DOUBLE_EQ(w->inner_product(x),   1);

    /* [-1,  1]^{T} [ 5, 12] =   7 */
    EXPECT_DOUBLE_EQ(w->inner_product(y),   7);

    /* [ 3,  4]^{T} [ 1,  1] =   7 */
    EXPECT_DOUBLE_EQ(x->inner_product(v),   7);

    /* [ 3,  4]^{T} [-1,  1] =   1 */
    EXPECT_DOUBLE_EQ(x->inner_product(w),   1);

    /* [ 3,  4]^{T} [ 3,  4] =  25 */
    EXPECT_DOUBLE_EQ(x->inner_product(x),  25);

    /* [ 3,  4]^{T} [ 5, 12] =  63 */
    EXPECT_DOUBLE_EQ(x->inner_product(y),  63);

    /* [ 5, 12]^{T} [ 1,  1] =  17 */
    EXPECT_DOUBLE_EQ(y->inner_product(v),  17);

    /* [ 5, 12]^{T} [-1,  1] =   7 */
    EXPECT_DOUBLE_EQ(y->inner_product(w),   7);

    /* [ 5, 12]^{T} [ 3,  4] =  63 */
    EXPECT_DOUBLE_EQ(y->inner_product(x),  63);

    /* [ 5, 12]^{T} [ 5, 12] = 169 */
    EXPECT_DOUBLE_EQ(y->inner_product(y), 169);

    delete v;
    delete w;
    delete x;
    delete y;
}

TEST(VectorSerialTest, Test_plus_const_reference)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    CAROM::Vector *result;

    /* ( 1,  1) + ( 1,  1) = ( 2,  2) */
    result = v.plus(v);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  2);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + (-1,  1) = ( 0,  2) */
    result = v.plus(w);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  0);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 3,  4) = ( 4,  5) */
    result = v.plus(x);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  4);
    EXPECT_DOUBLE_EQ((*result)(1),  5);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 5, 12) = ( 6, 13) */
    result = v.plus(y);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  6);
    EXPECT_DOUBLE_EQ((*result)(1), 13);
    delete result;
    result = NULL;
}

TEST(VectorSerialTest, Test_plus_const_pointer)
{
    CAROM::Vector *v = new CAROM::Vector(2, false);
    (*v)(0) =  1;
    (*v)(1) =  1;
    CAROM::Vector *w = new CAROM::Vector(2, false);
    (*w)(0) = -1;
    (*w)(1) =  1;
    CAROM::Vector *x = new CAROM::Vector(2, false);
    (*x)(0) =  3;
    (*x)(1) =  4;
    CAROM::Vector *y = new CAROM::Vector(2, false);
    (*y)(0) =  5;
    (*y)(1) = 12;

    CAROM::Vector *result;

    /* ( 1,  1) + ( 1,  1) = ( 2,  2) */
    result = v->plus(v);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  2);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + (-1,  1) = ( 0,  2) */
    result = v->plus(w);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  0);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 3,  4) = ( 4,  5) */
    result = v->plus(x);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  4);
    EXPECT_DOUBLE_EQ((*result)(1),  5);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 5, 12) = ( 6, 13) */
    result = v->plus(y);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  6);
    EXPECT_DOUBLE_EQ((*result)(1), 13);
    delete result;
    result = NULL;

    delete v;
    delete w;
    delete x;
    delete y;
}

/* TODO(oxberry1@llnl.gov): Test cases where pointer already allocated */
TEST(VectorSerialTest, Test_plus_const_reference_pointer)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    /*
      NOTE(oxberry1@llnl.gov): if assignment omitted, pointer has
      indeterminate value that is probably non-NULL, so
      CAROM::Vector::plus tries to assign to that memory, resulting in a
      segfault.
    */
    CAROM::Vector *result = NULL;

    /* ( 1,  1) + ( 1,  1) = ( 2,  2) */
    v.plus(v, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  2);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + (-1,  1) = ( 0,  2) */
    v.plus(w, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  0);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 3,  4) = ( 4,  5) */
    v.plus(x, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  4);
    EXPECT_DOUBLE_EQ((*result)(1),  5);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 5, 12) = ( 6, 13) */
    v.plus(y, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  6);
    EXPECT_DOUBLE_EQ((*result)(1), 13);
    delete result;
    result = NULL;
}

/*
  TODO(oxberry1@llnl.gov): Test cases where output vector must be resized.
*/
TEST(VectorSerialTest, Test_plus_const_reference_reference)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    CAROM::Vector result(2, false);

    /* ( 1,  1) + ( 1,  1) = ( 2,  2) */
    v.plus(v, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),  2);
    EXPECT_DOUBLE_EQ(result(1),  2);

    /* ( 1,  1) + (-1,  1) = ( 0,  2) */
    v.plus(w, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),  0);
    EXPECT_DOUBLE_EQ(result(1),  2);

    /* ( 1,  1) + ( 3,  4) = ( 4,  5) */
    v.plus(x, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),  4);
    EXPECT_DOUBLE_EQ(result(1),  5);

    /* ( 1,  1) + ( 5, 12) = ( 6, 13) */
    v.plus(y, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),  6);
    EXPECT_DOUBLE_EQ(result(1), 13);
}

/*
   TODO(oxberry1@llnl.gov): Test with double argument set to
   something other than 1.
*/
TEST(VectorSerialTest, Test_plusAx_const_reference)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    CAROM::Vector *result;

    /* ( 1,  1) + ( 1,  1) = ( 2,  2) */
    result = v.plusAx(1.0, v);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  2);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + (-1,  1) = ( 0,  2) */
    result = v.plusAx(1.0, w);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  0);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 3,  4) = ( 4,  5) */
    result = v.plusAx(1.0, x);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  4);
    EXPECT_DOUBLE_EQ((*result)(1),  5);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 5, 12) = ( 6, 13) */
    result = v.plusAx(1.0, y);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  6);
    EXPECT_DOUBLE_EQ((*result)(1), 13);
    delete result;
    result = NULL;
}

/*
   TODO(oxberry1@llnl.gov): Test with double argument set to
   something other than 1.
*/
TEST(VectorSerialTest, Test_plusAx_const_pointer)
{
    CAROM::Vector *v = new CAROM::Vector(2, false);
    (*v)(0) =  1;
    (*v)(1) =  1;
    CAROM::Vector *w = new CAROM::Vector(2, false);
    (*w)(0) = -1;
    (*w)(1) =  1;
    CAROM::Vector *x = new CAROM::Vector(2, false);
    (*x)(0) =  3;
    (*x)(1) =  4;
    CAROM::Vector *y = new CAROM::Vector(2, false);
    (*y)(0) =  5;
    (*y)(1) = 12;

    CAROM::Vector *result;

    /* ( 1,  1) + ( 1,  1) = ( 2,  2) */
    result = v->plusAx(1.0, v);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  2);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + (-1,  1) = ( 0,  2) */
    result = v->plusAx(1.0, w);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  0);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 3,  4) = ( 4,  5) */
    result = v->plusAx(1.0, x);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  4);
    EXPECT_DOUBLE_EQ((*result)(1),  5);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 5, 12) = ( 6, 13) */
    result = v->plusAx(1.0, y);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  6);
    EXPECT_DOUBLE_EQ((*result)(1), 13);
    delete result;
    result = NULL;

    delete v;
    delete w;
    delete x;
    delete y;
}

/*
   TODO(oxberry1@llnl.gov): Test with double argument set to
   something other than 1.
*/
/* TODO(oxberry1@llnl.gov): Test cases where pointer already allocated */
TEST(VectorSerialTest, Test_plusAx_const_reference_pointer)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    /*
      NOTE(oxberry1@llnl.gov): if assignment omitted, pointer has
      indeterminate value that is probably non-NULL, so
      CAROM::Vector::plus tries to assign to that memory, resulting in a
      segfault.
    */
    CAROM::Vector *result = NULL;

    /* ( 1,  1) + ( 1,  1) = ( 2,  2) */
    v.plusAx(1.0, v, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  2);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + (-1,  1) = ( 0,  2) */
    v.plusAx(1.0, w, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  0);
    EXPECT_DOUBLE_EQ((*result)(1),  2);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 3,  4) = ( 4,  5) */
    v.plusAx(1.0, x, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  4);
    EXPECT_DOUBLE_EQ((*result)(1),  5);
    delete result;
    result = NULL;

    /* ( 1,  1) + ( 5, 12) = ( 6, 13) */
    v.plusAx(1.0, y, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  6);
    EXPECT_DOUBLE_EQ((*result)(1), 13);
    delete result;
    result = NULL;
}

/*
   TODO(oxberry1@llnl.gov): Test with double argument set to
   something other than 1.
*/
/*
  TODO(oxberry1@llnl.gov): Test cases where output vector must be resized.
*/
TEST(VectorSerialTest, Test_plusAx_const_reference_reference)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    CAROM::Vector result(2, false);

    /* ( 1,  1) + ( 1,  1) = ( 2,  2) */
    v.plusAx(1.0, v, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),  2);
    EXPECT_DOUBLE_EQ(result(1),  2);

    /* ( 1,  1) + (-1,  1) = ( 0,  2) */
    v.plusAx(1.0, w, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),  0);
    EXPECT_DOUBLE_EQ(result(1),  2);

    /* ( 1,  1) + ( 3,  4) = ( 4,  5) */
    v.plusAx(1.0, x, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),  4);
    EXPECT_DOUBLE_EQ(result(1),  5);

    /* ( 1,  1) + ( 5, 12) = ( 6, 13) */
    v.plusAx(1.0, y, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),  6);
    EXPECT_DOUBLE_EQ(result(1), 13);
}

/*
   TODO(oxberry1@llnl.gov): Test with double argument set to
   something other than 1.
*/
TEST(VectorSerialTest, Test_plusEqAx_const_reference)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    /* ( 1,  1) + ( 1,  1) = ( 2,  2) */
    v.plusEqAx(1.0, v);
    EXPECT_FALSE(v.distributed());
    EXPECT_EQ(v.dim(), 2);
    EXPECT_DOUBLE_EQ(v(0),  2);
    EXPECT_DOUBLE_EQ(v(1),  2);

    /* ( 1,  1) + (-1,  1) = ( 0,  2) */
    v(0) =  1;
    v(1) =  1;
    v.plusEqAx(1.0, w);
    EXPECT_FALSE(v.distributed());
    EXPECT_EQ(v.dim(), 2);
    EXPECT_DOUBLE_EQ(v(0),  0);
    EXPECT_DOUBLE_EQ(v(1),  2);

    /* ( 1,  1) + ( 3,  4) = ( 4,  5) */
    v(0) =  1;
    v(1) =  1;
    v.plusEqAx(1.0, x);
    EXPECT_FALSE(v.distributed());
    EXPECT_EQ(v.dim(), 2);
    EXPECT_DOUBLE_EQ(v(0),  4);
    EXPECT_DOUBLE_EQ(v(1),  5);

    /* ( 1,  1) + ( 5, 12) = ( 6, 13) */
    v(0) =  1;
    v(1) =  1;
    v.plusEqAx(1.0, y);
    EXPECT_FALSE(v.distributed());
    EXPECT_EQ(v.dim(), 2);
    EXPECT_DOUBLE_EQ(v(0),  6);
    EXPECT_DOUBLE_EQ(v(1), 13);
}

/*
   TODO(oxberry1@llnl.gov): Test with double argument set to
   something other than 1.
*/
TEST(VectorSerialTest, Test_plusEqAx_const_pointer)
{
    CAROM::Vector *v = new CAROM::Vector(2, false);
    (*v)(0) =  1;
    (*v)(1) =  1;
    CAROM::Vector *w = new CAROM::Vector(2, false);
    (*w)(0) = -1;
    (*w)(1) =  1;
    CAROM::Vector *x = new CAROM::Vector(2, false);
    (*x)(0) =  3;
    (*x)(1) =  4;
    CAROM::Vector *y = new CAROM::Vector(2, false);
    (*y)(0) =  5;
    (*y)(1) = 12;

    /* ( 1,  1) + ( 1,  1) = ( 2,  2) */
    v->plusEqAx(1.0, v);
    EXPECT_FALSE(v->distributed());
    EXPECT_EQ(v->dim(), 2);
    EXPECT_DOUBLE_EQ((*v)(0),  2);
    EXPECT_DOUBLE_EQ((*v)(1),  2);

    /* ( 1,  1) + (-1,  1) = ( 0,  2) */
    (*v)(0) =  1;
    (*v)(1) =  1;
    v->plusEqAx(1.0, w);
    EXPECT_FALSE(v->distributed());
    EXPECT_EQ(v->dim(), 2);
    EXPECT_DOUBLE_EQ((*v)(0),  0);
    EXPECT_DOUBLE_EQ((*v)(1),  2);

    /* ( 1,  1) + ( 3,  4) = ( 4,  5) */
    (*v)(0) =  1;
    (*v)(1) =  1;
    v->plusEqAx(1.0, x);
    EXPECT_FALSE(v->distributed());
    EXPECT_EQ(v->dim(), 2);
    EXPECT_DOUBLE_EQ((*v)(0),  4);
    EXPECT_DOUBLE_EQ((*v)(1),  5);

    /* ( 1,  1) + ( 5, 12) = ( 6, 13) */
    (*v)(0) =  1;
    (*v)(1) =  1;
    v->plusEqAx(1.0, y);
    EXPECT_FALSE(v->distributed());
    EXPECT_EQ(v->dim(), 2);
    EXPECT_DOUBLE_EQ((*v)(0),  6);
    EXPECT_DOUBLE_EQ((*v)(1), 13);

    delete v;
    delete w;
    delete x;
    delete y;
}

TEST(VectorSerialTest, Test_minus_const_reference)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    CAROM::Vector *result;

    /* ( 1,  1) - ( 1,  1) = ( 0,   0) */
    result = v.minus(v);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  0);
    EXPECT_DOUBLE_EQ((*result)(1),  0);
    delete result;
    result = NULL;

    /* ( 1,  1) - (-1,  1) = ( 2,   0) */
    result = v.minus(w);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  2);
    EXPECT_DOUBLE_EQ((*result)(1),  0);
    delete result;
    result = NULL;

    /* ( 1,  1) - ( 3,  4) = (-2,  -3) */
    result = v.minus(x);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0), -2);
    EXPECT_DOUBLE_EQ((*result)(1), -3);
    delete result;
    result = NULL;

    /* ( 1,  1) - ( 5, 12) = (-4, -11) */
    result = v.minus(y);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0), -4);
    EXPECT_DOUBLE_EQ((*result)(1),-11);
    delete result;
    result = NULL;
}

TEST(VectorSerialTest, Test_minus_const_pointer)
{
    CAROM::Vector *v = new CAROM::Vector(2, false);
    (*v)(0) =  1;
    (*v)(1) =  1;
    CAROM::Vector *w = new CAROM::Vector(2, false);
    (*w)(0) = -1;
    (*w)(1) =  1;
    CAROM::Vector *x = new CAROM::Vector(2, false);
    (*x)(0) =  3;
    (*x)(1) =  4;
    CAROM::Vector *y = new CAROM::Vector(2, false);
    (*y)(0) =  5;
    (*y)(1) = 12;

    CAROM::Vector *result;

    /* ( 1,  1) - ( 1,  1) = ( 0,   0) */
    result = v->minus(v);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  0);
    EXPECT_DOUBLE_EQ((*result)(1),  0);
    delete result;
    result = NULL;

    /* ( 1,  1) - (-1,  1) = ( 2,   0) */
    result = v->minus(w);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  2);
    EXPECT_DOUBLE_EQ((*result)(1),  0);
    delete result;
    result = NULL;

    /* ( 1,  1) - ( 3,  4) = (-2,  -3) */
    result = v->minus(x);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0), -2);
    EXPECT_DOUBLE_EQ((*result)(1), -3);
    delete result;
    result = NULL;

    /* ( 1,  1) - ( 5, 12) = (-4, -11) */
    result = v->minus(y);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0), -4);
    EXPECT_DOUBLE_EQ((*result)(1),-11);
    delete result;
    result = NULL;

    delete v;
    delete w;
    delete x;
    delete y;
}

/* TODO(oxberry1@llnl.gov): Test cases where pointer already allocated */
TEST(VectorSerialTest, Test_minus_const_reference_pointer)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    /*
      NOTE(oxberry1@llnl.gov): if assignment omitted, pointer has
      indeterminate value that is probably non-NULL, so
      CAROM::Vector::minus tries to assign to that memory, resulting in a
      segfault.
    */
    CAROM::Vector *result = NULL;

    /* ( 1,  1) - ( 1,  1) = ( 0,   0) */
    v.minus(v, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  0);
    EXPECT_DOUBLE_EQ((*result)(1),  0);
    delete result;
    result = NULL;

    /* ( 1,  1) - (-1,  1) = ( 2,   0) */
    v.minus(w, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),  2);
    EXPECT_DOUBLE_EQ((*result)(1),  0);
    delete result;
    result = NULL;

    /* ( 1,  1) - ( 3,  4) = (-2,  -3) */
    v.minus(x, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0), -2);
    EXPECT_DOUBLE_EQ((*result)(1), -3);
    delete result;
    result = NULL;

    /* ( 1,  1) - ( 5, 12) = (-4, -11) */
    v.minus(y, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0), -4);
    EXPECT_DOUBLE_EQ((*result)(1),-11);
    delete result;
    result = NULL;
}

/*
  TODO(oxberry1@llnl.gov): Test cases where output vector must be resized.
*/
TEST(VectorSerialTest, Test_minus_const_reference_reference)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    CAROM::Vector result(2, false);

    /* ( 1,  1) - ( 1,  1) = ( 0,   0) */
    v.minus(v, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),  0);
    EXPECT_DOUBLE_EQ(result(1),  0);

    /* ( 1,  1) - (-1,  1) = ( 2,   0) */
    v.minus(w, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),  2);
    EXPECT_DOUBLE_EQ(result(1),  0);

    /* ( 1,  1) - ( 3,  4) = (-2,  -3) */
    v.minus(x, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0), -2);
    EXPECT_DOUBLE_EQ(result(1), -3);

    /* ( 1,  1) - ( 5, 12) = (-4, -11) */
    v.minus(y, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0), -4);
    EXPECT_DOUBLE_EQ(result(1),-11);
}

TEST(VectorSerialTest, Test_mult_double)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    CAROM::Vector *result;

    result = v.mult(2);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),   2);
    EXPECT_DOUBLE_EQ((*result)(1),   2);
    delete result;
    result = NULL;

    result = w.mult(-5);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),   5);
    EXPECT_DOUBLE_EQ((*result)(1),  -5);
    delete result;
    result = NULL;

    result = x.mult(3);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),   9);
    EXPECT_DOUBLE_EQ((*result)(1),  12);
    delete result;
    result = NULL;

    result = y.mult(0.5);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0), 2.5);
    EXPECT_DOUBLE_EQ((*result)(1),   6);
    delete result;
    result = NULL;
}

TEST(VectorSerialTest, Test_mult_double_pointer)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    CAROM::Vector *result = NULL;

    v.mult(2, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),   2);
    EXPECT_DOUBLE_EQ((*result)(1),   2);
    delete result;
    result = NULL;

    w.mult(-5, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),   5);
    EXPECT_DOUBLE_EQ((*result)(1),  -5);
    delete result;
    result = NULL;

    x.mult(3, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0),   9);
    EXPECT_DOUBLE_EQ((*result)(1),  12);
    delete result;
    result = NULL;

    y.mult(0.5, result);
    EXPECT_FALSE(result->distributed());
    EXPECT_EQ(result->dim(), 2);
    EXPECT_DOUBLE_EQ((*result)(0), 2.5);
    EXPECT_DOUBLE_EQ((*result)(1),   6);
    delete result;
    result = NULL;
}

TEST(VectorSerialTest, Test_mult_double_reference)
{
    CAROM::Vector v(2, false);
    v(0) =  1;
    v(1) =  1;
    CAROM::Vector w(2, false);
    w(0) = -1;
    w(1) =  1;
    CAROM::Vector x(2, false);
    x(0) =  3;
    x(1) =  4;
    CAROM::Vector y(2, false);
    y(0) =  5;
    y(1) = 12;

    CAROM::Vector result(2, false);

    v.mult(2, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),   2);
    EXPECT_DOUBLE_EQ(result(1),   2);

    w.mult(-5, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),   5);
    EXPECT_DOUBLE_EQ(result(1),  -5);

    x.mult(3, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0),   9);
    EXPECT_DOUBLE_EQ(result(1),  12);

    y.mult(0.5, result);
    EXPECT_FALSE(result.distributed());
    EXPECT_EQ(result.dim(), 2);
    EXPECT_DOUBLE_EQ(result(0), 2.5);
    EXPECT_DOUBLE_EQ(result(1),   6);
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#else // #ifndef CAROM_HAS_GTEST
int main()
{
    std::cout << "libROM was compiled without Google Test support, so unit "
              << "tests have been disabled. To enable unit tests, compile "
              << "libROM with Google Test support." << std::endl;
}
#endif // #endif CAROM_HAS_GTEST
