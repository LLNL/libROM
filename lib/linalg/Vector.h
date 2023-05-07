/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A simple, parallel Vector class with the utility needed to
//              support the basis generation methods of this library.  A
//              distributed Vector has its rows distributed across processors.

#ifndef included_Vector_h
#define included_Vector_h

#include "utils/Utilities.h"
#include <vector>
#include <functional>

namespace CAROM {

/**
 * Class Vector is a simple vector class in which the dimensions may be
 * distributed across multiple processes.  This class supports only the basic
 * operations that are needed by the SVD library.
 */
class Vector
{
public:
    Vector();

    /**
     * @brief Constructor creating a Vector with uninitialized values.
     *
     * @pre dim > 0
     *
     * @param[in] dim When undistributed, the total dimension of the Vector.
     *                When distributed, the part of the total dimension of
     *                the Vector on this processor.
     * @param[in] distributed If true the dimensions of the Vector are spread
     *                        over all processors.
     */
    Vector(
        int dim,
        bool distributed);

    /**
     * @brief Constructor in which the Vector is given its initial values.
     *
     * @pre vec != 0
     * @pre dim > 0
     *
     * @param[in] vec The initial values of the Vector.
     * @param[in] dim When undistributed, the total dimension of the Vector.
     *                When distributed, the part of the total dimension of
     *                the Vector on this processor.
     * @param[in] distributed If true the dimensions of the Vector are spread
     *                        over all processors.
     * @param[in] copy_data If true the vector allocates its own storage and
     *                      copies the contents of vec into its own storage.
     *                      Otherwise it uses vec as its storage.
     */
    Vector(
        double* vec,
        int dim,
        bool distributed,
        bool copy_data = true);

    /**
     * @brief Copy constructor.
     *
     * @param[in] other The Vector to copy.
     */
    Vector(
        const Vector& other);

    /**
     * @brief Destructor.
     */
    ~Vector();

    /**
     * @brief Assignment operator.
     *
     * @param[in] rhs The Vector to assign to this.
     *
     * @return This after rhs has been assigned to it.
     */
    Vector&
    operator = (
        const Vector& rhs);

    /**
     * @brief Addition operator.
     *
     * @param[in] rhs The Vector to add to this.
     *
     * @return This after rhs has been added to it.
     */
    Vector&
    operator += (
        const Vector& rhs);

    /**
     * @brief Subtraction operator.
     *
     * @param[in] rhs The Vector to subtract from this.
     *
     * @return This after rhs has been subtracted from it.
     */
    Vector&
    operator -= (
        const Vector& rhs);

    /**
     * @brief Equal operator.
     *
     * @param[in] a The double precision number to which every
     * Vector entry should be set.
     *
     * @return This with every element of the Vector set to a.
     */
    Vector&
    operator = (
        const double& a);

    /**
     * @brief Scaling operator.
     *
     * @param[in] a The double precision number by which every
     * Vector entry should be scaled.
     *
     * @return This with every element of the Vector scaled by a.
     */
    Vector&
    operator *= (
        const double& a);

    /**
     * @brief Transform the vector using a supplied function.
     *
     * @param[in] transformer A transformer function which takes in as input a
     *            size and a vector.
     *
     * @return The newly transformed vector.
     */
    Vector&
    transform(std::function<void(const int size, double* vector)> transformer);

    /**
     * @brief Transform a vector using a supplied function and store the
     *        results in another vector.
     *
     * @param[out] result A vector which will store the transformed result.
     *
     * @param[in] transformer A transformer function which takes in as input a
    *             size and transforms the vector.
     */
    void
    transform(Vector& result,
              std::function<void(const int size, double* vector)> transformer) const;

    /**
     * @brief Transform a vector using a supplied function and store the
     *        results in another vector.
     *
     * @param[out] result A vector which will store the transformed result.
     *
     * @param[in] transformer A transformer function which takes in as input a
     *            size and transforms the vector.
     */
    void
    transform(Vector*& result,
              std::function<void(const int size, double* vector)> transformer) const;

    /**
     * @brief Transform the vector using a supplied function.
     *
     * @param[in] transformer A transformer function which takes in as input a
     *            size and transforms the origVector and stores the result in
     *            resultVector.
     *
     * @return The newly transformed vector.
     */
    Vector&
    transform(
        std::function<void(const int size, double* origVector, double* resultVector)>
        transformer);

    /**
     * @brief Transform a vector using a supplied function and store the
     *        results in another vector.
     *
     * @param[out] result A vector which will store the transformed result.
     *
     * @param[in] transformer A transformer function which takes in as input a
     *            size and transforms the origVector and stores the result in
     *            resultVector.
     */
    void
    transform(Vector& result,
              std::function<void(const int size, double* origVector, double* resultVector)>
              transformer) const;

    /**
     * @brief Transform a vector using a supplied function and store the
     *        results in another vector.
     *
     * @param[out] result A vector which will store the transformed result.
     *
     * @param[in] transformer A transformer function which takes in as input a
     *            size and transforms the origVector and stores the result in
     *            resultVector.
     */
    void
    transform(Vector*& result,
              std::function<void(const int size, double* origVector, double* resultVector)>
              transformer) const;

    /**
     * @brief Sets the length of the vector and reallocates storage if
     * needed. All values are initialized to zero.
     *
     * @param[in] dim When undistributed, the total dimension of the Vector.
     *                When distributed, the part of the total dimension of
     *                the Vector on this processor.
     */
    void
    setSize(
        int dim)
    {
        if (dim > d_alloc_size) {
            if (!d_owns_data) {
                CAROM_ERROR("Can not reallocate externally owned storage.");
            }
            if (d_vec) {
                delete [] d_vec;
            }

            // Allocate new array and initialize all values to zero.
            d_vec = new double [dim] {0.0};
            d_alloc_size = dim;
        }
        d_dim = dim;
    }

    /**
     * @brief Returns true if the Vector is distributed.
     *
     * @return True if the Vector is distributed.
     */
    bool
    distributed() const
    {
        return d_distributed;
    }

    /**
     * @brief Returns the dimension of the Vector on this processor.
     *
     * @return The part of the Vector's dimension on this processor.
     */
    int
    dim() const
    {
        return d_dim;
    }

    /**
     * @brief Inner product, reference form.
     *
     * For distributed Vectors this is a parallel operation.
     *
     * @pre dim() == other.dim()
     * @pre distributed() == other.distributed()
     *
     * @param[in] other The Vector to form the inner product with this.
     *
     * @return The inner product of this and other.
     */
    double
    inner_product(
        const Vector& other) const;

    /**
     * @brief Inner product, pointer version.
     *
     * For distributed Vectors this is a parallel operation.
     *
     * @pre other != 0
     * @pre dim() == other->dim()
     * @pre distributed() == other->distributed()
     *
     * @param[in] other The Vector to form the inner product with this.
     *
     * @return The inner product of this and other.
     */
    double
    inner_product(
        const Vector* other) const
    {
        CAROM_VERIFY(other != 0);
        return inner_product(*other);
    }

    /**
     * @brief Form the norm of this.
     *
     * For a distributed Vector this is a parallel operation.
     *
     * @return The norm of this.
     */
    double
    norm() const;

    /**
     * @brief Form the squared norm of this.
     *
     * For a distributed Vector this is a parallel operation.
     *
     * @return The squared norm of this.
     */
    double
    norm2() const;

    /**
     * @brief Normalizes the Vector and returns its norm.
     *
     * For a distributed Vector this is a parallel operation.
     *
     * @return The norm of this.
     */
    double
    normalize();

    /**
     * @brief Adds other and this and returns the result, reference version.
     *
     * @pre distributed() == other.distributed()
     * @pre dim() == other.dim()
     *
     * @param[in] other The other summand.
     *
     * @return this + other
     */
    Vector*
    plus(
        const Vector& other) const
    {
        Vector* result = 0;
        plus(other, result);
        return result;
    }

    /**
     * @brief Adds other and this and returns the result, pointer version.
     *
     * @pre other != 0
     * @pre distributed() == other->distributed()
     * @pre dim() == other->dim()
     *
     * @param[in] other The other summand.
     *
     * @return this + other
     */
    Vector*
    plus(
        const Vector* other) const
    {
        CAROM_VERIFY(other != 0);
        return plus(*other);
    }

    /**
     * @brief Adds other and this and fills result with the answer.
     *
     * Result will be allocated if unallocated or resized appropriately if
     * already allocated.
     *
     * @pre result == 0 || result->distributed() == distributed()
     * @pre distributed() == other.distributed()
     * @pre dim() == other.dim()
     *
     * @param[in] other The other summand.
     * @param[out] result this + other
     */
    void
    plus(
        const Vector& other,
        Vector*& result) const;

    /**
     * @brief Adds other and this and fills result with the answer.
     *
     * Result will be resized appropriately.
     *
     * @pre result.distributed() == distributed()
     * @pre distributed() == other.distributed()
     * @pre dim() == other.dim()
     *
     * @param[in] other The other summand.
     * @param[out] result this + other
     */
    void
    plus(
        const Vector& other,
        Vector& result) const;

    /**
     * @brief Adds factor*other and this and returns the result, reference
     * version.
     *
     * @pre distributed() == other.distributed()
     * @pre dim() == other.dim()
     *
     * @param[in] factor Multiplicative factor applied to other.
     * @param[in] other The other summand.
     *
     * @return this + factor*other
     */
    Vector*
    plusAx(
        double factor,
        const Vector& other)
    {
        Vector* result = 0;
        plusAx(factor, other, result);
        return result;
    }

    /**
     * @brief Adds factor*other and this and returns the result, pointer
     * version.
     *
     * @pre distributed() == other->distributed()
     * @pre dim() == other->dim()
     *
     * @param[in] factor Multiplicative factor applied to other.
     * @param[in] other The other summand.
     *
     * @return this + factor*other
     */
    Vector*
    plusAx(
        double factor,
        const Vector* other)
    {
        CAROM_VERIFY(other != 0);
        return plusAx(factor, *other);
    }

    /**
     * @brief Adds factor*other and this and fills result with the answer.
     *
     * Result will be allocated if unallocated or resized appropriately if
     * already allocated.
     *
     * @pre result == 0 || result->distributed() == distributed()
     * @pre distributed() == other.distributed()
     * @pre dim() == other.dim()
     *
     * @param[in] factor Multiplicative factor applied to other.
     * @param[in] other The other summand.
     * @param[out] result this + factor*other
     */
    void
    plusAx(
        double factor,
        const Vector& other,
        Vector*& result) const;

    /**
     * @brief Adds factor*other and this and fills result with the answer.
     *
     * Result will be resized appropriately.
     *
     * @pre result.distributed() == distributed()
     * @pre distributed() == other.distributed()
     * @pre dim() == other.dim()
     *
     * @param[in] factor Multiplicative factor applied to other.
     * @param[in] other The other summand.
     * @param[out] result this + factor*other
     */
    void
    plusAx(
        double factor,
        const Vector& other,
        Vector& result) const;

    /**
     * @brief Adds factor*other to this, reference version.
     *
     * @pre distributed() == other.distributed()
     * @pre dim() == other.dim()
     *
     * @param[in] factor Multiplicative factor applied to other.
     * @param[in] other The other summand.
     */
    void
    plusEqAx(
        double factor,
        const Vector& other);

    /**
     * @brief Adds factor*other to this, pointer version.
     *
     * @pre other != 0
     * @pre distributed() == other->distributed()
     * @pre dim() == other->dim()
     *
     * @param[in] factor Multiplicative factor applied to other.
     * @param[in] other The other summand.
     */
    void
    plusEqAx(
        double factor,
        const Vector* other)
    {
        CAROM_VERIFY(other != 0);
        plusEqAx(factor, *other);
    }

    /**
     * @brief Subtracts other and this and returns the result, reference
     * version.
     *
     * @pre distributed() == other.distributed()
     * @pre dim() == other.dim()
     *
     * @param[in] other The other subtrahand.
     *
     * @return this - other
     */
    Vector*
    minus(
        const Vector& other) const
    {
        Vector* result = 0;
        minus(other, result);
        return result;
    }

    /**
     * @brief Subtracts other and this and returns the result, pointer
     * version.
     *
     * @pre other != 0
     * @pre distributed() == other->distributed()
     * @pre dim() == other->dim()
     *
     * @param[in] other The other subtrahand.
     *
     * @return this - other
     */
    Vector*
    minus(
        const Vector* other) const
    {
        CAROM_VERIFY(other != 0);
        return minus(*other);
    }

    /**
     * @brief Subtracts other and this and fills result with the answer.
     *
     * Result will be allocated if unallocated or resized appropriately if
     * already allocated.
     *
     * @pre result == 0 || result->distributed() == distributed()
     * @pre distributed() == other.distributed()
     * @pre dim() == other.dim()
     *
     * @param[in] other The other subtrahend.
     * @param[out] result this - other
     */
    void
    minus(
        const Vector& other,
        Vector*& result) const;

    /**
     * @brief Subtracts other and this and fills result with the answer.
     *
     * Result will be resized appropriately.
     *
     * @pre result.distributed() == distributed()
     * @pre distributed() == other.distributed()
     * @pre dim() == other.dim()
     *
     * @param[in] other The other subtrahend.
     * @param[out] result this - other
     */
    void
    minus(
        const Vector& other,
        Vector& result) const;

    /**
     * @brief Multiplies this by the supplied constant and returns the
     * result.
     *
     * @param[in] factor Factor to multiply by.
     *
     * @return factor*this
     */
    Vector*
    mult(
        double factor) const
    {
        Vector* result = 0;
        mult(factor, result);
        return result;
    }

    /**
     * @brief Multiplies this by the supplied constant and fills result with
     * the answer.
     *
     * @pre result == 0 || result->distributed() == distributed()
     *
     * @param[in] factor Factor to multiply by.
     * @param[out] result factor*this
     */
    void
    mult(
        double factor,
        Vector*& result) const;

    /**
     * @brief Multiplies this by the supplied constant and fills result with
     * the answer.
     *
     * @pre result.distributed() == distributed()
     *
     * @param[in] factor Factor to multiply by.
     * @param[out] result factor*this
     */
    void
    mult(
        double factor,
        Vector& result) const;

    /**
     * @brief Const Vector member access.
     *
     * @pre (0 <= i) && (i < dim())
     *
     * @param[in] i The component of the Vector on this processor requested.
     *
     * @return The requested component of the Vector on this processor.
     */
    const double&
    item(
        int i) const
    {
        CAROM_ASSERT((0 <= i) && (i < dim()));
        return d_vec[i];
    }

    /**
     * @brief Non-const Vector member access.
     *
     * Allows constructs of the form vec[i] = val;
     *
     * @pre (0 <= i) && (i < dim())
     *
     * @param[in] i The component of the Vector on this processor requested.
     *
     * @return The requested component of the Vector on this processor.
     */
    double&
    item(
        int i)
    {
        CAROM_ASSERT((0 <= i) && (i < dim()));
        return d_vec[i];
    }

    /**
     * @brief Const Vector member access.
     *
     * @pre (0 <= i) && (i < dim())
     *
     * @param[in] i The component of the Vector on this processor requested.
     *
     * @return The requested component of the Vector on this processor.
     */
    const double& operator() (int i) const
    {
        return item(i);
    }

    /**
     * @brief Non-const Vector member access.
     *
     * Allows constructs of the form vec[i] = val;
     *
     * @pre (0 <= i) && (i < dim())
     *
     * @param[in] i The component of the Vector on this processor requested.
     *
     * @return The requested component of the Vector on this processor.
     */
    double& operator() (int i)
    {
        return item(i);
    }

    /**
     * @brief print Vector into (a) ascii file(s).
     *
     * @param[in] prefix The name of the prefix of the file name.
     *
     */
    void print(const char * prefix) const;

    /**
     * @brief write Vector into (a) HDF file(s).
     *
     * @param[in] base_file_name The base part of the file name.
     *
     */
    void write(const std::string& base_file_name);

    /**
     * @brief read Vector from (a) HDF file(s).
     *
     * @param[in] base_file_name The base part of the file name.
     *
     */
    void read(const std::string& base_file_name);

    /**
     * @brief read read a single rank of a distributed Vector from (a) HDF file(s).
     *
     * @param[in] base_file_name The base part of the file name.
     * @param[in] rank The rank to read from.
     *
     */
    void local_read(const std::string& base_file_name, int rank);

    /**
     * @brief Get the vector data as a pointer.
     */
    double *getData() const
    {
        return d_vec;
    }

    /**
     * @brief Compute the local minimum of this.
     *
     * @param[in] nmax If positive, use only the first nmax entries of this.
     *
     * @return The local minimum of this.
     */
    double localMin(int nmax = 0);

private:
    /**
     * @brief The storage for the Vector's values on this processor.
     */
    double* d_vec;

    /**
     * @brief The part of the Vector's dimension on this processor.
     */
    int d_dim;

    /**
     * @brief The currently allocated size.
     *
     * d_dim <= d_alloc_size
     */
    int d_alloc_size;

    /**
     * @brief If true, the Vector's dimensions are distributed over all
     * processors.
     *
     * Each processor does not need to hold the same number of dimensions.
     */
    bool d_distributed;

    /**
     * @brief The number of processors being run on.
     */
    int d_num_procs;

    /**
     * @brief If true, this object owns its underlying data, d_vec, and
     * is responsible for its deletion.
     *
     * If d_owns_data is false, then the object may not reallocate d_vec.
     */
    bool d_owns_data;
};

/**
 * @brief Get center point of a group of points.

 */
int getCenterPoint(std::vector<Vector*>& points,
                   bool use_centroid);

/**
* @brief Get center point of a group of points.

*/
int getCenterPoint(std::vector<Vector>& points,
                   bool use_centroid);

/**
* @brief Get closest point to a test point among a group of points.

*/
int getClosestPoint(std::vector<Vector*>& points,
                    Vector* test_point);

/**
* @brief Get closest point to a test point among a group of points.

*/
int getClosestPoint(std::vector<Vector>& points,
                    Vector test_point);
}

#endif
