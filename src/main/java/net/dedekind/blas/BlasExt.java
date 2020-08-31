/*
 * Copyright 2020 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package net.dedekind.blas;

import java.util.Objects;

import net.dedekind.Order;
import net.frobenius.TTrans;

/**
 * {@code Intel MKL} (Math Kernel Library) native implementations of selected
 * BLAS extensions.
 * <p>
 * Matrix storage layout must be column-major as in Fortran. The number of
 * matrix rows and columns must always be strictly positive. All operations
 * throw a {@code NullPointerException} if any of the reference method arguments
 * is {@code null}.
 */
public class BlasExt {

    // The ordering is column-major
    private static final byte ORDERING = (byte) 'C';

    private static final boolean USE_CRITICAL = true;

    private static final BlasExt mkl;

    static {
        Blas.getInstance(true);
        mkl = new BlasExt();
    }

    public static BlasExt getInstance() {
        return mkl;
    }

    public final void cgemm3m(Trans transa, Trans transb, int m, int n, int k, float alphar, float alphai, float[] a,
            int lda, float[] b, int ldb, float betar, float betai, float[] c, int ldc) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(c, "c");
        Objects.requireNonNull(transa, "transa");
        Objects.requireNonNull(transb, "transb");
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(k, "k");
        cgemm3m_n(Order.COL.code(), transa.code(), transb.code(), m, n, k, alphar, alphai, a, lda, b, ldb, betar, betai,
                c, ldc, USE_CRITICAL);
    }

    public final void zgemm3m(Trans transa, Trans transb, int m, int n, int k, double alphar, double alphai, double[] a,
            int lda, double[] b, int ldb, double betar, double betai, double[] c, int ldc) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(c, "c");
        Objects.requireNonNull(transa, "transa");
        Objects.requireNonNull(transb, "transb");
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(k, "k");
        zgemm3m_n(Order.COL.code(), transa.code(), transb.code(), m, n, k, alphar, alphai, a, lda, b, ldb, betar, betai,
                c, ldc, USE_CRITICAL);
    }

    /**
     * Scaling and in-place transposition / copying of a float matrix
     * {@code AB := alpha *op( AB )} where the transposition operation
     * {@code op()} can be a normal matrix copy, a transposition, a conjugate
     * transposition, or just a conjugation.
     * 
     * @param trans
     *            specifies the operation type
     * @param rows
     *            the number of rows in matrix {@code AB} before the transpose
     *            operation
     * @param cols
     *            the number of columns in matrix {@code AB} before the
     *            transpose operation
     * @param alpha
     *            this parameter scales the input matrix by {@code alpha}
     * @param AB
     *            the input matrix modified in-place
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix; measured in the number of elements (this
     *            parameter must be at least {@code rows})
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            destination matrix; measured in the number of elements. If
     *            {@code trans} is {@link TTrans#TRANS} or
     *            {@link TTrans#CONJ_TRANS}, this parameter must be at least
     *            {@code max(1, cols)}; otherwise it must be at least
     *            {@code max(1, rows)}
     */
    public final void simatcopy(TTrans trans, int rows, int cols, float alpha, float[] AB, int lda, int ldb) {
        Objects.requireNonNull(trans, "trans");
        Objects.requireNonNull(AB, "AB");
        checkStrictlyPositive(rows, "rows");
        checkStrictlyPositive(cols, "cols");
        simatcopy_n(ORDERING, trans(trans), rows, cols, alpha, AB, lda, ldb, USE_CRITICAL);
    }

    /**
     * Scaling and in-place transposition / copying of a double matrix
     * {@code AB := alpha *op( AB )} where the transposition operation
     * {@code op()} can be a normal matrix copy, a transposition, a conjugate
     * transposition, or just a conjugation.
     * 
     * @param trans
     *            specifies the operation type
     * @param rows
     *            the number of rows in matrix {@code AB} before the transpose
     *            operation
     * @param cols
     *            the number of columns in matrix {@code AB} before the
     *            transpose operation
     * @param alpha
     *            this parameter scales the input matrix by {@code alpha}
     * @param AB
     *            the input matrix modified in-place
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix; measured in the number of elements (this
     *            parameter must be at least {@code rows})
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            destination matrix; measured in the number of elements. If
     *            {@code trans} is {@link TTrans#TRANS} or
     *            {@link TTrans#CONJ_TRANS}, this parameter must be at least
     *            {@code max(1, cols)}; otherwise it must be at least
     *            {@code max(1, rows)}
     */
    public final void dimatcopy(TTrans trans, int rows, int cols, double alpha, double[] AB, int lda, int ldb) {
        Objects.requireNonNull(trans, "trans");
        Objects.requireNonNull(AB, "AB");
        checkStrictlyPositive(rows, "rows");
        checkStrictlyPositive(cols, "cols");
        dimatcopy_n(ORDERING, trans(trans), rows, cols, alpha, AB, lda, ldb, USE_CRITICAL);
    }

    /**
     * Scaling and in-place transposition / copying of a single precision
     * complex matrix {@code AB := alpha *op( AB )} where the transposition
     * operation {@code op()} can be a normal matrix copy, a transposition, a
     * conjugate transposition, or just a conjugation.
     * 
     * @param trans
     *            specifies the operation type
     * @param rows
     *            the number of rows in matrix {@code AB} before the transpose
     *            operation
     * @param cols
     *            the number of columns in matrix {@code AB} before the
     *            transpose operation
     * @param alphar
     *            the real part of the scale factor for the input matrix
     * @param alphai
     *            the imaginary part of the scale factor for the input matrix
     * @param AB
     *            the input matrix modified in-place
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix; measured in the number of elements (this
     *            parameter must be at least {@code rows})
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            destination matrix; measured in the number of elements. If
     *            {@code trans} is {@link TTrans#TRANS} or
     *            {@link TTrans#CONJ_TRANS}, this parameter must be at least
     *            {@code max(1, cols)}; otherwise it must be at least
     *            {@code max(1, rows)}
     */
    public final void cimatcopy(TTrans trans, int rows, int cols, float alphar, float alphai, float[] AB, int lda,
            int ldb) {
        Objects.requireNonNull(trans, "trans");
        Objects.requireNonNull(AB, "AB");
        checkStrictlyPositive(rows, "rows");
        checkStrictlyPositive(cols, "cols");
        cimatcopy_n(ORDERING, trans(trans), rows, cols, alphar, alphai, AB, lda, ldb, USE_CRITICAL);
    }

    /**
     * Scaling and in-place transposition / copying of a double precision
     * complex matrix {@code AB := alpha *op( AB )} where the transposition
     * operation {@code op()} can be a normal matrix copy, a transposition, a
     * conjugate transposition, or just a conjugation.
     * 
     * @param trans
     *            specifies the operation type
     * @param rows
     *            the number of rows in matrix {@code AB} before the transpose
     *            operation
     * @param cols
     *            the number of columns in matrix {@code AB} before the
     *            transpose operation
     * @param alphar
     *            the real part of the scale factor for the input matrix
     * @param alphai
     *            the imaginary part of the scale factor for the input matrix
     * @param AB
     *            the input matrix modified in-place
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix; measured in the number of elements (this
     *            parameter must be at least {@code rows})
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            destination matrix; measured in the number of elements. If
     *            {@code trans} is {@link TTrans#TRANS} or
     *            {@link TTrans#CONJ_TRANS}, this parameter must be at least
     *            {@code max(1, cols)}; otherwise it must be at least
     *            {@code max(1, rows)}
     */
    public final void zimatcopy(TTrans trans, int rows, int cols, double alphar, double alphai, double[] AB, int lda,
            int ldb) {
        Objects.requireNonNull(trans, "trans");
        Objects.requireNonNull(AB, "AB");
        checkStrictlyPositive(rows, "rows");
        checkStrictlyPositive(cols, "cols");
        zimatcopy_n(ORDERING, trans(trans), rows, cols, alphar, alphai, AB, lda, ldb, USE_CRITICAL);
    }

    /**
     * Scaling and out-of-place transposition / copying of a float matrix
     * {@code B := alpha *op( A )} where the transposition operation
     * {@code op()} can be a normal matrix copy, a transposition, a conjugate
     * transposition, or just a conjugation.
     * 
     * @param trans
     *            specifies the operation type
     * @param rows
     *            the number of rows in matrix {@code A} (the source matrix)
     * @param cols
     *            the number of columns in matrix {@code A} (the source matrix)
     * @param alpha
     *            this parameter scales the input matrix by {@code alpha}
     * @param A
     *            the input matrix, unmodified
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix ({@code A}); measured in the number of elements
     *            (this parameter must be at least {@code rows})
     * @param B
     *            the output matrix, modified
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            destination matrix ({@code B}); measured in the number of
     *            elements. If {@code trans} is {@link TTrans#TRANS} or
     *            {@link TTrans#CONJ_TRANS}, this parameter must be at least
     *            {@code max(1, cols)}; otherwise it must be at least
     *            {@code max(1, rows)}
     */
    public final void somatcopy(TTrans trans, int rows, int cols, float alpha, float[] A, int lda, float[] B, int ldb) {
        Objects.requireNonNull(trans, "trans");
        Objects.requireNonNull(A, "A");
        Objects.requireNonNull(B, "B");
        checkStrictlyPositive(rows, "rows");
        checkStrictlyPositive(cols, "cols");
        somatcopy_n(ORDERING, trans(trans), rows, cols, alpha, A, lda, B, ldb, USE_CRITICAL);
    }

    /**
     * Scaling and out-of-place transposition / copying of a double matrix
     * {@code B := alpha *op( A )} where the transposition operation
     * {@code op()} can be a normal matrix copy, a transposition, a conjugate
     * transposition, or just a conjugation.
     * 
     * @param trans
     *            specifies the operation type
     * @param rows
     *            the number of rows in matrix {@code A} (the source matrix)
     * @param cols
     *            the number of columns in matrix {@code A} (the source matrix)
     * @param alpha
     *            this parameter scales the input matrix by {@code alpha}
     * @param A
     *            the input matrix, unmodified
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix ({@code A}); measured in the number of elements
     *            (this parameter must be at least {@code rows})
     * @param B
     *            the output matrix, modified
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            destination matrix ({@code B}); measured in the number of
     *            elements. If {@code trans} is {@link TTrans#TRANS} or
     *            {@link TTrans#CONJ_TRANS}, this parameter must be at least
     *            {@code max(1, cols)}; otherwise it must be at least
     *            {@code max(1, rows)}
     */
    public final void domatcopy(TTrans trans, int rows, int cols, double alpha, double[] A, int lda, double[] B,
            int ldb) {
        Objects.requireNonNull(trans, "trans");
        Objects.requireNonNull(A, "A");
        Objects.requireNonNull(B, "B");
        checkStrictlyPositive(rows, "rows");
        checkStrictlyPositive(cols, "cols");
        domatcopy_n(ORDERING, trans(trans), rows, cols, alpha, A, lda, B, ldb, USE_CRITICAL);
    }

    /**
     * Scaling and out-of-place transposition / copying of a single precision
     * complex matrix {@code B := alpha *op( A )} where the transposition
     * operation {@code op()} can be a normal matrix copy, a transposition, a
     * conjugate transposition, or just a conjugation.
     * 
     * @param trans
     *            specifies the operation type
     * @param rows
     *            the number of rows in matrix {@code A} (the source matrix)
     * @param cols
     *            the number of columns in matrix {@code A} (the source matrix)
     * @param alphar
     *            the real part of the scale factor for the input matrix
     * @param alphai
     *            the imaginary part of the scale factor for the input matrix
     * @param A
     *            the input matrix, unmodified
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix ({@code A}); measured in the number of elements
     *            (this parameter must be at least {@code rows})
     * @param B
     *            the output matrix, modified
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            destination matrix ({@code B}); measured in the number of
     *            elements. If {@code trans} is {@link TTrans#TRANS} or
     *            {@link TTrans#CONJ_TRANS}, this parameter must be at least
     *            {@code max(1, cols)}; otherwise it must be at least
     *            {@code max(1, rows)}
     */
    public final void comatcopy(TTrans trans, int rows, int cols, float alphar, float alphai, float[] A, int lda,
            float[] B, int ldb) {
        Objects.requireNonNull(trans, "trans");
        Objects.requireNonNull(A, "A");
        Objects.requireNonNull(B, "B");
        checkStrictlyPositive(rows, "rows");
        checkStrictlyPositive(cols, "cols");
        comatcopy_n(ORDERING, trans(trans), rows, cols, alphar, alphai, A, lda, B, ldb, USE_CRITICAL);
    }

    /**
     * Scaling and out-of-place transposition / copying of a double precision
     * complex matrix {@code B := alpha *op( A )} where the transposition
     * operation {@code op()} can be a normal matrix copy, a transposition, a
     * conjugate transposition, or just a conjugation.
     * 
     * @param trans
     *            specifies the operation type
     * @param rows
     *            the number of rows in matrix {@code A} (the source matrix)
     * @param cols
     *            the number of columns in matrix {@code A} (the source matrix)
     * @param alphar
     *            the real part of the scale factor for the input matrix
     * @param alphai
     *            the imaginary part of the scale factor for the input matrix
     * @param A
     *            the input matrix, unmodified
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix ({@code A}); measured in the number of elements
     *            (this parameter must be at least {@code rows})
     * @param B
     *            the output matrix, modified
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            destination matrix ({@code B}); measured in the number of
     *            elements. If {@code trans} is {@link TTrans#TRANS} or
     *            {@link TTrans#CONJ_TRANS}, this parameter must be at least
     *            {@code max(1, cols)}; otherwise it must be at least
     *            {@code max(1, rows)}
     */
    public final void zomatcopy(TTrans trans, int rows, int cols, double alphar, double alphai, double[] A, int lda,
            double[] B, int ldb) {
        Objects.requireNonNull(trans, "trans");
        Objects.requireNonNull(A, "A");
        Objects.requireNonNull(B, "B");
        checkStrictlyPositive(rows, "rows");
        checkStrictlyPositive(cols, "cols");
        zomatcopy_n(ORDERING, trans(trans), rows, cols, alphar, alphai, A, lda, B, ldb, USE_CRITICAL);
    }

    /**
     * Scales and adds two float matrices, as well as performing out-of-place
     * transposition operations {@code C := alpha *op(A) + beta *op(B)} where
     * the {@code op()} operations are transpose, conjugate-transpose, conjugate
     * (no transpose), or no transpose, depending on the values of
     * {@code transa} and {@code transb}. If no transposition of the source
     * matrices is required, {@code m} is the number of rows and {@code n} is
     * the number of columns in the source matrices {@code A} and {@code B}. In
     * this case, the output matrix {@code C} is {@code m x n}.
     * 
     * @param transa
     *            specifies the operation type for matrix {@code A}
     * @param transb
     *            specifies the operation type for matrix {@code B}
     * @param m
     *            the number of matrix rows in {@code op(A), op(B)}, and
     *            {@code C}
     * @param n
     *            the number of matrix columns in {@code op(A), op(B)}, and
     *            {@code C}
     * @param alpha
     *            this parameter scales the input matrix {@code A} by
     *            {@code alpha}
     * @param A
     *            the input matrix {@code A}, unmodified
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix {@code A}; measured in the number of elements.
     *            If {@code transa} is {@link TTrans#NO_TRANS} or
     *            {@link TTrans#CONJ}, this parameter must be at least
     *            {@code max(1, m)}; otherwise it must be {@code max(1, n)}
     * @param beta
     *            this parameter scales the input matrix {@code B} by
     *            {@code beta}
     * @param B
     *            the input matrix {@code B}, unmodified
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            source matrix {@code B}; measured in the number of elements.
     *            If {@code transb} is {@link TTrans#NO_TRANS} or
     *            {@link TTrans#CONJ}, this parameter must be at least
     *            {@code max(1, m)}; otherwise it must be {@code max(1, n)}
     * @param C
     *            the output matrix, modified
     * @param ldc
     *            distance between the first elements in adjacent columns in the
     *            destination matrix {@code C}; measured in the number of
     *            elements (it must be at least {@code max(1, m)})
     */
    public final void somatadd(TTrans transa, TTrans transb, int m, int n, float alpha, float[] A, int lda, float beta,
            float[] B, int ldb, float[] C, int ldc) {
        Objects.requireNonNull(transa, "transa");
        Objects.requireNonNull(transb, "transb");
        Objects.requireNonNull(A, "A");
        Objects.requireNonNull(B, "B");
        Objects.requireNonNull(C, "C");
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        somatadd_n(ORDERING, trans(transa), trans(transb), m, n, alpha, A, lda, beta, B, ldb, C, ldc, USE_CRITICAL);
    }

    /**
     * Scales and adds two doubles matrices, as well as performing out-of-place
     * transposition operations {@code C := alpha *op(A) + beta *op(B)} where
     * the {@code op()} operations are transpose, conjugate-transpose, conjugate
     * (no transpose), or no transpose, depending on the values of
     * {@code transa} and {@code transb}. If no transposition of the source
     * matrices is required, {@code m} is the number of rows and {@code n} is
     * the number of columns in the source matrices {@code A} and {@code B}. In
     * this case, the output matrix {@code C} is {@code m x n}.
     * 
     * @param transa
     *            specifies the operation type for matrix {@code A}
     * @param transb
     *            specifies the operation type for matrix {@code B}
     * @param m
     *            the number of matrix rows in {@code op(A), op(B)}, and
     *            {@code C}
     * @param n
     *            the number of matrix columns in {@code op(A), op(B)}, and
     *            {@code C}
     * @param alpha
     *            this parameter scales the input matrix {@code A} by
     *            {@code alpha}
     * @param A
     *            the input matrix {@code A}, unmodified
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix {@code A}; measured in the number of elements.
     *            If {@code transa} is {@link TTrans#NO_TRANS} or
     *            {@link TTrans#CONJ}, this parameter must be at least
     *            {@code max(1, m)}; otherwise it must be {@code max(1, n)}
     * @param beta
     *            this parameter scales the input matrix {@code B} by
     *            {@code beta}
     * @param B
     *            the input matrix {@code B}, unmodified
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            source matrix {@code B}; measured in the number of elements.
     *            If {@code transb} is {@link TTrans#NO_TRANS} or
     *            {@link TTrans#CONJ}, this parameter must be at least
     *            {@code max(1, m)}; otherwise it must be {@code max(1, n)}
     * @param C
     *            the output matrix, modified
     * @param ldc
     *            distance between the first elements in adjacent columns in the
     *            destination matrix {@code C}; measured in the number of
     *            elements (it must be at least {@code max(1, m)})
     */
    public final void domatadd(TTrans transa, TTrans transb, int m, int n, double alpha, double[] A, int lda,
            double beta, double[] B, int ldb, double[] C, int ldc) {
        Objects.requireNonNull(transa, "transa");
        Objects.requireNonNull(transb, "transb");
        Objects.requireNonNull(A, "A");
        Objects.requireNonNull(B, "B");
        Objects.requireNonNull(C, "C");
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        domatadd_n(ORDERING, trans(transa), trans(transb), m, n, alpha, A, lda, beta, B, ldb, C, ldc, USE_CRITICAL);
    }

    /**
     * Scales and adds two single precision complex matrices, as well as
     * performing out-of-place transposition operations
     * {@code C := alpha *op(A) + beta *op(B)} where the {@code op()} operations
     * are transpose, conjugate-transpose, conjugate (no transpose), or no
     * transpose, depending on the values of {@code transa} and {@code transb}.
     * If no transposition of the source matrices is required, {@code m} is the
     * number of rows and {@code n} is the number of columns in the source
     * matrices {@code A} and {@code B}. In this case, the output matrix
     * {@code C} is {@code m x n}.
     * 
     * @param transa
     *            specifies the operation type for matrix {@code A}
     * @param transb
     *            specifies the operation type for matrix {@code B}
     * @param m
     *            the number of matrix rows in {@code op(A), op(B)}, and
     *            {@code C}
     * @param n
     *            the number of matrix columns in {@code op(A), op(B)}, and
     *            {@code C}
     * @param alphar
     *            the real part of the scale factor for the input matrix
     *            {@code A}
     * @param alphai
     *            the imaginary part of the scale factor for the input matrix
     *            {@code A}
     * @param A
     *            the input matrix {@code A}, unmodified
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix {@code A}; measured in the number of elements.
     *            If {@code transa} is {@link TTrans#NO_TRANS} or
     *            {@link TTrans#CONJ}, this parameter must be at least
     *            {@code max(1, m)}; otherwise it must be {@code max(1, n)}
     * @param betar
     *            the real part of the scale factor for the input matrix
     *            {@code B}
     * @param betai
     *            the imaginary part of the scale factor for the input matrix
     *            {@code B}
     * @param B
     *            the input matrix {@code B}, unmodified
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            source matrix {@code B}; measured in the number of elements.
     *            If {@code transb} is {@link TTrans#NO_TRANS} or
     *            {@link TTrans#CONJ}, this parameter must be at least
     *            {@code max(1, m)}; otherwise it must be {@code max(1, n)}
     * @param C
     *            the output matrix, modified
     * @param ldc
     *            distance between the first elements in adjacent columns in the
     *            destination matrix {@code C}; measured in the number of
     *            elements (it must be at least {@code max(1, m)})
     */
    public final void comatadd(TTrans transa, TTrans transb, int m, int n, float alphar, float alphai, float[] A,
            int lda, float betar, float betai, float[] B, int ldb, float[] C, int ldc) {
        Objects.requireNonNull(transa, "transa");
        Objects.requireNonNull(transb, "transb");
        Objects.requireNonNull(A, "A");
        Objects.requireNonNull(B, "B");
        Objects.requireNonNull(C, "C");
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        comatadd_n(ORDERING, trans(transa), trans(transb), m, n, alphar, alphai, A, lda, betar, betai, B, ldb, C, ldc,
                USE_CRITICAL);
    }

    /**
     * Scales and adds two double precision complex matrices, as well as
     * performing out-of-place transposition operations
     * {@code C := alpha *op(A) + beta *op(B)} where the {@code op()} operations
     * are transpose, conjugate-transpose, conjugate (no transpose), or no
     * transpose, depending on the values of {@code transa} and {@code transb}.
     * If no transposition of the source matrices is required, {@code m} is the
     * number of rows and {@code n} is the number of columns in the source
     * matrices {@code A} and {@code B}. In this case, the output matrix
     * {@code C} is {@code m x n}.
     * 
     * @param transa
     *            specifies the operation type for matrix {@code A}
     * @param transb
     *            specifies the operation type for matrix {@code B}
     * @param m
     *            the number of matrix rows in {@code op(A), op(B)}, and
     *            {@code C}
     * @param n
     *            the number of matrix columns in {@code op(A), op(B)}, and
     *            {@code C}
     * @param alphar
     *            the real part of the scale factor for the input matrix
     *            {@code A}
     * @param alphai
     *            the imaginary part of the scale factor for the input matrix
     *            {@code A}
     * @param A
     *            the input matrix {@code A}, unmodified
     * @param lda
     *            distance between the first elements in adjacent columns in the
     *            source matrix {@code A}; measured in the number of elements.
     *            If {@code transa} is {@link TTrans#NO_TRANS} or
     *            {@link TTrans#CONJ}, this parameter must be at least
     *            {@code max(1, m)}; otherwise it must be {@code max(1, n)}
     * @param betar
     *            the real part of the scale factor for the input matrix
     *            {@code B}
     * @param betai
     *            the imaginary part of the scale factor for the input matrix
     *            {@code B}
     * @param B
     *            the input matrix {@code B}, unmodified
     * @param ldb
     *            distance between the first elements in adjacent columns in the
     *            source matrix {@code B}; measured in the number of elements.
     *            If {@code transb} is {@link TTrans#NO_TRANS} or
     *            {@link TTrans#CONJ}, this parameter must be at least
     *            {@code max(1, m)}; otherwise it must be {@code max(1, n)}
     * @param C
     *            the output matrix, modified
     * @param ldc
     *            distance between the first elements in adjacent columns in the
     *            destination matrix {@code C}; measured in the number of
     *            elements (it must be at least {@code max(1, m)})
     */
    public final void zomatadd(TTrans transa, TTrans transb, int m, int n, double alphar, double alphai, double[] A,
            int lda, double betar, double betai, double[] B, int ldb, double[] C, int ldc) {
        Objects.requireNonNull(transa, "transa");
        Objects.requireNonNull(transb, "transb");
        Objects.requireNonNull(A, "A");
        Objects.requireNonNull(B, "B");
        Objects.requireNonNull(C, "C");
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        zomatadd_n(ORDERING, trans(transa), trans(transb), m, n, alphar, alphai, A, lda, betar, betai, B, ldb, C, ldc,
                USE_CRITICAL);
    }

    private static native void cgemm3m_n(int order, int transa, int transb, int m, int n, int k, float alphar,
            float alphai, float[] a, int lda, float[] b, int ldb, float betar, float betai, float[] c, int ldc,
            boolean useCriticalRegion);

    private static native void zgemm3m_n(int order, int transa, int transb, int m, int n, int k, double alphar,
            double alphai, double[] a, int lda, double[] b, int ldb, double betar, double betai, double[] c, int ldc,
            boolean useCriticalRegion);

    private static native void simatcopy_n(byte ordering, byte trans, int rows, int cols, float alpha, float[] AB,
            int lda, int ldb, boolean useCriticalRegion);

    private static native void dimatcopy_n(byte ordering, byte trans, int rows, int cols, double alpha, double[] AB,
            int lda, int ldb, boolean useCriticalRegion);

    private static native void cimatcopy_n(byte ordering, byte trans, int rows, int cols, float alphar, float alphai,
            float[] AB, int lda, int ldb, boolean useCriticalRegion);

    private static native void zimatcopy_n(byte ordering, byte trans, int rows, int cols, double alphar, double alphai,
            double[] AB, int lda, int ldb, boolean useCriticalRegion);

    private static native void somatcopy_n(byte ordering, byte trans, int rows, int cols, float alpha, float[] A,
            int lda, float[] B, int ldb, boolean useCriticalRegion);

    private static native void domatcopy_n(byte ordering, byte trans, int rows, int cols, double alpha, double[] A,
            int lda, double[] B, int ldb, boolean useCriticalRegion);

    private static native void comatcopy_n(byte ordering, byte trans, int rows, int cols, float alphar, float alphai,
            float[] A, int lda, float[] B, int ldb, boolean useCriticalRegion);

    private static native void zomatcopy_n(byte ordering, byte trans, int rows, int cols, double alphar, double alphai,
            double[] A, int lda, double[] B, int ldb, boolean useCriticalRegion);

    private static native void somatadd_n(byte ordering, byte transa, byte transb, int m, int n, float alpha, float[] A,
            int lda, float beta, float[] B, int ldb, float[] C, int ldc, boolean useCriticalRegion);

    private static native void domatadd_n(byte ordering, byte transa, byte transb, int m, int n, double alpha,
            double[] A, int lda, double beta, double[] B, int ldb, double[] C, int ldc, boolean useCriticalRegion);

    private static native void comatadd_n(byte ordering, byte transa, byte transb, int m, int n, float alphar,
            float alphai, float[] A, int lda, float betar, float betai, float[] B, int ldb, float[] C, int ldc,
            boolean useCriticalRegion);

    private static native void zomatadd_n(byte ordering, byte transa, byte transb, int m, int n, double alphar,
            double alphai, double[] A, int lda, double betar, double betai, double[] B, int ldb, double[] C, int ldc,
            boolean useCriticalRegion);

    private static byte trans(TTrans trans) {
        return (byte) Character.toUpperCase(trans.val().charAt(0));
    }

    private static void checkStrictlyPositive(int value, String name) {
        if (value <= 0) {
            throw new IllegalArgumentException(
                    "Parameter " + name + " must be strictly positive (value = " + value + ")");
        }
    }

    protected BlasExt() {
    }
}
