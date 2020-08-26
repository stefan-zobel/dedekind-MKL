/*
 * Copyright 2019, 2020 Stefan Zobel
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

/**
 * Matrix storage layout must be column-major as in Fortran. All operations
 * throw a {@code NullPointerException} if any of the reference method arguments
 * is {@code null}.
 */
public abstract class Blas {

    private static final Blas mkl;
    private static final Blas netlib;

    static {
        netlib = new BlasJ();
        mkl = loadNative();
    }

    private static Blas loadNative() {
        try {
            System.loadLibrary("dedekind-mkl");
            BlasN.redirect_xerbla_n();
            return new BlasN();
        } catch (Throwable t) {
            t.printStackTrace();
        }
        return null;
    }

    /**
     * Returns an instance of a {@code Blas} implementation. Prefers the
     * {@code Intel MKL} implementation if that is available, otherwise the
     * {@code Netlib} implementation is used.
     * 
     * @return an instance of a {@code Blas} implementation
     */
    public static Blas getInstance() {
        Blas impl = mkl;
        if (impl != null) {
            return impl;
        }
        return netlib;
    }

    /**
     * Returns an instance of a specific {@code Blas} implementation depending
     * on the the value of the {@code useMKL} parameter. If {@code useMKL} is
     * {@code false} then the {@code Netlib} implementation (which is always
     * available) is used. If {@code useMKL} is {@code true} and the
     * {@code Intel MKL} implementation is available then that one will be used,
     * otherwise a RuntimeException is thrown.
     * 
     * @param useMKL
     *            whether to use the {@code Intel MKL} implementation (if
     *            {@code true}), or not
     * @return the requested {@code Blas} implementation if possible
     * @throws RuntimeException
     *             when {@code useMKL} is {@code true} and the {@code Intel MKL}
     *             implementation wasn't found
     */
    public static Blas getInstance(boolean useMKL) {
        if (!useMKL) {
            return netlib;
        }
        Blas impl = mkl;
        if (impl != null) {
            return impl;
        }
        throw new RuntimeException("MKL not loaded");
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGBMV  performs one of the matrix-vector operations
     *
     *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
     *
     *  where alpha and beta are scalars, x and y are vectors and A is an
     *  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
     *
     *  Arguments
     *  ==========
     *
     *  TRANS  - CHARACTER*1.
     *           On entry, TRANS specifies the operation to be performed as
     *           follows:
     *
     *              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
     *
     *              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
     *
     *              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
     *
     *           Unchanged on exit.
     *
     *  M      - INTEGER.
     *           On entry, M specifies the number of rows of the matrix A.
     *           M must be at least zero.
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the number of columns of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  KL     - INTEGER.
     *           On entry, KL specifies the number of sub-diagonals of the
     *           matrix A. KL must satisfy  0 .le. KL.
     *           Unchanged on exit.
     *
     *  KU     - INTEGER.
     *           On entry, KU specifies the number of super-diagonals of the
     *           matrix A. KU must satisfy  0 .le. KU.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
     *           Before entry, the leading ( kl + ku + 1 ) by n part of the
     *           array A must contain the matrix of coefficients, supplied
     *           column by column, with the leading diagonal of the matrix in
     *           row ( ku + 1 ) of the array, the first super-diagonal
     *           starting at position 2 in row ku, the first sub-diagonal
     *           starting at position 1 in row ( ku + 2 ), and so on.
     *           Elements in the array A that do not correspond to elements
     *           in the band matrix (such as the top left ku by ku triangle)
     *           are not referenced.
     *           The following program segment will transfer a band matrix
     *           from conventional full matrix storage to band storage:
     *
     *                 DO 20, J = 1, N
     *                    K = KU + 1 - J
     *                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
     *                       A( K + I, J ) = matrix( I, J )
     *              10    CONTINUE
     *              20 CONTINUE
     *
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. LDA must be at least
     *           ( kl + ku + 1 ).
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
     *           Before entry, the incremented array X must contain the
     *           vector x.
     *           Unchanged on exit.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *  BETA   - DOUBLE PRECISION.
     *           On entry, BETA specifies the scalar beta. When BETA is
     *           supplied as zero then Y need not be set on input.
     *           Unchanged on exit.
     *
     *  Y      - DOUBLE PRECISION array of DIMENSION at least
     *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
     *           Before entry, the incremented array Y must contain the
     *           vector y. On exit, Y is overwritten by the updated vector y.
     *
     *  INCY   - INTEGER.
     *           On entry, INCY specifies the increment for the elements of
     *           Y. INCY must not be zero.
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param trans
     * @param m
     * @param n
     * @param kl
     * @param ku
     * @param alpha
     * @param a
     * @param lda
     * @param x
     * @param incx
     * @param beta
     * @param y
     * @param incy
     */
    public final void dgbmv(String trans, int m, int n, int kl, int ku, double alpha, double[] a, int lda, double[] x,
            int incx, double beta, double[] y, int incy) {
        dgbmv(trans, m, n, kl, ku, alpha, a, 0, lda, x, 0, incx, beta, y, 0, incy);
    }

    public abstract void dgbmv(String trans, int m, int n, int kl, int ku, double alpha, double[] a, int aOffset,
            int lda, double[] x, int xOffset, int incx, double beta, double[] y, int yOffset, int incy);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGEMM  performs one of the matrix-matrix operations
     *
     *     C := alpha*op( A )*op( B ) + beta*C,
     *
     *  where  op( X ) is one of
     *
     *     op( X ) = X   or   op( X ) = X',
     *
     *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
     *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
     *
     *  Arguments
     *  ==========
     *
     *  TRANSA - CHARACTER*1.
     *           On entry, TRANSA specifies the form of op( A ) to be used in
     *           the matrix multiplication as follows:
     *
     *              TRANSA = 'N' or 'n',  op( A ) = A.
     *
     *              TRANSA = 'T' or 't',  op( A ) = A'.
     *
     *              TRANSA = 'C' or 'c',  op( A ) = A'.
     *
     *           Unchanged on exit.
     *
     *  TRANSB - CHARACTER*1.
     *           On entry, TRANSB specifies the form of op( B ) to be used in
     *           the matrix multiplication as follows:
     *
     *              TRANSB = 'N' or 'n',  op( B ) = B.
     *
     *              TRANSB = 'T' or 't',  op( B ) = B'.
     *
     *              TRANSB = 'C' or 'c',  op( B ) = B'.
     *
     *           Unchanged on exit.
     *
     *  M      - INTEGER.
     *           On entry,  M  specifies  the number  of rows  of the  matrix
     *           op( A )  and of the  matrix  C.  M  must  be at least  zero.
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry,  N  specifies the number  of columns of the matrix
     *           op( B ) and the number of columns of the matrix C. N must be
     *           at least zero.
     *           Unchanged on exit.
     *
     *  K      - INTEGER.
     *           On entry,  K  specifies  the number of columns of the matrix
     *           op( A ) and the number of rows of the matrix op( B ). K must
     *           be at least  zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
     *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by m  part of the array  A  must contain  the
     *           matrix A.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
     *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
     *           least  max( 1, k ).
     *           Unchanged on exit.
     *
     *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
     *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
     *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  n by k  part of the array  B  must contain  the
     *           matrix B.
     *           Unchanged on exit.
     *
     *  LDB    - INTEGER.
     *           On entry, LDB specifies the first dimension of B as declared
     *           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
     *           LDB must be at least  max( 1, k ), otherwise  LDB must be at
     *           least  max( 1, n ).
     *           Unchanged on exit.
     *
     *  BETA   - DOUBLE PRECISION.
     *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
     *           supplied as zero then C need not be set on input.
     *           Unchanged on exit.
     *
     *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     *           On exit, the array  C  is overwritten by the  m by n  matrix
     *           ( alpha*op( A )*op( B ) + beta*C ).
     *
     *  LDC    - INTEGER.
     *           On entry, LDC specifies the first dimension of C as declared
     *           in  the  calling  (sub)  program.   LDC  must  be  at  least
     *           max( 1, m ).
     *           Unchanged on exit.
     *
     *
     *  Level 3 Blas routine.
     *
     *  -- Written on 8-February-1989.
     *     Jack Dongarra, Argonne National Laboratory.
     *     Iain Duff, AERE Harwell.
     *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
     *     Sven Hammarling, Numerical Algorithms Group Ltd.
     *
     * </code>
     * </pre>
     *
     * @param transa
     * @param transb
     * @param m
     * @param n
     * @param k
     * @param alpha
     * @param a
     * @param lda
     * @param b
     * @param ldb
     * @param beta
     * @param c
     * @param ldc
     */
    public final void dgemm(String transa, String transb, int m, int n, int k, double alpha, double[] a, int lda,
            double[] b, int ldb, double beta, double[] c, int ldc) {
        dgemm(transa, transb, m, n, k, alpha, a, 0, lda, b, 0, ldb, beta, c, 0, ldc);
    }

    public abstract void dgemm(String transa, String transb, int m, int n, int k, double alpha, double[] a, int aOffset,
            int lda, double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGEMV  performs one of the matrix-vector operations
     *
     *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
     *
     *  where alpha and beta are scalars, x and y are vectors and A is an
     *  m by n matrix.
     *
     *  Arguments
     *  ==========
     *
     *  TRANS  - CHARACTER*1.
     *           On entry, TRANS specifies the operation to be performed as
     *           follows:
     *
     *              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
     *
     *              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
     *
     *              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
     *
     *           Unchanged on exit.
     *
     *  M      - INTEGER.
     *           On entry, M specifies the number of rows of the matrix A.
     *           M must be at least zero.
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the number of columns of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
     *           Before entry, the leading m by n part of the array A must
     *           contain the matrix of coefficients.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. LDA must be at least
     *           max( 1, m ).
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
     *           Before entry, the incremented array X must contain the
     *           vector x.
     *           Unchanged on exit.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *  BETA   - DOUBLE PRECISION.
     *           On entry, BETA specifies the scalar beta. When BETA is
     *           supplied as zero then Y need not be set on input.
     *           Unchanged on exit.
     *
     *  Y      - DOUBLE PRECISION array of DIMENSION at least
     *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
     *           Before entry with BETA non-zero, the incremented array Y
     *           must contain the vector y. On exit, Y is overwritten by the
     *           updated vector y.
     *
     *  INCY   - INTEGER.
     *           On entry, INCY specifies the increment for the elements of
     *           Y. INCY must not be zero.
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param trans
     * @param m
     * @param n
     * @param alpha
     * @param a
     * @param lda
     * @param x
     * @param incx
     * @param beta
     * @param y
     * @param incy
     */
    public final void dgemv(String trans, int m, int n, double alpha, double[] a, int lda, double[] x, int incx,
            double beta, double[] y, int incy) {
        dgemv(trans, m, n, alpha, a, 0, lda, x, 0, incx, beta, y, 0, incy);
    }

    public abstract void dgemv(String trans, int m, int n, double alpha, double[] a, int aOffset, int lda, double[] x,
            int xOffset, int incx, double beta, double[] y, int yOffset, int incy);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGER   performs the rank 1 operation
     *
     *     A := alpha*x*y' + A,
     *
     *  where alpha is a scalar, x is an m element vector, y is an n element
     *  vector and A is an m by n matrix.
     *
     *  Arguments
     *  ==========
     *
     *  M      - INTEGER.
     *           On entry, M specifies the number of rows of the matrix A.
     *           M must be at least zero.
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the number of columns of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( m - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the m
     *           element vector x.
     *           Unchanged on exit.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *  Y      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *           Before entry, the incremented array Y must contain the n
     *           element vector y.
     *           Unchanged on exit.
     *
     *  INCY   - INTEGER.
     *           On entry, INCY specifies the increment for the elements of
     *           Y. INCY must not be zero.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
     *           Before entry, the leading m by n part of the array A must
     *           contain the matrix of coefficients. On exit, A is
     *           overwritten by the updated matrix.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. LDA must be at least
     *           max( 1, m ).
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param alpha
     * @param x
     * @param incx
     * @param y
     * @param incy
     * @param a
     * @param lda
     */
    public final void dger(int m, int n, double alpha, double[] x, int incx, double[] y, int incy, double[] a,
            int lda) {
        dger(m, n, alpha, x, 0, incx, y, 0, incy, a, 0, lda);
    }

    public abstract void dger(int m, int n, double alpha, double[] x, int xOffset, int incx, double[] y, int yOffset,
            int incy, double[] a, int aOffset, int lda);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSBMV  performs the matrix-vector  operation
     *
     *     y := alpha*A*x + beta*y,
     *
     *  where alpha and beta are scalars, x and y are n element vectors and
     *  A is an n by n symmetric band matrix, with k super-diagonals.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the upper or lower
     *           triangular part of the band matrix A is being supplied as
     *           follows:
     *
     *              UPLO = 'U' or 'u'   The upper triangular part of A is
     *                                  being supplied.
     *
     *              UPLO = 'L' or 'l'   The lower triangular part of A is
     *                                  being supplied.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  K      - INTEGER.
     *           On entry, K specifies the number of super-diagonals of the
     *           matrix A. K must satisfy  0 .le. K.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
     *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
     *           by n part of the array A must contain the upper triangular
     *           band part of the symmetric matrix, supplied column by
     *           column, with the leading diagonal of the matrix in row
     *           ( k + 1 ) of the array, the first super-diagonal starting at
     *           position 2 in row k, and so on. The top left k by k triangle
     *           of the array A is not referenced.
     *           The following program segment will transfer the upper
     *           triangular part of a symmetric band matrix from conventional
     *           full matrix storage to band storage:
     *
     *                 DO 20, J = 1, N
     *                    M = K + 1 - J
     *                    DO 10, I = MAX( 1, J - K ), J
     *                       A( M + I, J ) = matrix( I, J )
     *              10    CONTINUE
     *              20 CONTINUE
     *
     *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
     *           by n part of the array A must contain the lower triangular
     *           band part of the symmetric matrix, supplied column by
     *           column, with the leading diagonal of the matrix in row 1 of
     *           the array, the first sub-diagonal starting at position 1 in
     *           row 2, and so on. The bottom right k by k triangle of the
     *           array A is not referenced.
     *           The following program segment will transfer the lower
     *           triangular part of a symmetric band matrix from conventional
     *           full matrix storage to band storage:
     *
     *                 DO 20, J = 1, N
     *                    M = 1 - J
     *                    DO 10, I = J, MIN( N, J + K )
     *                       A( M + I, J ) = matrix( I, J )
     *              10    CONTINUE
     *              20 CONTINUE
     *
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. LDA must be at least
     *           ( k + 1 ).
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the
     *           vector x.
     *           Unchanged on exit.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *  BETA   - DOUBLE PRECISION.
     *           On entry, BETA specifies the scalar beta.
     *           Unchanged on exit.
     *
     *  Y      - DOUBLE PRECISION array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *           Before entry, the incremented array Y must contain the
     *           vector y. On exit, Y is overwritten by the updated vector y.
     *
     *  INCY   - INTEGER.
     *           On entry, INCY specifies the increment for the elements of
     *           Y. INCY must not be zero.
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param k
     * @param alpha
     * @param a
     * @param lda
     * @param x
     * @param incx
     * @param beta
     * @param y
     * @param incy
     */
    public final void dsbmv(String uplo, int n, int k, double alpha, double[] a, int lda, double[] x, int incx,
            double beta, double[] y, int incy) {
        dsbmv(uplo, n, k, alpha, a, 0, lda, x, 0, incx, beta, y, 0, incy);
    }

    public abstract void dsbmv(String uplo, int n, int k, double alpha, double[] a, int aOffset, int lda, double[] x,
            int xOffset, int incx, double beta, double[] y, int yOffset, int incy);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSPMV  performs the matrix-vector operation
     *
     *     y := alpha*A*x + beta*y,
     *
     *  where alpha and beta are scalars, x and y are n element vectors and
     *  A is an n by n symmetric matrix, supplied in packed form.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the upper or lower
     *           triangular part of the matrix A is supplied in the packed
     *           array AP as follows:
     *
     *              UPLO = 'U' or 'u'   The upper triangular part of A is
     *                                  supplied in AP.
     *
     *              UPLO = 'L' or 'l'   The lower triangular part of A is
     *                                  supplied in AP.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  AP     - DOUBLE PRECISION array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     *           Before entry with UPLO = 'U' or 'u', the array AP must
     *           contain the upper triangular part of the symmetric matrix
     *           packed sequentially, column by column, so that AP( 1 )
     *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
     *           and a( 2, 2 ) respectively, and so on.
     *           Before entry with UPLO = 'L' or 'l', the array AP must
     *           contain the lower triangular part of the symmetric matrix
     *           packed sequentially, column by column, so that AP( 1 )
     *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
     *           and a( 3, 1 ) respectively, and so on.
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element vector x.
     *           Unchanged on exit.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *  BETA   - DOUBLE PRECISION.
     *           On entry, BETA specifies the scalar beta. When BETA is
     *           supplied as zero then Y need not be set on input.
     *           Unchanged on exit.
     *
     *  Y      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *           Before entry, the incremented array Y must contain the n
     *           element vector y. On exit, Y is overwritten by the updated
     *           vector y.
     *
     *  INCY   - INTEGER.
     *           On entry, INCY specifies the increment for the elements of
     *           Y. INCY must not be zero.
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param alpha
     * @param ap
     * @param x
     * @param incx
     * @param beta
     * @param y
     * @param incy
     */
    public final void dspmv(String uplo, int n, double alpha, double[] ap, double[] x, int incx, double beta,
            double[] y, int incy) {
        dspmv(uplo, n, alpha, ap, 0, x, 0, incx, beta, y, 0, incy);
    }

    public abstract void dspmv(String uplo, int n, double alpha, double[] ap, int apOffset, double[] x, int xOffset,
            int incx, double beta, double[] y, int yOffset, int incy);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSPR    performs the symmetric rank 1 operation
     *
     *     A := alpha*x*x' + A,
     *
     *  where alpha is a real scalar, x is an n element vector and A is an
     *  n by n symmetric matrix, supplied in packed form.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the upper or lower
     *           triangular part of the matrix A is supplied in the packed
     *           array AP as follows:
     *
     *              UPLO = 'U' or 'u'   The upper triangular part of A is
     *                                  supplied in AP.
     *
     *              UPLO = 'L' or 'l'   The lower triangular part of A is
     *                                  supplied in AP.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element vector x.
     *           Unchanged on exit.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *  AP     - DOUBLE PRECISION array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     *           Before entry with  UPLO = 'U' or 'u', the array AP must
     *           contain the upper triangular part of the symmetric matrix
     *           packed sequentially, column by column, so that AP( 1 )
     *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
     *           and a( 2, 2 ) respectively, and so on. On exit, the array
     *           AP is overwritten by the upper triangular part of the
     *           updated matrix.
     *           Before entry with UPLO = 'L' or 'l', the array AP must
     *           contain the lower triangular part of the symmetric matrix
     *           packed sequentially, column by column, so that AP( 1 )
     *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
     *           and a( 3, 1 ) respectively, and so on. On exit, the array
     *           AP is overwritten by the lower triangular part of the
     *           updated matrix.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param alpha
     * @param x
     * @param incx
     * @param ap
     */
    public final void dspr(String uplo, int n, double alpha, double[] x, int incx, double[] ap) {
        dspr(uplo, n, alpha, x, 0, incx, ap, 0);
    }

    public abstract void dspr(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] ap,
            int apOffset);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSPR2  performs the symmetric rank 2 operation
     *
     *     A := alpha*x*y' + alpha*y*x' + A,
     *
     *  where alpha is a scalar, x and y are n element vectors and A is an
     *  n by n symmetric matrix, supplied in packed form.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the upper or lower
     *           triangular part of the matrix A is supplied in the packed
     *           array AP as follows:
     *
     *              UPLO = 'U' or 'u'   The upper triangular part of A is
     *                                  supplied in AP.
     *
     *              UPLO = 'L' or 'l'   The lower triangular part of A is
     *                                  supplied in AP.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element vector x.
     *           Unchanged on exit.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *  Y      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *           Before entry, the incremented array Y must contain the n
     *           element vector y.
     *           Unchanged on exit.
     *
     *  INCY   - INTEGER.
     *           On entry, INCY specifies the increment for the elements of
     *           Y. INCY must not be zero.
     *           Unchanged on exit.
     *
     *  AP     - DOUBLE PRECISION array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     *           Before entry with  UPLO = 'U' or 'u', the array AP must
     *           contain the upper triangular part of the symmetric matrix
     *           packed sequentially, column by column, so that AP( 1 )
     *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
     *           and a( 2, 2 ) respectively, and so on. On exit, the array
     *           AP is overwritten by the upper triangular part of the
     *           updated matrix.
     *           Before entry with UPLO = 'L' or 'l', the array AP must
     *           contain the lower triangular part of the symmetric matrix
     *           packed sequentially, column by column, so that AP( 1 )
     *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
     *           and a( 3, 1 ) respectively, and so on. On exit, the array
     *           AP is overwritten by the lower triangular part of the
     *           updated matrix.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param alpha
     * @param x
     * @param incx
     * @param y
     * @param incy
     * @param ap
     */
    public final void dspr2(String uplo, int n, double alpha, double[] x, int incx, double[] y, int incy, double[] ap) {
        dspr2(uplo, n, alpha, x, 0, incx, y, 0, incy, ap, 0);
    }

    public abstract void dspr2(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] y,
            int yOffset, int incy, double[] ap, int apOffset);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSYMM  performs one of the matrix-matrix operations
     *
     *     C := alpha*A*B + beta*C,
     *
     *  or
     *
     *     C := alpha*B*A + beta*C,
     *
     *  where alpha and beta are scalars,  A is a symmetric matrix and  B and
     *  C are  m by n matrices.
     *
     *  Arguments
     *  ==========
     *
     *  SIDE   - CHARACTER*1.
     *           On entry,  SIDE  specifies whether  the  symmetric matrix  A
     *           appears on the  left or right  in the  operation as follows:
     *
     *              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
     *
     *              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
     *
     *           Unchanged on exit.
     *
     *  UPLO   - CHARACTER*1.
     *           On  entry,   UPLO  specifies  whether  the  upper  or  lower
     *           triangular  part  of  the  symmetric  matrix   A  is  to  be
     *           referenced as follows:
     *
     *              UPLO = 'U' or 'u'   Only the upper triangular part of the
     *                                  symmetric matrix is to be referenced.
     *
     *              UPLO = 'L' or 'l'   Only the lower triangular part of the
     *                                  symmetric matrix is to be referenced.
     *
     *           Unchanged on exit.
     *
     *  M      - INTEGER.
     *           On entry,  M  specifies the number of rows of the matrix  C.
     *           M  must be at least zero.
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the number of columns of the matrix C.
     *           N  must be at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
     *           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
     *           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
     *           the array  A  must contain the  symmetric matrix,  such that
     *           when  UPLO = 'U' or 'u', the leading m by m upper triangular
     *           part of the array  A  must contain the upper triangular part
     *           of the  symmetric matrix and the  strictly  lower triangular
     *           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
     *           the leading  m by m  lower triangular part  of the  array  A
     *           must  contain  the  lower triangular part  of the  symmetric
     *           matrix and the  strictly upper triangular part of  A  is not
     *           referenced.
     *           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
     *           the array  A  must contain the  symmetric matrix,  such that
     *           when  UPLO = 'U' or 'u', the leading n by n upper triangular
     *           part of the array  A  must contain the upper triangular part
     *           of the  symmetric matrix and the  strictly  lower triangular
     *           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
     *           the leading  n by n  lower triangular part  of the  array  A
     *           must  contain  the  lower triangular part  of the  symmetric
     *           matrix and the  strictly upper triangular part of  A  is not
     *           referenced.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
     *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
     *           least  max( 1, n ).
     *           Unchanged on exit.
     *
     *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
     *           Before entry, the leading  m by n part of the array  B  must
     *           contain the matrix B.
     *           Unchanged on exit.
     *
     *  LDB    - INTEGER.
     *           On entry, LDB specifies the first dimension of B as declared
     *           in  the  calling  (sub)  program.   LDB  must  be  at  least
     *           max( 1, m ).
     *           Unchanged on exit.
     *
     *  BETA   - DOUBLE PRECISION.
     *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
     *           supplied as zero then C need not be set on input.
     *           Unchanged on exit.
     *
     *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     *           On exit, the array  C  is overwritten by the  m by n updated
     *           matrix.
     *
     *  LDC    - INTEGER.
     *           On entry, LDC specifies the first dimension of C as declared
     *           in  the  calling  (sub)  program.   LDC  must  be  at  least
     *           max( 1, m ).
     *           Unchanged on exit.
     *
     *
     *  Level 3 Blas routine.
     *
     *  -- Written on 8-February-1989.
     *     Jack Dongarra, Argonne National Laboratory.
     *     Iain Duff, AERE Harwell.
     *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
     *     Sven Hammarling, Numerical Algorithms Group Ltd.
     *
     * </code>
     * </pre>
     *
     * @param side
     * @param uplo
     * @param m
     * @param n
     * @param alpha
     * @param a
     * @param lda
     * @param b
     * @param ldb
     * @param beta
     * @param c
     * @param ldc
     */
    public final void dsymm(String side, String uplo, int m, int n, double alpha, double[] a, int lda, double[] b,
            int ldb, double beta, double[] c, int ldc) {
        dsymm(side, uplo, m, n, alpha, a, 0, lda, b, 0, ldb, beta, c, 0, ldc);
    }

    public abstract void dsymm(String side, String uplo, int m, int n, double alpha, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSYMV  performs the matrix-vector  operation
     *
     *     y := alpha*A*x + beta*y,
     *
     *  where alpha and beta are scalars, x and y are n element vectors and
     *  A is an n by n symmetric matrix.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the upper or lower
     *           triangular part of the array A is to be referenced as
     *           follows:
     *
     *              UPLO = 'U' or 'u'   Only the upper triangular part of A
     *                                  is to be referenced.
     *
     *              UPLO = 'L' or 'l'   Only the lower triangular part of A
     *                                  is to be referenced.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     *           triangular part of the symmetric matrix and the strictly
     *           lower triangular part of A is not referenced.
     *           Before entry with UPLO = 'L' or 'l', the leading n by n
     *           lower triangular part of the array A must contain the lower
     *           triangular part of the symmetric matrix and the strictly
     *           upper triangular part of A is not referenced.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. LDA must be at least
     *           max( 1, n ).
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element vector x.
     *           Unchanged on exit.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *  BETA   - DOUBLE PRECISION.
     *           On entry, BETA specifies the scalar beta. When BETA is
     *           supplied as zero then Y need not be set on input.
     *           Unchanged on exit.
     *
     *  Y      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *           Before entry, the incremented array Y must contain the n
     *           element vector y. On exit, Y is overwritten by the updated
     *           vector y.
     *
     *  INCY   - INTEGER.
     *           On entry, INCY specifies the increment for the elements of
     *           Y. INCY must not be zero.
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param alpha
     * @param a
     * @param lda
     * @param x
     * @param incx
     * @param beta
     * @param y
     * @param incy
     */
    public final void dsymv(String uplo, int n, double alpha, double[] a, int lda, double[] x, int incx, double beta,
            double[] y, int incy) {
        dsymv(uplo, n, alpha, a, 0, lda, x, 0, incx, beta, y, 0, incy);
    }

    public abstract void dsymv(String uplo, int n, double alpha, double[] a, int aOffset, int lda, double[] x,
            int xOffset, int incx, double beta, double[] y, int yOffset, int incy);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSYR   performs the symmetric rank 1 operation
     *
     *     A := alpha*x*x' + A,
     *
     *  where alpha is a real scalar, x is an n element vector and A is an
     *  n by n symmetric matrix.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the upper or lower
     *           triangular part of the array A is to be referenced as
     *           follows:
     *
     *              UPLO = 'U' or 'u'   Only the upper triangular part of A
     *                                  is to be referenced.
     *
     *              UPLO = 'L' or 'l'   Only the lower triangular part of A
     *                                  is to be referenced.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element vector x.
     *           Unchanged on exit.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     *           triangular part of the symmetric matrix and the strictly
     *           lower triangular part of A is not referenced. On exit, the
     *           upper triangular part of the array A is overwritten by the
     *           upper triangular part of the updated matrix.
     *           Before entry with UPLO = 'L' or 'l', the leading n by n
     *           lower triangular part of the array A must contain the lower
     *           triangular part of the symmetric matrix and the strictly
     *           upper triangular part of A is not referenced. On exit, the
     *           lower triangular part of the array A is overwritten by the
     *           lower triangular part of the updated matrix.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. LDA must be at least
     *           max( 1, n ).
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param alpha
     * @param x
     * @param incx
     * @param a
     * @param lda
     */
    public final void dsyr(String uplo, int n, double alpha, double[] x, int incx, double[] a, int lda) {
        dsyr(uplo, n, alpha, x, 0, incx, a, 0, lda);
    }

    public abstract void dsyr(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] a,
            int aOffset, int lda);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSYR2  performs the symmetric rank 2 operation
     *
     *     A := alpha*x*y' + alpha*y*x' + A,
     *
     *  where alpha is a scalar, x and y are n element vectors and A is an n
     *  by n symmetric matrix.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the upper or lower
     *           triangular part of the array A is to be referenced as
     *           follows:
     *
     *              UPLO = 'U' or 'u'   Only the upper triangular part of A
     *                                  is to be referenced.
     *
     *              UPLO = 'L' or 'l'   Only the lower triangular part of A
     *                                  is to be referenced.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element vector x.
     *           Unchanged on exit.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *  Y      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *           Before entry, the incremented array Y must contain the n
     *           element vector y.
     *           Unchanged on exit.
     *
     *  INCY   - INTEGER.
     *           On entry, INCY specifies the increment for the elements of
     *           Y. INCY must not be zero.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     *           triangular part of the symmetric matrix and the strictly
     *           lower triangular part of A is not referenced. On exit, the
     *           upper triangular part of the array A is overwritten by the
     *           upper triangular part of the updated matrix.
     *           Before entry with UPLO = 'L' or 'l', the leading n by n
     *           lower triangular part of the array A must contain the lower
     *           triangular part of the symmetric matrix and the strictly
     *           upper triangular part of A is not referenced. On exit, the
     *           lower triangular part of the array A is overwritten by the
     *           lower triangular part of the updated matrix.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. LDA must be at least
     *           max( 1, n ).
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param alpha
     * @param x
     * @param incx
     * @param y
     * @param incy
     * @param a
     * @param lda
     */
    public final void dsyr2(String uplo, int n, double alpha, double[] x, int incx, double[] y, int incy, double[] a,
            int lda) {
        dsyr2(uplo, n, alpha, x, 0, incx, y, 0, incy, a, 0, lda);
    }

    public abstract void dsyr2(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] y,
            int yOffset, int incy, double[] a, int aOffset, int lda);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSYR2K  performs one of the symmetric rank 2k operations
     *
     *     C := alpha*A*B' + alpha*B*A' + beta*C,
     *
     *  or
     *
     *     C := alpha*A'*B + alpha*B'*A + beta*C,
     *
     *  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     *  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     *  matrices in the second case.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On  entry,   UPLO  specifies  whether  the  upper  or  lower
     *           triangular  part  of the  array  C  is to be  referenced  as
     *           follows:
     *
     *              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
     *                                  is to be referenced.
     *
     *              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
     *                                  is to be referenced.
     *
     *           Unchanged on exit.
     *
     *  TRANS  - CHARACTER*1.
     *           On entry,  TRANS  specifies the operation to be performed as
     *           follows:
     *
     *              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
     *                                        beta*C.
     *
     *              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
     *                                        beta*C.
     *
     *              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
     *                                        beta*C.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry,  N specifies the order of the matrix C.  N must be
     *           at least zero.
     *           Unchanged on exit.
     *
     *  K      - INTEGER.
     *           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
     *           of  columns  of the  matrices  A and B,  and on  entry  with
     *           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
     *           of rows of the matrices  A and B.  K must be at least  zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
     *           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by n  part of the array  A  must contain  the
     *           matrix A.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
     *           then  LDA must be at least  max( 1, n ), otherwise  LDA must
     *           be at least  max( 1, k ).
     *           Unchanged on exit.
     *
     *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
     *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
     *           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  k by n  part of the array  B  must contain  the
     *           matrix B.
     *           Unchanged on exit.
     *
     *  LDB    - INTEGER.
     *           On entry, LDB specifies the first dimension of B as declared
     *           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
     *           then  LDB must be at least  max( 1, n ), otherwise  LDB must
     *           be at least  max( 1, k ).
     *           Unchanged on exit.
     *
     *  BETA   - DOUBLE PRECISION.
     *           On entry, BETA specifies the scalar beta.
     *           Unchanged on exit.
     *
     *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
     *           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
     *           upper triangular part of the array C must contain the upper
     *           triangular part  of the  symmetric matrix  and the strictly
     *           lower triangular part of C is not referenced.  On exit, the
     *           upper triangular part of the array  C is overwritten by the
     *           upper triangular part of the updated matrix.
     *           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
     *           lower triangular part of the array C must contain the lower
     *           triangular part  of the  symmetric matrix  and the strictly
     *           upper triangular part of C is not referenced.  On exit, the
     *           lower triangular part of the array  C is overwritten by the
     *           lower triangular part of the updated matrix.
     *
     *  LDC    - INTEGER.
     *           On entry, LDC specifies the first dimension of C as declared
     *           in  the  calling  (sub)  program.   LDC  must  be  at  least
     *           max( 1, n ).
     *           Unchanged on exit.
     *
     *
     *  Level 3 Blas routine.
     *
     *
     *  -- Written on 8-February-1989.
     *     Jack Dongarra, Argonne National Laboratory.
     *     Iain Duff, AERE Harwell.
     *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
     *     Sven Hammarling, Numerical Algorithms Group Ltd.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param trans
     * @param n
     * @param k
     * @param alpha
     * @param a
     * @param lda
     * @param b
     * @param ldb
     * @param beta
     * @param c
     * @param ldc
     */
    public final void dsyr2k(String uplo, String trans, int n, int k, double alpha, double[] a, int lda, double[] b,
            int ldb, double beta, double[] c, int ldc) {
        dsyr2k(uplo, trans, n, k, alpha, a, 0, lda, b, 0, ldb, beta, c, 0, ldc);
    }

    public abstract void dsyr2k(String uplo, String trans, int n, int k, double alpha, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSYRK  performs one of the symmetric rank k operations
     *
     *     C := alpha*A*A' + beta*C,
     *
     *  or
     *
     *     C := alpha*A'*A + beta*C,
     *
     *  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     *  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     *  in the second case.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On  entry,   UPLO  specifies  whether  the  upper  or  lower
     *           triangular  part  of the  array  C  is to be  referenced  as
     *           follows:
     *
     *              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
     *                                  is to be referenced.
     *
     *              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
     *                                  is to be referenced.
     *
     *           Unchanged on exit.
     *
     *  TRANS  - CHARACTER*1.
     *           On entry,  TRANS  specifies the operation to be performed as
     *           follows:
     *
     *              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
     *
     *              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
     *
     *              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry,  N specifies the order of the matrix C.  N must be
     *           at least zero.
     *           Unchanged on exit.
     *
     *  K      - INTEGER.
     *           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
     *           of  columns   of  the   matrix   A,   and  on   entry   with
     *           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
     *           of rows of the matrix  A.  K must be at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
     *           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by n  part of the array  A  must contain  the
     *           matrix A.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
     *           then  LDA must be at least  max( 1, n ), otherwise  LDA must
     *           be at least  max( 1, k ).
     *           Unchanged on exit.
     *
     *  BETA   - DOUBLE PRECISION.
     *           On entry, BETA specifies the scalar beta.
     *           Unchanged on exit.
     *
     *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
     *           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
     *           upper triangular part of the array C must contain the upper
     *           triangular part  of the  symmetric matrix  and the strictly
     *           lower triangular part of C is not referenced.  On exit, the
     *           upper triangular part of the array  C is overwritten by the
     *           upper triangular part of the updated matrix.
     *           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
     *           lower triangular part of the array C must contain the lower
     *           triangular part  of the  symmetric matrix  and the strictly
     *           upper triangular part of C is not referenced.  On exit, the
     *           lower triangular part of the array  C is overwritten by the
     *           lower triangular part of the updated matrix.
     *
     *  LDC    - INTEGER.
     *           On entry, LDC specifies the first dimension of C as declared
     *           in  the  calling  (sub)  program.   LDC  must  be  at  least
     *           max( 1, n ).
     *           Unchanged on exit.
     *
     *
     *  Level 3 Blas routine.
     *
     *  -- Written on 8-February-1989.
     *     Jack Dongarra, Argonne National Laboratory.
     *     Iain Duff, AERE Harwell.
     *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
     *     Sven Hammarling, Numerical Algorithms Group Ltd.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param trans
     * @param n
     * @param k
     * @param alpha
     * @param a
     * @param lda
     * @param beta
     * @param c
     * @param ldc
     */
    public final void dsyrk(String uplo, String trans, int n, int k, double alpha, double[] a, int lda, double beta,
            double[] c, int ldc) {
        dsyrk(uplo, trans, n, k, alpha, a, 0, lda, beta, c, 0, ldc);
    }

    public abstract void dsyrk(String uplo, String trans, int n, int k, double alpha, double[] a, int aOffset, int lda,
            double beta, double[] c, int cOffset, int ldc);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DTBMV  performs one of the matrix-vector operations
     *
     *     x := A*x,   or   x := A'*x,
     *
     *  where x is an n element vector and  A is an n by n unit, or non-unit,
     *  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the matrix is an upper or
     *           lower triangular matrix as follows:
     *
     *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
     *
     *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
     *
     *           Unchanged on exit.
     *
     *  TRANS  - CHARACTER*1.
     *           On entry, TRANS specifies the operation to be performed as
     *           follows:
     *
     *              TRANS = 'N' or 'n'   x := A*x.
     *
     *              TRANS = 'T' or 't'   x := A'*x.
     *
     *              TRANS = 'C' or 'c'   x := A'*x.
     *
     *           Unchanged on exit.
     *
     *  DIAG   - CHARACTER*1.
     *           On entry, DIAG specifies whether or not A is unit
     *           triangular as follows:
     *
     *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
     *
     *              DIAG = 'N' or 'n'   A is not assumed to be unit
     *                                  triangular.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  K      - INTEGER.
     *           On entry with UPLO = 'U' or 'u', K specifies the number of
     *           super-diagonals of the matrix A.
     *           On entry with UPLO = 'L' or 'l', K specifies the number of
     *           sub-diagonals of the matrix A.
     *           K must satisfy  0 .le. K.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
     *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
     *           by n part of the array A must contain the upper triangular
     *           band part of the matrix of coefficients, supplied column by
     *           column, with the leading diagonal of the matrix in row
     *           ( k + 1 ) of the array, the first super-diagonal starting at
     *           position 2 in row k, and so on. The top left k by k triangle
     *           of the array A is not referenced.
     *           The following program segment will transfer an upper
     *           triangular band matrix from conventional full matrix storage
     *           to band storage:
     *
     *                 DO 20, J = 1, N
     *                    M = K + 1 - J
     *                    DO 10, I = MAX( 1, J - K ), J
     *                       A( M + I, J ) = matrix( I, J )
     *              10    CONTINUE
     *              20 CONTINUE
     *
     *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
     *           by n part of the array A must contain the lower triangular
     *           band part of the matrix of coefficients, supplied column by
     *           column, with the leading diagonal of the matrix in row 1 of
     *           the array, the first sub-diagonal starting at position 1 in
     *           row 2, and so on. The bottom right k by k triangle of the
     *           array A is not referenced.
     *           The following program segment will transfer a lower
     *           triangular band matrix from conventional full matrix storage
     *           to band storage:
     *
     *                 DO 20, J = 1, N
     *                    M = 1 - J
     *                    DO 10, I = J, MIN( N, J + K )
     *                       A( M + I, J ) = matrix( I, J )
     *              10    CONTINUE
     *              20 CONTINUE
     *
     *           Note that when DIAG = 'U' or 'u' the elements of the array A
     *           corresponding to the diagonal elements of the matrix are not
     *           referenced, but are assumed to be unity.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. LDA must be at least
     *           ( k + 1 ).
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element vector x. On exit, X is overwritten with the
     *           tranformed vector x.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param trans
     * @param diag
     * @param n
     * @param k
     * @param a
     * @param lda
     * @param x
     * @param incx
     */
    public final void dtbmv(String uplo, String trans, String diag, int n, int k, double[] a, int lda, double[] x,
            int incx) {
        dtbmv(uplo, trans, diag, n, k, a, 0, lda, x, 0, incx);
    }

    public abstract void dtbmv(String uplo, String trans, String diag, int n, int k, double[] a, int aOffset, int lda,
            double[] x, int xOffset, int incx);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DTPMV  performs one of the matrix-vector operations
     *
     *     x := A*x,   or   x := A'*x,
     *
     *  where x is an n element vector and  A is an n by n unit, or non-unit,
     *  upper or lower triangular matrix, supplied in packed form.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the matrix is an upper or
     *           lower triangular matrix as follows:
     *
     *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
     *
     *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
     *
     *           Unchanged on exit.
     *
     *  TRANS  - CHARACTER*1.
     *           On entry, TRANS specifies the operation to be performed as
     *           follows:
     *
     *              TRANS = 'N' or 'n'   x := A*x.
     *
     *              TRANS = 'T' or 't'   x := A'*x.
     *
     *              TRANS = 'C' or 'c'   x := A'*x.
     *
     *           Unchanged on exit.
     *
     *  DIAG   - CHARACTER*1.
     *           On entry, DIAG specifies whether or not A is unit
     *           triangular as follows:
     *
     *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
     *
     *              DIAG = 'N' or 'n'   A is not assumed to be unit
     *                                  triangular.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  AP     - DOUBLE PRECISION array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     *           Before entry with  UPLO = 'U' or 'u', the array AP must
     *           contain the upper triangular matrix packed sequentially,
     *           column by column, so that AP( 1 ) contains a( 1, 1 ),
     *           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
     *           respectively, and so on.
     *           Before entry with UPLO = 'L' or 'l', the array AP must
     *           contain the lower triangular matrix packed sequentially,
     *           column by column, so that AP( 1 ) contains a( 1, 1 ),
     *           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
     *           respectively, and so on.
     *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
     *           A are not referenced, but are assumed to be unity.
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element vector x. On exit, X is overwritten with the
     *           tranformed vector x.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param trans
     * @param diag
     * @param n
     * @param ap
     * @param x
     * @param incx
     */
    public final void dtpmv(String uplo, String trans, String diag, int n, double[] ap, double[] x, int incx) {
        dtpmv(uplo, trans, diag, n, ap, 0, x, 0, incx);
    }

    public abstract void dtpmv(String uplo, String trans, String diag, int n, double[] ap, int apOffset, double[] x,
            int xOffset, int incx);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DTRMM  performs one of the matrix-matrix operations
     *
     *     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
     *
     *  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
     *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     *
     *     op( A ) = A   or   op( A ) = A'.
     *
     *  Arguments
     *  ==========
     *
     *  SIDE   - CHARACTER*1.
     *           On entry,  SIDE specifies whether  op( A ) multiplies B from
     *           the left or right as follows:
     *
     *              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
     *
     *              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
     *
     *           Unchanged on exit.
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the matrix A is an upper or
     *           lower triangular matrix as follows:
     *
     *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
     *
     *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
     *
     *           Unchanged on exit.
     *
     *  TRANSA - CHARACTER*1.
     *           On entry, TRANSA specifies the form of op( A ) to be used in
     *           the matrix multiplication as follows:
     *
     *              TRANSA = 'N' or 'n'   op( A ) = A.
     *
     *              TRANSA = 'T' or 't'   op( A ) = A'.
     *
     *              TRANSA = 'C' or 'c'   op( A ) = A'.
     *
     *           Unchanged on exit.
     *
     *  DIAG   - CHARACTER*1.
     *           On entry, DIAG specifies whether or not A is unit triangular
     *           as follows:
     *
     *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
     *
     *              DIAG = 'N' or 'n'   A is not assumed to be unit
     *                                  triangular.
     *
     *           Unchanged on exit.
     *
     *  M      - INTEGER.
     *           On entry, M specifies the number of rows of B. M must be at
     *           least zero.
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the number of columns of B.  N must be
     *           at least zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - DOUBLE PRECISION.
     *           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
     *           zero then  A is not referenced and  B need not be set before
     *           entry.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
     *           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
     *           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
     *           upper triangular part of the array  A must contain the upper
     *           triangular matrix  and the strictly lower triangular part of
     *           A is not referenced.
     *           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
     *           lower triangular part of the array  A must contain the lower
     *           triangular matrix  and the strictly upper triangular part of
     *           A is not referenced.
     *           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
     *           A  are not referenced either,  but are assumed to be  unity.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
     *           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
     *           then LDA must be at least max( 1, n ).
     *           Unchanged on exit.
     *
     *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
     *           Before entry,  the leading  m by n part of the array  B must
     *           contain the matrix  B,  and  on exit  is overwritten  by the
     *           transformed matrix.
     *
     *  LDB    - INTEGER.
     *           On entry, LDB specifies the first dimension of B as declared
     *           in  the  calling  (sub)  program.   LDB  must  be  at  least
     *           max( 1, m ).
     *           Unchanged on exit.
     *
     *
     *  Level 3 Blas routine.
     *
     *  -- Written on 8-February-1989.
     *     Jack Dongarra, Argonne National Laboratory.
     *     Iain Duff, AERE Harwell.
     *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
     *     Sven Hammarling, Numerical Algorithms Group Ltd.
     *
     * </code>
     * </pre>
     *
     * @param side
     * @param uplo
     * @param transa
     * @param diag
     * @param m
     * @param n
     * @param alpha
     * @param a
     * @param lda
     * @param b
     * @param ldb
     */
    public final void dtrmm(String side, String uplo, String transa, String diag, int m, int n, double alpha,
            double[] a, int lda, double[] b, int ldb) {
        dtrmm(side, uplo, transa, diag, m, n, alpha, a, 0, lda, b, 0, ldb);
    }

    public abstract void dtrmm(String side, String uplo, String transa, String diag, int m, int n, double alpha,
            double[] a, int aOffset, int lda, double[] b, int bOffset, int ldb);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DTRMV  performs one of the matrix-vector operations
     *
     *     x := A*x,   or   x := A'*x,
     *
     *  where x is an n element vector and  A is an n by n unit, or non-unit,
     *  upper or lower triangular matrix.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the matrix is an upper or
     *           lower triangular matrix as follows:
     *
     *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
     *
     *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
     *
     *           Unchanged on exit.
     *
     *  TRANS  - CHARACTER*1.
     *           On entry, TRANS specifies the operation to be performed as
     *           follows:
     *
     *              TRANS = 'N' or 'n'   x := A*x.
     *
     *              TRANS = 'T' or 't'   x := A'*x.
     *
     *              TRANS = 'C' or 'c'   x := A'*x.
     *
     *           Unchanged on exit.
     *
     *  DIAG   - CHARACTER*1.
     *           On entry, DIAG specifies whether or not A is unit
     *           triangular as follows:
     *
     *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
     *
     *              DIAG = 'N' or 'n'   A is not assumed to be unit
     *                                  triangular.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     *           triangular matrix and the strictly lower triangular part of
     *           A is not referenced.
     *           Before entry with UPLO = 'L' or 'l', the leading n by n
     *           lower triangular part of the array A must contain the lower
     *           triangular matrix and the strictly upper triangular part of
     *           A is not referenced.
     *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
     *           A are not referenced either, but are assumed to be unity.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. LDA must be at least
     *           max( 1, n ).
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element vector x. On exit, X is overwritten with the
     *           tranformed vector x.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param trans
     * @param diag
     * @param n
     * @param a
     * @param lda
     * @param x
     * @param incx
     */
    public final void dtrmv(String uplo, String trans, String diag, int n, double[] a, int lda, double[] x, int incx) {
        dtrmv(uplo, trans, diag, n, a, 0, lda, x, 0, incx);
    }

    public abstract void dtrmv(String uplo, String trans, String diag, int n, double[] a, int aOffset, int lda,
            double[] x, int xOffset, int incx);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DTBSV  solves one of the systems of equations
     *
     *     A*x = b,   or   A'*x = b,
     *
     *  where b and x are n element vectors and A is an n by n unit, or
     *  non-unit, upper or lower triangular band matrix, with ( k + 1 )
     *  diagonals.
     *
     *  No test for singularity or near-singularity is included in this
     *  routine. Such tests must be performed before calling this routine.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the matrix is an upper or
     *           lower triangular matrix as follows:
     *
     *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
     *
     *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
     *
     *           Unchanged on exit.
     *
     *  TRANS  - CHARACTER*1.
     *           On entry, TRANS specifies the equations to be solved as
     *           follows:
     *
     *              TRANS = 'N' or 'n'   A*x = b.
     *
     *              TRANS = 'T' or 't'   A'*x = b.
     *
     *              TRANS = 'C' or 'c'   A'*x = b.
     *
     *           Unchanged on exit.
     *
     *  DIAG   - CHARACTER*1.
     *           On entry, DIAG specifies whether or not A is unit
     *           triangular as follows:
     *
     *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
     *
     *              DIAG = 'N' or 'n'   A is not assumed to be unit
     *                                  triangular.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  K      - INTEGER.
     *           On entry with UPLO = 'U' or 'u', K specifies the number of
     *           super-diagonals of the matrix A.
     *           On entry with UPLO = 'L' or 'l', K specifies the number of
     *           sub-diagonals of the matrix A.
     *           K must satisfy  0 .le. K.
     *           Unchanged on exit.
     *
     *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
     *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
     *           by n part of the array A must contain the upper triangular
     *           band part of the matrix of coefficients, supplied column by
     *           column, with the leading diagonal of the matrix in row
     *           ( k + 1 ) of the array, the first super-diagonal starting at
     *           position 2 in row k, and so on. The top left k by k triangle
     *           of the array A is not referenced.
     *           The following program segment will transfer an upper
     *           triangular band matrix from conventional full matrix storage
     *           to band storage:
     *
     *                 DO 20, J = 1, N
     *                    M = K + 1 - J
     *                    DO 10, I = MAX( 1, J - K ), J
     *                       A( M + I, J ) = matrix( I, J )
     *              10    CONTINUE
     *              20 CONTINUE
     *
     *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
     *           by n part of the array A must contain the lower triangular
     *           band part of the matrix of coefficients, supplied column by
     *           column, with the leading diagonal of the matrix in row 1 of
     *           the array, the first sub-diagonal starting at position 1 in
     *           row 2, and so on. The bottom right k by k triangle of the
     *           array A is not referenced.
     *           The following program segment will transfer a lower
     *           triangular band matrix from conventional full matrix storage
     *           to band storage:
     *
     *                 DO 20, J = 1, N
     *                    M = 1 - J
     *                    DO 10, I = J, MIN( N, J + K )
     *                       A( M + I, J ) = matrix( I, J )
     *              10    CONTINUE
     *              20 CONTINUE
     *
     *           Note that when DIAG = 'U' or 'u' the elements of the array A
     *           corresponding to the diagonal elements of the matrix are not
     *           referenced, but are assumed to be unity.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. LDA must be at least
     *           ( k + 1 ).
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element right-hand side vector b. On exit, X is overwritten
     *           with the solution vector x.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param trans
     * @param diag
     * @param n
     * @param k
     * @param a
     * @param lda
     * @param x
     * @param incx
     */
    public final void dtbsv(String uplo, String trans, String diag, int n, int k, double[] a, int lda, double[] x,
            int incx) {
        dtbsv(uplo, trans, diag, n, k, a, 0, lda, x, 0, incx);
    }

    public abstract void dtbsv(String uplo, String trans, String diag, int n, int k, double[] a, int aOffset, int lda,
            double[] x, int xOffset, int incx);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DTPSV  solves one of the systems of equations
     *
     *     A*x = b,   or   A'*x = b,
     *
     *  where b and x are n element vectors and A is an n by n unit, or
     *  non-unit, upper or lower triangular matrix, supplied in packed form.
     *
     *  No test for singularity or near-singularity is included in this
     *  routine. Such tests must be performed before calling this routine.
     *
     *  Arguments
     *  ==========
     *
     *  UPLO   - CHARACTER*1.
     *           On entry, UPLO specifies whether the matrix is an upper or
     *           lower triangular matrix as follows:
     *
     *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
     *
     *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
     *
     *           Unchanged on exit.
     *
     *  TRANS  - CHARACTER*1.
     *           On entry, TRANS specifies the equations to be solved as
     *           follows:
     *
     *              TRANS = 'N' or 'n'   A*x = b.
     *
     *              TRANS = 'T' or 't'   A'*x = b.
     *
     *              TRANS = 'C' or 'c'   A'*x = b.
     *
     *           Unchanged on exit.
     *
     *  DIAG   - CHARACTER*1.
     *           On entry, DIAG specifies whether or not A is unit
     *           triangular as follows:
     *
     *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
     *
     *              DIAG = 'N' or 'n'   A is not assumed to be unit
     *                                  triangular.
     *
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry, N specifies the order of the matrix A.
     *           N must be at least zero.
     *           Unchanged on exit.
     *
     *  AP     - DOUBLE PRECISION array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     *           Before entry with  UPLO = 'U' or 'u', the array AP must
     *           contain the upper triangular matrix packed sequentially,
     *           column by column, so that AP( 1 ) contains a( 1, 1 ),
     *           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
     *           respectively, and so on.
     *           Before entry with UPLO = 'L' or 'l', the array AP must
     *           contain the lower triangular matrix packed sequentially,
     *           column by column, so that AP( 1 ) contains a( 1, 1 ),
     *           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
     *           respectively, and so on.
     *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
     *           A are not referenced, but are assumed to be unity.
     *           Unchanged on exit.
     *
     *  X      - DOUBLE PRECISION array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *           Before entry, the incremented array X must contain the n
     *           element right-hand side vector b. On exit, X is overwritten
     *           with the solution vector x.
     *
     *  INCX   - INTEGER.
     *           On entry, INCX specifies the increment for the elements of
     *           X. INCX must not be zero.
     *           Unchanged on exit.
     *
     *
     *  Level 2 Blas routine.
     *
     *  -- Written on 22-October-1986.
     *     Jack Dongarra, Argonne National Lab.
     *     Jeremy Du Croz, Nag Central Office.
     *     Sven Hammarling, Nag Central Office.
     *     Richard Hanson, Sandia National Labs.
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param trans
     * @param diag
     * @param n
     * @param ap
     * @param x
     * @param incx
     */
    public final void dtpsv(String uplo, String trans, String diag, int n, double[] ap, double[] x, int incx) {
        dtpsv(uplo, trans, diag, n, ap, 0, x, 0, incx);
    }

    public abstract void dtpsv(String uplo, String trans, String diag, int n, double[] ap, int apOffset, double[] x,
            int xOffset, int incx);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *     constant times a vector plus a vector.
     *     uses unrolled loops for increments equal to one.
     *     jack dongarra, linpack, 3/11/78.
     *     modified 12/3/93, array(1) declarations changed to array(*)
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param da
     * @param dx
     * @param incx
     * @param dy
     * @param incy
     */
    public final void daxpy(int n, double da, double[] dx, int incx, double[] dy, int incy) {
        daxpy(n, da, dx, 0, incx, dy, 0, incy);
    }

    public abstract void daxpy(int n, double da, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset,
            int incy);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *     copies a vector, x, to a vector, y.
     *     uses unrolled loops for increments equal to one.
     *     jack dongarra, linpack, 3/11/78.
     *     modified 12/3/93, array(1) declarations changed to array(*)
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param dx
     * @param incx
     * @param dy
     * @param incy
     */
    public final void dcopy(int n, double[] dx, int incx, double[] dy, int incy) {
        dcopy(n, dx, 0, incx, dy, 0, incy);
    }

    public abstract void dcopy(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     **
     *     scales a vector by a constant.
     *     uses unrolled loops for increment equal to one.
     *     jack dongarra, linpack, 3/11/78.
     *     modified 3/93 to return if incx .le. 0.
     *     modified 12/3/93, array(1) declarations changed to array(*)
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param da
     * @param dx
     * @param incx
     */
    public final void dscal(int n, double da, double[] dx, int incx) {
        dscal(n, da, dx, 0, incx);
    }

    public abstract void dscal(int n, double da, double[] dx, int dxOffset, int incx);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *     interchanges two vectors.
     *     uses unrolled loops for increments equal one.
     *     jack dongarra, linpack, 3/11/78.
     *     modified 12/3/93, array(1) declarations changed to array(*)
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param dx
     * @param incx
     * @param dy
     * @param incy
     */
    public final void dswap(int n, double[] dx, int incx, double[] dy, int incy) {
        dswap(n, dx, 0, incx, dy, 0, incy);
    }

    public abstract void dswap(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *     forms the dot product of two vectors.
     *     uses unrolled loops for increments equal to one.
     *     jack dongarra, linpack, 3/11/78.
     *     modified 12/3/93, array(1) declarations changed to array(*)
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param dx
     * @param incx
     * @param dy
     * @param incy
     * @return the dot product of {@code dx} and {@code dy}
     */
    public final double ddot(int n, double[] dx, int incx, double[] dy, int incy) {
        return ddot(n, dx, 0, incx, dy, 0, incy);
    }

    public abstract double ddot(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *     applies a plane rotation.
     *     jack dongarra, linpack, 3/11/78.
     *     modified 12/3/93, array(1) declarations changed to array(*)
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param dx
     * @param incx
     * @param dy
     * @param incy
     * @param c
     * @param s
     */
    public final void drot(int n, double[] dx, int incx, double[] dy, int incy, double c, double s) {
        drot(n, dx, 0, incx, dy, 0, incy, c, s);
    }

    public abstract void drot(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy, double c,
            double s);

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *     finds the index of element having max. absolute value.
     *     jack dongarra, linpack, 3/11/78.
     *     modified 3/93 to return if incx .le. 0.
     *     modified 12/3/93, array(1) declarations changed to array(*)
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param dx
     * @param incx
     * @return index of the element having max. absolute value
     */
    public final int idamax(int n, double[] dx, int incx) {
        return idamax(n, dx, 0, incx);
    }

    public abstract int idamax(int n, double[] dx, int dxOffset, int incx);

    // miscellaneous float routines

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  SGEMM  performs one of the matrix-matrix operations
     *
     *     C := alpha*op( A )*op( B ) + beta*C,
     *
     *  where  op( X ) is one of
     *
     *     op( X ) = X   or   op( X ) = X',
     *
     *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
     *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
     *
     *  Arguments
     *  ==========
     *
     *  TRANSA - CHARACTER*1.
     *           On entry, TRANSA specifies the form of op( A ) to be used in
     *           the matrix multiplication as follows:
     *
     *              TRANSA = 'N' or 'n',  op( A ) = A.
     *
     *              TRANSA = 'T' or 't',  op( A ) = A'.
     *
     *              TRANSA = 'C' or 'c',  op( A ) = A'.
     *
     *           Unchanged on exit.
     *
     *  TRANSB - CHARACTER*1.
     *           On entry, TRANSB specifies the form of op( B ) to be used in
     *           the matrix multiplication as follows:
     *
     *              TRANSB = 'N' or 'n',  op( B ) = B.
     *
     *              TRANSB = 'T' or 't',  op( B ) = B'.
     *
     *              TRANSB = 'C' or 'c',  op( B ) = B'.
     *
     *           Unchanged on exit.
     *
     *  M      - INTEGER.
     *           On entry,  M  specifies  the number  of rows  of the  matrix
     *           op( A )  and of the  matrix  C.  M  must  be at least  zero.
     *           Unchanged on exit.
     *
     *  N      - INTEGER.
     *           On entry,  N  specifies the number  of columns of the matrix
     *           op( B ) and the number of columns of the matrix C. N must be
     *           at least zero.
     *           Unchanged on exit.
     *
     *  K      - INTEGER.
     *           On entry,  K  specifies  the number of columns of the matrix
     *           op( A ) and the number of rows of the matrix op( B ). K must
     *           be at least  zero.
     *           Unchanged on exit.
     *
     *  ALPHA  - REAL            .
     *           On entry, ALPHA specifies the scalar alpha.
     *           Unchanged on exit.
     *
     *  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
     *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by m  part of the array  A  must contain  the
     *           matrix A.
     *           Unchanged on exit.
     *
     *  LDA    - INTEGER.
     *           On entry, LDA specifies the first dimension of A as declared
     *           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
     *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
     *           least  max( 1, k ).
     *           Unchanged on exit.
     *
     *  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
     *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
     *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  n by k  part of the array  B  must contain  the
     *           matrix B.
     *           Unchanged on exit.
     *
     *  LDB    - INTEGER.
     *           On entry, LDB specifies the first dimension of B as declared
     *           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
     *           LDB must be at least  max( 1, k ), otherwise  LDB must be at
     *           least  max( 1, n ).
     *           Unchanged on exit.
     *
     *  BETA   - REAL            .
     *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
     *           supplied as zero then C need not be set on input.
     *           Unchanged on exit.
     *
     *  C      - REAL             array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     *           On exit, the array  C  is overwritten by the  m by n  matrix
     *           ( alpha*op( A )*op( B ) + beta*C ).
     *
     *  LDC    - INTEGER.
     *           On entry, LDC specifies the first dimension of C as declared
     *           in  the  calling  (sub)  program.   LDC  must  be  at  least
     *           max( 1, m ).
     *           Unchanged on exit.
     *
     *
     *  Level 3 Blas routine.
     *
     *  -- Written on 8-February-1989.
     *     Jack Dongarra, Argonne National Laboratory.
     *     Iain Duff, AERE Harwell.
     *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
     *     Sven Hammarling, Numerical Algorithms Group Ltd.
     *
     * </code>
     * </pre>
     *
     * @param transa
     * @param transb
     * @param m
     * @param n
     * @param k
     * @param alpha
     * @param a
     * @param lda
     * @param b
     * @param ldb
     * @param beta
     * @param c
     * @param ldc
     */
    public final void sgemm(String transa, String transb, int m, int n, int k, float alpha, float[] a, int lda,
            float[] b, int ldb, float beta, float[] c, int ldc) {
        sgemm(transa, transb, m, n, k, alpha, a, 0, lda, b, 0, ldb, beta, c, 0, ldc);
    }

    public abstract void sgemm(String transa, String transb, int m, int n, int k, float alpha, float[] a, int aOffset,
            int lda, float[] b, int bOffset, int ldb, float beta, float[] c, int cOffset, int ldc);

    // miscellaneous complex routines

    // Note:
    // An array of 3 complex numbers a0, a1 and a2 is always represented as an
    // array of 6 real numbers [a0.re, a0.im, a1.re, a1.im, a2.re, a2.im]

    public abstract void cgemm(Trans transa, Trans transb, int m, int n, int k, float alphar, float alphai, float[] a,
            int lda, float[] b, int ldb, float betar, float betai, float[] c, int ldc);

    public abstract void zgemm(Trans transa, Trans transb, int m, int n, int k, double alphar, double alphai,
            double[] a, int lda, double[] b, int ldb, double betar, double betai, double[] c, int ldc);

    protected Blas() {
    }
}
