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
package net.frobenius.lapack;

import org.netlib.util.doubleW;
import org.netlib.util.intW;

import net.dedekind.lapack.Lapack;
import net.frobenius.ComputationTruncatedException;
import net.frobenius.NotConvergedException;
import net.frobenius.TDiag;
import net.frobenius.TEigJob;
import net.frobenius.TNorm;
import net.frobenius.TRange;
import net.frobenius.TSide;
import net.frobenius.TSvdJob;
import net.frobenius.TTrans;
import net.frobenius.TUpLo;

/**
 * Matrix storage layout must be column-major as in Fortran. The number of
 * matrix rows and columns, if required, must be strictly positive. All
 * operations throw a {@code NullPointerException} if any of the reference
 * method arguments is {@code null}.
 */
public final class PlainLapack {
    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGBCON estimates the reciprocal of the condition number of a real
     *  general band matrix A, in either the 1-norm or the infinity-norm,
     *  using the LU factorization computed by DGBTRF.
     *
     *  An estimate is obtained for norm(inv(A)), and the reciprocal of the
     *  condition number is computed as
     *     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
     *
     *  Arguments
     *  =========
     *
     *  NORM    (input) CHARACTER*1
     *          Specifies whether the 1-norm condition number or the
     *          infinity-norm condition number is required:
     *          = '1' or 'O':  1-norm;
     *          = 'I':         Infinity-norm.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  KL      (input) INTEGER
     *          The number of subdiagonals within the band of A.  KL >= 0.
     *
     *  KU      (input) INTEGER
     *          The number of superdiagonals within the band of A.  KU >= 0.
     *
     *  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
     *          Details of the LU factorization of the band matrix A, as
     *          computed by DGBTRF.  U is stored as an upper triangular band
     *          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
     *          the multipliers used during the factorization are stored in
     *          rows KL+KU+2 to 2*KL+KU+1.
     *
     *  LDAB    (input) INTEGER
     *          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
     *
     *  IPIV    (input) INTEGER array, dimension (N)
     *          The pivot indices; for 1 <= i <= N, row i of the matrix was
     *          interchanged with row IPIV(i).
     *
     *  ANORM   (input) DOUBLE PRECISION
     *          If NORM = '1' or 'O', the 1-norm of the original matrix A.
     *          If NORM = 'I', the infinity-norm of the original matrix A.
     *
     *  RCOND   (output) DOUBLE PRECISION
     *          The reciprocal of the condition number of the matrix A,
     *          computed as RCOND = 1/(norm(A) * norm(inv(A))).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param norm
     * @param n
     * @param kl
     * @param ku
     * @param ab
     * @param indices
     * @param normA
     */
    public static double dgbcon(Lapack la, TNorm norm, int n, int kl, int ku, double[] ab, int[] indices,
            double normA) {
        checkStrictlyPositive(n, "n");
        checkNonNegative(kl, "kl");
        checkNonNegative(ku, "ku");
        checkMinLen(indices, n, "indices");
        int ldab = Math.max(1, 2 * kl + ku + 1);
        checkMinLen(ab, ldab * n, "ab");

        intW info = new intW(0);
        doubleW rcond = new doubleW(0.0);
        double[] work = new double[3 * n];
        int[] iwork = new int[n];
        la.dgbcon(norm.val(), n, kl, ku, ab, ldab, indices, normA, rcond, work, iwork, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        return rcond.val;
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGBSV computes the solution to a real system of linear equations
     *  A * X = B, where A is a band matrix of order N with KL subdiagonals
     *  and KU superdiagonals, and X and B are N-by-NRHS matrices.
     *
     *  The LU decomposition with partial pivoting and row interchanges is
     *  used to factor A as A = L * U, where L is a product of permutation
     *  and unit lower triangular matrices with KL subdiagonals, and U is
     *  upper triangular with KL+KU superdiagonals.  The factored form of A
     *  is then used to solve the system of equations A * X = B.
     *
     *  Arguments
     *  =========
     *
     *  N       (input) INTEGER
     *          The number of linear equations, i.e., the order of the
     *          matrix A.  N >= 0.
     *
     *  KL      (input) INTEGER
     *          The number of subdiagonals within the band of A.  KL >= 0.
     *
     *  KU      (input) INTEGER
     *          The number of superdiagonals within the band of A.  KU >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
     *          On entry, the matrix A in band storage, in rows KL+1 to
     *          2*KL+KU+1; rows 1 to KL of the array need not be set.
     *          The j-th column of A is stored in the j-th column of the
     *          array AB as follows:
     *          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
     *          On exit, details of the factorization: U is stored as an
     *          upper triangular band matrix with KL+KU superdiagonals in
     *          rows 1 to KL+KU+1, and the multipliers used during the
     *          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
     *          See below for further details.
     *
     *  LDAB    (input) INTEGER
     *          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
     *
     *  IPIV    (output) INTEGER array, dimension (N)
     *          The pivot indices that define the permutation matrix P;
     *          row i of the matrix was interchanged with row IPIV(i).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the N-by-NRHS right hand side matrix B.
     *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
     *                has been completed, but the factor U is exactly
     *                singular, and the solution has not been computed.
     *
     *  Further Details
     *  ===============
     *
     *  The band storage scheme is illustrated by the following example, when
     *  M = N = 6, KL = 2, KU = 1:
     *
     *  On entry:                       On exit:
     *
     *      *    *    *    +    +    +       *    *    *   u14  u25  u36
     *      *    *    +    +    +    +       *    *   u13  u24  u35  u46
     *      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
     *     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
     *     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
     *     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
     *
     *  Array elements marked * are not used by the routine; elements marked
     *  + need not be set on entry, but are required by the routine to store
     *  elements of U because of fill-in resulting from the row interchanges.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param kl
     * @param ku
     * @param rhsCount
     * @param ab
     * @param indices
     * @param b
     * @param ldb
     */
    public static void dgbsv(Lapack la, int n, int kl, int ku, int rhsCount, double[] ab, int[] indices, double[] b,
            int ldb) {
        checkStrictlyPositive(n, "n");
        checkNonNegative(kl, "kl");
        checkNonNegative(ku, "ku");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkMinLen(indices, n, "indices");
        int ldab = Math.max(1, 2 * kl + ku + 1);
        checkMinLen(ab, ldab * n, "ab");
        checkValueAtLeast(ldb, n, "ldb");
        checkMinLen(b, ldb * rhsCount, "b");

        intW info = new intW(0);
        la.dgbsv(n, kl, ku, rhsCount, ab, ldab, indices, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "Factor U in the LU decomposition is exactly singular. Solution could not be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGBTRF computes an LU factorization of a real m-by-n band matrix A
     *  using partial pivoting with row interchanges.
     *
     *  This is the blocked version of the algorithm, calling Level 3 BLAS.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     *
     *  KL      (input) INTEGER
     *          The number of subdiagonals within the band of A.  KL >= 0.
     *
     *  KU      (input) INTEGER
     *          The number of superdiagonals within the band of A.  KU >= 0.
     *
     *  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
     *          On entry, the matrix A in band storage, in rows KL+1 to
     *          2*KL+KU+1; rows 1 to KL of the array need not be set.
     *          The j-th column of A is stored in the j-th column of the
     *          array AB as follows:
     *          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
     *
     *          On exit, details of the factorization: U is stored as an
     *          upper triangular band matrix with KL+KU superdiagonals in
     *          rows 1 to KL+KU+1, and the multipliers used during the
     *          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
     *          See below for further details.
     *
     *  LDAB    (input) INTEGER
     *          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
     *
     *  IPIV    (output) INTEGER array, dimension (min(M,N))
     *          The pivot indices; for 1 <= i <= min(M,N), row i of the
     *          matrix was interchanged with row IPIV(i).
     *
     *  INFO    (output) INTEGER
     *          = 0: successful exit
     *          < 0: if INFO = -i, the i-th argument had an illegal value
     *          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
     *               has been completed, but the factor U is exactly
     *               singular, and division by zero will occur if it is used
     *               to solve a system of equations.
     *
     *  Further Details
     *  ===============
     *
     *  The band storage scheme is illustrated by the following example, when
     *  M = N = 6, KL = 2, KU = 1:
     *
     *  On entry:                       On exit:
     *
     *      *    *    *    +    +    +       *    *    *   u14  u25  u36
     *      *    *    +    +    +    +       *    *   u13  u24  u35  u46
     *      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
     *     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
     *     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
     *     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
     *
     *  Array elements marked * are not used by the routine; elements marked
     *  + need not be set on entry, but are required by the routine to store
     *  elements of U because of fill-in resulting from the row interchanges.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param kl
     * @param ku
     * @param ab
     * @param indices
     */
    public static void dgbtrf(Lapack la, int m, int n, int kl, int ku, double[] ab, int[] indices) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkNonNegative(kl, "kl");
        checkNonNegative(ku, "ku");
        checkMinLen(indices, Math.min(m, n), "indices");
        int ldab = Math.max(1, 2 * kl + ku + 1);
        checkMinLen(ab, ldab * n, "ab");

        intW info = new intW(0);
        la.dgbtrf(m, n, kl, ku, ab, ldab, indices, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("Factor U in the LU decomposition is exactly singular");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGBTRS solves a system of linear equations
     *     A * X = B  or  A' * X = B
     *  with a general band matrix A using the LU factorization computed
     *  by DGBTRF.
     *
     *  Arguments
     *  =========
     *
     *  TRANS   (input) CHARACTER*1
     *          Specifies the form of the system of equations.
     *          = 'N':  A * X = B  (No transpose)
     *          = 'T':  A'* X = B  (Transpose)
     *          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  KL      (input) INTEGER
     *          The number of subdiagonals within the band of A.  KL >= 0.
     *
     *  KU      (input) INTEGER
     *          The number of superdiagonals within the band of A.  KU >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
     *          Details of the LU factorization of the band matrix A, as
     *          computed by DGBTRF.  U is stored as an upper triangular band
     *          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
     *          the multipliers used during the factorization are stored in
     *          rows KL+KU+2 to 2*KL+KU+1.
     *
     *  LDAB    (input) INTEGER
     *          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
     *
     *  IPIV    (input) INTEGER array, dimension (N)
     *          The pivot indices; for 1 <= i <= N, row i of the matrix was
     *          interchanged with row IPIV(i).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the right hand side matrix B.
     *          On exit, the solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param trans
     * @param n
     * @param kl
     * @param ku
     * @param rhsCount
     * @param ab
     * @param indices
     * @param b
     * @param ldb
     */
    public static void dgbtrs(Lapack la, TTrans trans, int n, int kl, int ku, int rhsCount, double[] ab, int[] indices,
            double[] b, int ldb) {
        checkStrictlyPositive(n, "n");
        checkNonNegative(kl, "kl");
        checkNonNegative(ku, "ku");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkMinLen(indices, n, "indices");
        int ldab = Math.max(1, 2 * kl + ku + 1);
        checkMinLen(ab, ldab * n, "ab");
        checkValueAtLeast(ldb, n, "ldb");
        checkMinLen(b, ldb * rhsCount, "b");

        intW info = new intW(0);
        la.dgbtrs(trans.val(), n, kl, ku, rhsCount, ab, ldab, indices, b, ldb, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGECON estimates the reciprocal of the condition number of a general
     *  real matrix A, in either the 1-norm or the infinity-norm, using
     *  the LU factorization computed by DGETRF.
     *
     *  An estimate is obtained for norm(inv(A)), and the reciprocal of the
     *  condition number is computed as
     *     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
     *
     *  Arguments
     *  =========
     *
     *  NORM    (input) CHARACTER*1
     *          Specifies whether the 1-norm condition number or the
     *          infinity-norm condition number is required:
     *          = '1' or 'O':  1-norm;
     *          = 'I':         Infinity-norm.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
     *          The factors L and U from the factorization A = P*L*U
     *          as computed by DGETRF.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  ANORM   (input) DOUBLE PRECISION
     *          If NORM = '1' or 'O', the 1-norm of the original matrix A.
     *          If NORM = 'I', the infinity-norm of the original matrix A.
     *
     *  RCOND   (output) DOUBLE PRECISION
     *          The reciprocal of the condition number of the matrix A,
     *          computed as RCOND = 1/(norm(A) * norm(inv(A))).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param norm
     * @param n
     * @param a
     * @param lda
     * @param normA
     */
    public static double dgecon(Lapack la, TNorm norm, int n, double[] a, int lda, double normA) {
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, n, "lda");
        checkMinLen(a, lda * n, "a");

        intW info = new intW(0);
        doubleW rcond = new doubleW(0.0);
        double[] work = new double[4 * n];
        int[] iwork = new int[n];
        la.dgecon(norm.val(), n, a, lda, normA, rcond, work, iwork, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        return rcond.val;
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGEEV computes for an N-by-N real nonsymmetric matrix A, the
     *  eigenvalues and, optionally, the left and/or right eigenvectors.
     *
     *  The right eigenvector v(j) of A satisfies
     *                   A * v(j) = lambda(j) * v(j)
     *  where lambda(j) is its eigenvalue.
     *  The left eigenvector u(j) of A satisfies
     *                u(j)**H * A = lambda(j) * u(j)**H
     *  where u(j)**H denotes the conjugate transpose of u(j).
     *
     *  The computed eigenvectors are normalized to have Euclidean norm
     *  equal to 1 and largest component real.
     *
     *  Arguments
     *  =========
     *
     *  JOBVL   (input) CHARACTER*1
     *          = 'N': left eigenvectors of A are not computed;
     *          = 'V': left eigenvectors of A are computed.
     *
     *  JOBVR   (input) CHARACTER*1
     *          = 'N': right eigenvectors of A are not computed;
     *          = 'V': right eigenvectors of A are computed.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A. N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the N-by-N matrix A.
     *          On exit, A has been overwritten.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  WR      (output) DOUBLE PRECISION array, dimension (N)
     *  WI      (output) DOUBLE PRECISION array, dimension (N)
     *          WR and WI contain the real and imaginary parts,
     *          respectively, of the computed eigenvalues.  Complex
     *          conjugate pairs of eigenvalues appear consecutively
     *          with the eigenvalue having the positive imaginary part
     *          first.
     *
     *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
     *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
     *          after another in the columns of VL, in the same order
     *          as their eigenvalues.
     *          If JOBVL = 'N', VL is not referenced.
     *          If the j-th eigenvalue is real, then u(j) = VL(:,j),
     *          the j-th column of VL.
     *          If the j-th and (j+1)-st eigenvalues form a complex
     *          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
     *          u(j+1) = VL(:,j) - i*VL(:,j+1).
     *
     *  LDVL    (input) INTEGER
     *          The leading dimension of the array VL.  LDVL >= 1; if
     *          JOBVL = 'V', LDVL >= N.
     *
     *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
     *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
     *          after another in the columns of VR, in the same order
     *          as their eigenvalues.
     *          If JOBVR = 'N', VR is not referenced.
     *          If the j-th eigenvalue is real, then v(j) = VR(:,j),
     *          the j-th column of VR.
     *          If the j-th and (j+1)-st eigenvalues form a complex
     *          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
     *          v(j+1) = VR(:,j) - i*VR(:,j+1).
     *
     *  LDVR    (input) INTEGER
     *          The leading dimension of the array VR.  LDVR >= 1; if
     *          JOBVR = 'V', LDVR >= N.
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value.
     *          > 0:  if INFO = i, the QR algorithm failed to compute all the
     *                eigenvalues, and no eigenvectors have been computed;
     *                elements i+1:N of WR and WI contain eigenvalues which
     *                have converged.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param jobvl
     * @param jobvr
     * @param n
     * @param a
     * @param lda
     * @param wr
     * @param wi
     * @param vl
     * @param ldvl
     * @param vr
     * @param ldvr
     */
    public static void dgeev(Lapack la, TEigJob jobvl, TEigJob jobvr, int n, double[] a, int lda, double[] wr,
            double[] wi, double[] vl, int ldvl, double[] vr, int ldvr) {
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, n, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(wr, n, "wr");
        checkMinLen(wi, n, "wi");
        if (jobvl == TEigJob.ALL) {
            checkValueAtLeast(ldvl, n, "ldvl");
            checkMinLen(vl, ldvl * n, "vl");
        } else {
            checkValueAtLeast(ldvl, 1, "ldvl");
        }
        if (jobvr == TEigJob.ALL) {
            checkValueAtLeast(ldvr, n, "ldvr");
            checkMinLen(vr, ldvr * n, "vr");
        } else {
            checkValueAtLeast(ldvr, 1, "ldvr");
        }

        intW info = new intW(0);
        double[] work = new double[1];
        la.dgeev(jobvl.val(), jobvr.val(), n, new double[0], lda, new double[0], new double[0], new double[0], ldvl,
                new double[0], ldvr, work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dgeev(jobvl.val(), jobvr.val(), n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, work.length, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "QR failed to compute all eigenvalues, eigenvectors haven't been computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  SGEEV computes for an N-by-N real nonsymmetric matrix A, the
     *  eigenvalues and, optionally, the left and/or right eigenvectors.
     *
     *  The right eigenvector v(j) of A satisfies
     *                   A * v(j) = lambda(j) * v(j)
     *  where lambda(j) is its eigenvalue.
     *  The left eigenvector u(j) of A satisfies
     *                u(j)**H * A = lambda(j) * u(j)**H
     *  where u(j)**H denotes the conjugate transpose of u(j).
     *
     *  The computed eigenvectors are normalized to have Euclidean norm
     *  equal to 1 and largest component real.
     *
     *  Arguments
     *  =========
     *
     *  JOBVL   (input) CHARACTER*1
     *          = 'N': left eigenvectors of A are not computed;
     *          = 'V': left eigenvectors of A are computed.
     *
     *  JOBVR   (input) CHARACTER*1
     *          = 'N': right eigenvectors of A are not computed;
     *          = 'V': right eigenvectors of A are computed.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A. N >= 0.
     *
     *  A       (input/output) REAL array, dimension (LDA,N)
     *          On entry, the N-by-N matrix A.
     *          On exit, A has been overwritten.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  WR      (output) REAL array, dimension (N)
     *  WI      (output) REAL array, dimension (N)
     *          WR and WI contain the real and imaginary parts,
     *          respectively, of the computed eigenvalues.  Complex
     *          conjugate pairs of eigenvalues appear consecutively
     *          with the eigenvalue having the positive imaginary part
     *          first.
     *
     *  VL      (output) REAL array, dimension (LDVL,N)
     *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
     *          after another in the columns of VL, in the same order
     *          as their eigenvalues.
     *          If JOBVL = 'N', VL is not referenced.
     *          If the j-th eigenvalue is real, then u(j) = VL(:,j),
     *          the j-th column of VL.
     *          If the j-th and (j+1)-st eigenvalues form a complex
     *          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
     *          u(j+1) = VL(:,j) - i*VL(:,j+1).
     *
     *  LDVL    (input) INTEGER
     *          The leading dimension of the array VL.  LDVL >= 1; if
     *          JOBVL = 'V', LDVL >= N.
     *
     *  VR      (output) REAL array, dimension (LDVR,N)
     *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
     *          after another in the columns of VR, in the same order
     *          as their eigenvalues.
     *          If JOBVR = 'N', VR is not referenced.
     *          If the j-th eigenvalue is real, then v(j) = VR(:,j),
     *          the j-th column of VR.
     *          If the j-th and (j+1)-st eigenvalues form a complex
     *          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
     *          v(j+1) = VR(:,j) - i*VR(:,j+1).
     *
     *  LDVR    (input) INTEGER
     *          The leading dimension of the array VR.  LDVR >= 1; if
     *          JOBVR = 'V', LDVR >= N.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param jobvl
     * @param jobvr
     * @param n
     * @param a
     * @param lda
     * @param wr
     * @param wi
     * @param vl
     * @param ldvl
     * @param vr
     * @param ldvr
     */
    public static void sgeev(Lapack la, TEigJob jobvl, TEigJob jobvr, int n, float[] a, int lda, float[] wr, float[] wi,
            float[] vl, int ldvl, float[] vr, int ldvr) {
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, n, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(wr, n, "wr");
        checkMinLen(wi, n, "wi");
        if (jobvl == TEigJob.ALL) {
            checkValueAtLeast(ldvl, n, "ldvl");
            checkMinLen(vl, ldvl * n, "vl");
        } else {
            checkValueAtLeast(ldvl, 1, "ldvl");
        }
        if (jobvr == TEigJob.ALL) {
            checkValueAtLeast(ldvr, n, "ldvr");
            checkMinLen(vr, ldvr * n, "vr");
        } else {
            checkValueAtLeast(ldvr, 1, "ldvr");
        }

        intW info = new intW(0);
        float[] work = new float[1];
        la.sgeev(jobvl.val(), jobvr.val(), n, new float[0], lda, new float[0], new float[0], new float[0], ldvl,
                new float[0], ldvr, work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new float[(int) work[0]];
        la.sgeev(jobvl.val(), jobvr.val(), n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, work.length, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "QR failed to compute all eigenvalues, eigenvectors haven't been computed.");
            }
        }
    }

    // note: there is no need to make a copy of 'a' before calling cgeev!
    public static void cgeev(Lapack la, TEigJob jobvl, TEigJob jobvr, int n, float[] a, int lda, float[] w, float[] vl,
            int ldvl, float[] vr, int ldvr) {
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, n, "lda");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(w, 2 * n, "w");
        if (jobvl == TEigJob.ALL) {
            checkValueAtLeast(ldvl, n, "ldvl");
            checkMinLen(vl, 2 * ldvl * n, "vl");
        } else {
            checkValueAtLeast(ldvl, 1, "ldvl");
        }
        if (jobvr == TEigJob.ALL) {
            checkValueAtLeast(ldvr, n, "ldvr");
            checkMinLen(vr, 2 * ldvr * n, "vr");
        } else {
            checkValueAtLeast(ldvr, 1, "ldvr");
        }

        int info = la.cgeev(jobvl.val(), jobvr.val(), n, a, lda, w, vl, ldvl, vr, ldvr);
        if (info != 0) {
            if (info < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "QR failed to compute all eigenvalues, eigenvectors haven't been computed.");
            }
        }
    }

    // note: there is no need to make a copy of 'a' before calling zgeev!
    public static void zgeev(Lapack la, TEigJob jobvl, TEigJob jobvr, int n, double[] a, int lda, double[] w,
            double[] vl, int ldvl, double[] vr, int ldvr) {
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, n, "lda");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(w, 2 * n, "w");
        if (jobvl == TEigJob.ALL) {
            checkValueAtLeast(ldvl, n, "ldvl");
            checkMinLen(vl, 2 * ldvl * n, "vl");
        } else {
            checkValueAtLeast(ldvl, 1, "ldvl");
        }
        if (jobvr == TEigJob.ALL) {
            checkValueAtLeast(ldvr, n, "ldvr");
            checkMinLen(vr, 2 * ldvr * n, "vr");
        } else {
            checkValueAtLeast(ldvr, 1, "ldvr");
        }

        int info = la.zgeev(jobvl.val(), jobvr.val(), n, a, lda, w, vl, ldvl, vr, ldvr);
        if (info != 0) {
            if (info < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "QR failed to compute all eigenvalues, eigenvectors haven't been computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGELQF computes an LQ factorization of a real M-by-N matrix A:
     *  A = L * Q.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit, the elements on and below the diagonal of the array
     *          contain the m-by-min(m,n) lower trapezoidal matrix L (L is
     *          lower triangular if m <= n); the elements above the diagonal,
     *          with the array TAU, represent the orthogonal matrix Q as a
     *          product of elementary reflectors (see Further Details).
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     *
     *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
     *          The scalar factors of the elementary reflectors (see Further
     *          Details).
    *
     *  Further Details
     *  ===============
     *
     *  The matrix Q is represented as a product of elementary reflectors
     *
     *     Q = H(k) . . . H(2) H(1), where k = min(m,n).
     *
     *  Each H(i) has the form
     *
     *     H(i) = I - tau * v * v'
     *
     *  where tau is a real scalar, and v is a real vector with
     *  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
     *  and tau in TAU(i).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param a
     * @param lda
     * @param tau
     */
    public static void dgelqf(Lapack la, int m, int n, double[] a, int lda, double[] tau) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "m");
        checkMinLen(a, lda * n, "a");
        checkMinLen(tau, Math.min(n, m), "tau");

        intW info = new intW(0);
        double[] work = new double[1];
        la.dgelqf(m, n, new double[0], lda, new double[0], work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dgelqf(m, n, a, lda, tau, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGELS solves overdetermined or underdetermined real linear systems
     *  involving an M-by-N matrix A, or its transpose, using a QR or LQ
     *  factorization of A.  It is assumed that A has full rank.
     *
     *  The following options are provided:
     *
     *  1. If TRANS = 'N' and m >= n:  find the least squares solution of
     *     an overdetermined system, i.e., solve the least squares problem
     *                  minimize || B - A*X ||.
     *
     *  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     *     an underdetermined system A * X = B.
     *
     *  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
     *     an undetermined system A**T * X = B.
     *
     *  4. If TRANS = 'T' and m < n:  find the least squares solution of
     *     an overdetermined system, i.e., solve the least squares problem
     *                  minimize || B - A**T * X ||.
     *
     *  Several right hand side vectors b and solution vectors x can be
     *  handled in a single call; they are stored as the columns of the
     *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     *  matrix X.
     *
     *  Arguments
     *  =========
     *
     *  TRANS   (input) CHARACTER*1
     *          = 'N': the linear system involves A;
     *          = 'T': the linear system involves A**T.
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of
     *          columns of the matrices B and X. NRHS >=0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit,
     *            if M >= N, A is overwritten by details of its QR
     *                       factorization as returned by DGEQRF;
     *            if M <  N, A is overwritten by details of its LQ
     *                       factorization as returned by DGELQF.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the matrix B of right hand side vectors, stored
     *          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
     *          if TRANS = 'T'.
     *          On exit, if INFO = 0, B is overwritten by the solution
     *          vectors, stored columnwise:
     *          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
     *          squares solution vectors; the residual sum of squares for the
     *          solution in each column is given by the sum of squares of
     *          elements N+1 to M in that column;
     *          if TRANS = 'N' and m < n, rows 1 to N of B contain the
     *          minimum norm solution vectors;
     *          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
     *          minimum norm solution vectors;
     *          if TRANS = 'T' and m < n, rows 1 to M of B contain the
     *          least squares solution vectors; the residual sum of squares
     *          for the solution in each column is given by the sum of
     *          squares of elements M+1 to N in that column.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B. LDB >= MAX(1,M,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO =  i, the i-th diagonal element of the
     *                triangular factor of A is zero, so that A does not have
     *                full rank; the least squares solution could not be
     *                computed.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param trans
     * @param m
     * @param n
     * @param rhsCount
     * @param a
     * @param lda
     * @param b
     * @param ldb
     */
    public static void dgels(Lapack la, TTrans trans, int m, int n, int rhsCount, double[] a, int lda, double[] b,
            int ldb) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkValueAtLeast(lda, m, "lda");
        checkValueAtLeast(ldb, Math.max(n, m), "ldb");
        checkMinLen(a, lda * n, "a");
        checkMinLen(b, ldb * rhsCount, "b");

        intW info = new intW(0);
        double[] work = new double[1];
        la.dgels(trans.val(), m, n, rhsCount, new double[0], lda, new double[0], ldb, work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dgels(trans.val(), m, n, rhsCount, a, lda, b, ldb, work, work.length, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "A does not have full rank. Least squares solution could not be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     * 
     *  Purpose
     *  =======
     * 
     *  SGELS solves overdetermined or underdetermined real linear systems
     *  involving an M-by-N matrix A, or its transpose, using a QR or LQ
     *  factorization of A.  It is assumed that A has full rank.
     * 
     *  The following options are provided:
     * 
     *  1. If TRANS = 'N' and m >= n:  find the least squares solution of
     *     an overdetermined system, i.e., solve the least squares problem
     *                  minimize || B - A*X ||.
     * 
     *  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
     *     an underdetermined system A * X = B.
     * 
     *  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
     *     an undetermined system A**T * X = B.
     * 
     *  4. If TRANS = 'T' and m < n:  find the least squares solution of
     *     an overdetermined system, i.e., solve the least squares problem
     *                  minimize || B - A**T * X ||.
     * 
     *  Several right hand side vectors b and solution vectors x can be
     *  handled in a single call; they are stored as the columns of the
     *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     *  matrix X.
     * 
     *  Arguments
     *  =========
     * 
     *  TRANS   (input) CHARACTER*1
     *          = 'N': the linear system involves A;
     *          = 'T': the linear system involves A**T.
     * 
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M >= 0.
     * 
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     * 
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of
     *          columns of the matrices B and X. NRHS >=0.
     * 
     *  A       (input/output) REAL array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit,
     *            if M >= N, A is overwritten by details of its QR
     *                       factorization as returned by SGEQRF;
     *            if M <  N, A is overwritten by details of its LQ
     *                       factorization as returned by SGELQF.
     * 
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     * 
     *  B       (input/output) REAL array, dimension (LDB,NRHS)
     *          On entry, the matrix B of right hand side vectors, stored
     *          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
     *          if TRANS = 'T'.
     *          On exit, if INFO = 0, B is overwritten by the solution
     *          vectors, stored columnwise:
     *          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
     *          squares solution vectors; the residual sum of squares for the
     *          solution in each column is given by the sum of squares of
     *          elements N+1 to M in that column;
     *          if TRANS = 'N' and m < n, rows 1 to N of B contain the
     *          minimum norm solution vectors;
     *          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
     *          minimum norm solution vectors;
     *          if TRANS = 'T' and m < n, rows 1 to M of B contain the
     *          least squares solution vectors; the residual sum of squares
     *          for the solution in each column is given by the sum of
     *          squares of elements M+1 to N in that column.
     * 
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B. LDB >= MAX(1,M,N).
     * 
     *  =====================================================================
     * 
     * </code>
     * </pre>
     *
     * @param trans
     * @param m
     * @param n
     * @param rhsCount
     * @param a
     * @param lda
     * @param b
     * @param ldb
     */
    public static void sgels(Lapack la, TTrans trans, int m, int n, int rhsCount, float[] a, int lda, float[] b, int ldb) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkValueAtLeast(lda, m, "lda");
        checkValueAtLeast(ldb, Math.max(n, m), "ldb");
        checkMinLen(a, lda * n, "a");
        checkMinLen(b, ldb * rhsCount, "b");

        intW info = new intW(0);
        float[] work = new float[1];
        la.sgels(trans.val(), m, n, rhsCount, new float[0], lda, new float[0], ldb, work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new float[(int) work[0]];
        la.sgels(trans.val(), m, n, rhsCount, a, lda, b, ldb, work, work.length, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "A does not have full rank. Least squares solution could not be computed.");
            }
        }
    }

    public static void cgels(Lapack la, TTrans trans, int m, int n, int rhsCount, float[] a, int lda, float[] b,
            int ldb) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkValueAtLeast(lda, m, "lda");
        checkValueAtLeast(ldb, Math.max(n, m), "ldb");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(b, 2 * ldb * rhsCount, "b");

        int info = la.cgels(trans.val(), m, n, rhsCount, a, lda, b, ldb);
        if (info != 0) {
            if (info < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "A does not have full rank. Least squares solution could not be computed.");
            }
        }
    }

    public static void zgels(Lapack la, TTrans trans, int m, int n, int rhsCount, double[] a, int lda, double[] b,
            int ldb) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkValueAtLeast(lda, m, "lda");
        checkValueAtLeast(ldb, Math.max(n, m), "ldb");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(b, 2 * ldb * rhsCount, "b");

        int info = la.zgels(trans.val(), m, n, rhsCount, a, lda, b, ldb);
        if (info != 0) {
            if (info < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "A does not have full rank. Least squares solution could not be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGEQLF computes a QL factorization of a real M-by-N matrix A:
     *  A = Q * L.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit,
     *          if m >= n, the lower triangle of the subarray
     *          A(m-n+1:m,1:n) contains the N-by-N lower triangular matrix L;
     *          if m <= n, the elements on and below the (n-m)-th
     *          superdiagonal contain the M-by-N lower trapezoidal matrix L;
     *          the remaining elements, with the array TAU, represent the
     *          orthogonal matrix Q as a product of elementary reflectors
     *          (see Further Details).
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     *
     *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
     *          The scalar factors of the elementary reflectors (see Further
     *          Details).
     *
     *  Further Details
     *  ===============
     *
     *  The matrix Q is represented as a product of elementary reflectors
     *
     *     Q = H(k) . . . H(2) H(1), where k = min(m,n).
     *
     *  Each H(i) has the form
     *
     *     H(i) = I - tau * v * v'
     *
     *  where tau is a real scalar, and v is a real vector with
     *  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in
     *  A(1:m-k+i-1,n-k+i), and tau in TAU(i).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param a
     * @param lda
     * @param tau
     */
    public static void dgeqlf(Lapack la, int m, int n, double[] a, int lda, double[] tau) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(tau, Math.min(m, n), "tau");

        intW info = new intW(0);
        double[] work = new double[1];
        la.dgeqlf(m, n, new double[0], lda, new double[0], work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dgeqlf(m, n, a, lda, tau, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGEQP3 computes a QR factorization with column pivoting of a
     *  matrix A:  A*P = Q*R  using Level 3 BLAS.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A. M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit, the upper triangle of the array contains the
     *          min(M,N)-by-N upper trapezoidal matrix R; the elements below
     *          the diagonal, together with the array TAU, represent the
     *          orthogonal matrix Q as a product of min(M,N) elementary
     *          reflectors.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A. LDA >= max(1,M).
     *
     *  JPVT    (input/output) INTEGER array, dimension (N)
     *          On entry, if JPVT(J).ne.0, the J-th column of A is permuted
     *          to the front of A*P (a leading column); if JPVT(J)=0,
     *          the J-th column of A is a free column.
     *          On exit, if JPVT(J)=K, then the J-th column of A*P was the
     *          the K-th column of A.
     *
     *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
     *          The scalar factors of the elementary reflectors.
     *
     *  Further Details
     *  ===============
     *
     *  The matrix Q is represented as a product of elementary reflectors
     *
     *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
     *
     *  Each H(i) has the form
     *
     *     H(i) = I - tau * v * v'
     *
     *  where tau is a real/complex scalar, and v is a real/complex vector
     *  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
     *  A(i+1:m,i), and tau in TAU(i).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param a
     * @param lda
     * @param jPivot
     * @param tau
     */
    public static void dgeqp3(Lapack la, int m, int n, double[] a, int lda, int[] jPivot, double[] tau) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(jPivot, n, "jPivot");
        checkMinLen(tau, Math.min(m, n), "tau");

        intW info = new intW(0);
        double[] work = new double[1];
        la.dgeqp3(m, n, new double[0], lda, new int[0], new double[0], work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dgeqp3(m, n, a, lda, jPivot, tau, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGEQRF computes a QR factorization of a real M-by-N matrix A:
     *  A = Q * R.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit, the elements on and above the diagonal of the array
     *          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
     *          upper triangular if m >= n); the elements below the diagonal,
     *          with the array TAU, represent the orthogonal matrix Q as a
     *          product of min(m,n) elementary reflectors (see Further
     *          Details).
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     *
     *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
     *          The scalar factors of the elementary reflectors (see Further
     *          Details).
     *
     *  Further Details
     *  ===============
     *
     *  The matrix Q is represented as a product of elementary reflectors
     *
     *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
     *
     *  Each H(i) has the form
     *
     *     H(i) = I - tau * v * v'
     *
     *  where tau is a real scalar, and v is a real vector with
     *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
     *  and tau in TAU(i).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param a
     * @param lda
     * @param tau
     */
    public static void dgeqrf(Lapack la, int m, int n, double[] a, int lda, double[] tau) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(tau, Math.min(m, n), "tau");

        intW info = new intW(0);
        double[] work = new double[1];
        la.dgeqrf(m, n, new double[0], lda, new double[0], work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dgeqrf(m, n, a, lda, tau, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  SGEQRF computes a QR factorization of a real M-by-N matrix A:
     *  A = Q * R.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     *
     *  A       (input/output) REAL array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit, the elements on and above the diagonal of the array
     *          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
     *          upper triangular if m >= n); the elements below the diagonal,
     *          with the array TAU, represent the orthogonal matrix Q as a
     *          product of min(m,n) elementary reflectors (see Further
     *          Details).
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     *
     *  TAU     (output) REAL array, dimension (min(M,N))
     *          The scalar factors of the elementary reflectors (see Further
     *          Details).
     *
     *  Further Details
     *  ===============
     *
     *  The matrix Q is represented as a product of elementary reflectors
     *
     *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
     *
     *  Each H(i) has the form
     *
     *     H(i) = I - tau * v * v'
     *
     *  where tau is a real scalar, and v is a real vector with
     *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
     *  and tau in TAU(i).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param a
     * @param lda
     * @param tau
     */
    public static void sgeqrf(Lapack la, int m, int n, float[] a, int lda, float[] tau) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(tau, Math.min(m, n), "tau");

        intW info = new intW(0);
        float[] work = new float[1];
        la.sgeqrf(m, n, new float[0], lda, new float[0], work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new float[(int) work[0]];
        la.sgeqrf(m, n, a, lda, tau, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    public static void cgeqrf(Lapack la, int m, int n, float[] a, int lda, float[] tau) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(tau, 2 * Math.min(m, n), "tau");

        int info = la.cgeqrf(m, n, a, lda, tau);
        if (info != 0) {
            throwIAEPosition(info);
        }
    }

    public static void zgeqrf(Lapack la, int m, int n, double[] a, int lda, double[] tau) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(tau, 2 * Math.min(m, n), "tau");

        int info = la.zgeqrf(m, n, a, lda, tau);
        if (info != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGERQF computes an RQ factorization of a real M-by-N matrix A:
     *  A = R * Q.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit,
     *          if m <= n, the upper triangle of the subarray
     *          A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R;
     *          if m >= n, the elements on and above the (m-n)-th subdiagonal
     *          contain the M-by-N upper trapezoidal matrix R;
     *          the remaining elements, with the array TAU, represent the
     *          orthogonal matrix Q as a product of min(m,n) elementary
     *          reflectors (see Further Details).
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     *
     *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
     *          The scalar factors of the elementary reflectors (see Further
     *          Details).
     *
     *  Further Details
     *  ===============
     *
     *  The matrix Q is represented as a product of elementary reflectors
     *
     *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
     *
     *  Each H(i) has the form
     *
     *     H(i) = I - tau * v * v'
     *
     *  where tau is a real scalar, and v is a real vector with
     *  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in
     *  A(m-k+i,1:n-k+i-1), and tau in TAU(i).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param a
     * @param lda
     * @param tau
     */
    public static void dgerqf(Lapack la, int m, int n, double[] a, int lda, double[] tau) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(tau, Math.min(m, n), "tau");

        intW info = new intW(0);
        double[] work = new double[1];
        la.dgerqf(m, n, new double[0], lda, new double[0], work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dgerqf(m, n, a, lda, tau, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGESDD computes the singular value decomposition (SVD) of a real
     *  M-by-N matrix A, optionally computing the left and right singular
     *  vectors.  If singular vectors are desired, it uses a
     *  divide-and-conquer algorithm.
     *
     *  The SVD is written
     *
     *       A = U * SIGMA * transpose(V)
     *
     *  where SIGMA is an M-by-N matrix which is zero except for its
     *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
     *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
     *  are the singular values of A; they are real and non-negative, and
     *  are returned in descending order.  The first min(m,n) columns of
     *  U and V are the left and right singular vectors of A.
     *
     *  Note that the routine returns VT = V**T, not V.
     *
     *  The divide and conquer algorithm makes very mild assumptions about
     *  floating point arithmetic. It will work on machines with a guard
     *  digit in add/subtract, or on those binary machines without guard
     *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
     *  without guard digits, but we know of none.
     *
     *  Arguments
     *  =========
     *
     *  JOBZ    (input) CHARACTER*1
     *          Specifies options for computing all or part of the matrix U:
     *          = 'A':  all M columns of U and all N rows of V**T are
     *                  returned in the arrays U and VT;
     *          = 'S':  the first min(M,N) columns of U and the first
     *                  min(M,N) rows of V**T are returned in the arrays U
     *                  and VT;
     *          = 'O':  If M >= N, the first N columns of U are overwritten
     *                  on the array A and all rows of V**T are returned in
     *                  the array VT;
     *                  otherwise, all columns of U are returned in the
     *                  array U and the first M rows of V**T are overwritten
     *                  in the array A;
     *          = 'N':  no columns of U or rows of V**T are computed.
     *
     *  M       (input) INTEGER
     *          The number of rows of the input matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the input matrix A.  N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit,
     *          if JOBZ = 'O',  A is overwritten with the first N columns
     *                          of U (the left singular vectors, stored
     *                          columnwise) if M >= N;
     *                          A is overwritten with the first M rows
     *                          of V**T (the right singular vectors, stored
     *                          rowwise) otherwise.
     *          if JOBZ .ne. 'O', the contents of A are destroyed.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     *
     *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
     *          The singular values of A, sorted so that S(i) >= S(i+1).
     *
     *  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
     *          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
     *          UCOL = min(M,N) if JOBZ = 'S'.
     *          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
     *          orthogonal matrix U;
     *          if JOBZ = 'S', U contains the first min(M,N) columns of U
     *          (the left singular vectors, stored columnwise);
     *          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
     *
     *  LDU     (input) INTEGER
     *          The leading dimension of the array U.  LDU >= 1; if
     *          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
     *
     *  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
     *          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
     *          N-by-N orthogonal matrix V**T;
     *          if JOBZ = 'S', VT contains the first min(M,N) rows of
     *          V**T (the right singular vectors, stored rowwise);
     *          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
     *
     *  LDVT    (input) INTEGER
     *          The leading dimension of the array VT.  LDVT >= 1; if
     *          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
     *          if JOBZ = 'S', LDVT >= min(M,N).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param jobz
     * @param m
     * @param n
     * @param a
     * @param lda
     * @param s
     * @param u
     * @param ldu
     * @param vt
     * @param ldvt
     */
    public static void dgesdd(Lapack la, TSvdJob jobz, int m, int n, double[] a, int lda, double[] s, double[] u,
            int ldu, double[] vt, int ldvt) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(s, Math.min(m, n), "s");
        checkStrictlyPositive(ldu, "ldu");
        checkStrictlyPositive(ldvt, "ldvt");
        // ldu + u
        if (jobz == TSvdJob.ALL || jobz == TSvdJob.PART || (m < n && jobz == TSvdJob.OVERWRITE)) {
            checkValueAtLeast(ldu, m, "ldu");
            int ucol = (jobz == TSvdJob.PART) ? Math.min(m, n) : m;
            checkMinLen(u, ldu * ucol, "u");
        }
        // ldvt + vt
        if (jobz == TSvdJob.ALL || (m >= n && jobz == TSvdJob.OVERWRITE)) {
            checkValueAtLeast(ldvt, n, "ldvt");
            checkMinLen(vt, ldvt * n, "vt");
        } else if (jobz == TSvdJob.PART) {
            checkValueAtLeast(ldvt, Math.min(m, n), "ldvt");
            checkMinLen(vt, ldvt * n, "vt");
        }

        intW info = new intW(0);
        int[] iwork = new int[8 * Math.min(m, n)];
        double[] work = new double[1];
        la.dgesdd(jobz.val(), m, n, new double[0], lda, new double[0], new double[0], ldu, new double[0], ldvt, work,
                -1, new int[0], info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dgesdd(jobz.val(), m, n, a, lda, s, u, ldu, vt, ldvt, work, work.length, iwork, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new NotConvergedException("Did not converge. Update failed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  SGESDD computes the singular value decomposition (SVD) of a real
     *  M-by-N matrix A, optionally computing the left and right singular
     *  vectors.  If singular vectors are desired, it uses a
     *  divide-and-conquer algorithm.
     *
     *  The SVD is written
     *
     *       A = U * SIGMA * transpose(V)
     *
     *  where SIGMA is an M-by-N matrix which is zero except for its
     *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
     *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
     *  are the singular values of A; they are real and non-negative, and
     *  are returned in descending order.  The first min(m,n) columns of
     *  U and V are the left and right singular vectors of A.
     *
     *  Note that the routine returns VT = V**T, not V.
     *
     *  The divide and conquer algorithm makes very mild assumptions about
     *  floating point arithmetic. It will work on machines with a guard
     *  digit in add/subtract, or on those binary machines without guard
     *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
     *  without guard digits, but we know of none.
     *
     *  Arguments
     *  =========
     *
     *  JOBZ    (input) CHARACTER*1
     *          Specifies options for computing all or part of the matrix U:
     *          = 'A':  all M columns of U and all N rows of V**T are
     *                  returned in the arrays U and VT;
     *          = 'S':  the first min(M,N) columns of U and the first
     *                  min(M,N) rows of V**T are returned in the arrays U
     *                  and VT;
     *          = 'O':  If M >= N, the first N columns of U are overwritten
     *                  on the array A and all rows of V**T are returned in
     *                  the array VT;
     *                  otherwise, all columns of U are returned in the
     *                  array U and the first M rows of V**T are overwritten
     *                  in the array A;
     *          = 'N':  no columns of U or rows of V**T are computed.
     *
     *  M       (input) INTEGER
     *          The number of rows of the input matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the input matrix A.  N >= 0.
     *
     *  A       (input/output) REAL array, dimension (LDA,N)
     *          On entry, the M-by-N matrix A.
     *          On exit,
     *          if JOBZ = 'O',  A is overwritten with the first N columns
     *                          of U (the left singular vectors, stored
     *                          columnwise) if M >= N;
     *                          A is overwritten with the first M rows
     *                          of V**T (the right singular vectors, stored
     *                          rowwise) otherwise.
     *          if JOBZ .ne. 'O', the contents of A are destroyed.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     *
     *  S       (output) REAL array, dimension (min(M,N))
     *          The singular values of A, sorted so that S(i) >= S(i+1).
     *
     *  U       (output) REAL array, dimension (LDU,UCOL)
     *          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
     *          UCOL = min(M,N) if JOBZ = 'S'.
     *          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
     *          orthogonal matrix U;
     *          if JOBZ = 'S', U contains the first min(M,N) columns of U
     *          (the left singular vectors, stored columnwise);
     *          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
     *
     *  LDU     (input) INTEGER
     *          The leading dimension of the array U.  LDU >= 1; if
     *          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
     *
     *  VT      (output) REAL array, dimension (LDVT,N)
     *          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
     *          N-by-N orthogonal matrix V**T;
     *          if JOBZ = 'S', VT contains the first min(M,N) rows of
     *          V**T (the right singular vectors, stored rowwise);
     *          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
     *
     *  LDVT    (input) INTEGER
     *          The leading dimension of the array VT.  LDVT >= 1; if
     *          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
     *          if JOBZ = 'S', LDVT >= min(M,N).
     *
     *  Further Details
     *  ===============
     *
     *  Based on contributions by
     *     Ming Gu and Huan Ren, Computer Science Division, University of
     *     California at Berkeley, USA
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param jobz
     * @param m
     * @param n
     * @param a
     * @param lda
     * @param s
     * @param u
     * @param ldu
     * @param vt
     * @param ldvt
     */
    public static void sgesdd(Lapack la, TSvdJob jobz, int m, int n, float[] a, int lda, float[] s, float[] u, int ldu,
            float[] vt, int ldvt) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(s, Math.min(m, n), "s");
        checkStrictlyPositive(ldu, "ldu");
        checkStrictlyPositive(ldvt, "ldvt");
        // ldu + u
        if (jobz == TSvdJob.ALL || jobz == TSvdJob.PART || (m < n && jobz == TSvdJob.OVERWRITE)) {
            checkValueAtLeast(ldu, m, "ldu");
            int ucol = (jobz == TSvdJob.PART) ? Math.min(m, n) : m;
            checkMinLen(u, ldu * ucol, "u");
        }
        // ldvt + vt
        if (jobz == TSvdJob.ALL || (m >= n && jobz == TSvdJob.OVERWRITE)) {
            checkValueAtLeast(ldvt, n, "ldvt");
            checkMinLen(vt, ldvt * n, "vt");
        } else if (jobz == TSvdJob.PART) {
            checkValueAtLeast(ldvt, Math.min(m, n), "ldvt");
            checkMinLen(vt, ldvt * n, "vt");
        }

        intW info = new intW(0);
        int[] iwork = new int[8 * Math.min(m, n)];
        float[] work = new float[1];
        la.sgesdd(jobz.val(), m, n, new float[0], lda, new float[0], new float[0], ldu, new float[0], ldvt, work, -1,
                new int[0], info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new float[(int) work[0]];
        la.sgesdd(jobz.val(), m, n, a, lda, s, u, ldu, vt, ldvt, work, work.length, iwork, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new NotConvergedException("Did not converge. Update failed.");
            }
        }
    }

    // note: 's' (the singular values) are real (a float[]) wheras 'a[]', 'u[]'
    // and 'vt[]' are arrays of complex numbers; the only case where 'a' must
    // be copied before calling cgesdd is when jobz == 'O'
    public static void cgesdd(Lapack la, TSvdJob jobz, int m, int n, float[] a, int lda, float[] s, float[] u, int ldu,
            float[] vt, int ldvt) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(s, Math.min(m, n), "s");
        checkStrictlyPositive(ldu, "ldu");
        checkStrictlyPositive(ldvt, "ldvt");
        // ldu + u
        if (jobz == TSvdJob.ALL || jobz == TSvdJob.PART || (m < n && jobz == TSvdJob.OVERWRITE)) {
            checkValueAtLeast(ldu, m, "ldu");
            int ucol = (jobz == TSvdJob.PART) ? Math.min(m, n) : m;
            checkMinLen(u, 2 * ldu * ucol, "u");
        }
        // ldvt + vt
        if (jobz == TSvdJob.ALL || (m >= n && jobz == TSvdJob.OVERWRITE)) {
            checkValueAtLeast(ldvt, n, "ldvt");
            checkMinLen(vt, 2 * ldvt * n, "vt");
        } else if (jobz == TSvdJob.PART) {
            checkValueAtLeast(ldvt, Math.min(m, n), "ldvt");
            checkMinLen(vt, 2 * ldvt * n, "vt");
        }

        int info = la.cgesdd(jobz.val(), m, n, a, lda, s, u, ldu, vt, ldvt);
        if (info != 0) {
            if (info < 0) {
                throwIAEPosition(info);
            } else {
                throw new NotConvergedException("Did not converge. Update failed.");
            }
        }
    }

    // note: 's' (the singular values) are real (a double[]) wheras 'a[]', 'u[]'
    // and 'vt[]' are arrays of complex numbers; the only case where 'a' must
    // be copied before calling zgesdd is when jobz == 'O'
    public static void zgesdd(Lapack la, TSvdJob jobz, int m, int n, double[] a, int lda, double[] s, double[] u,
            int ldu, double[] vt, int ldvt) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(s, Math.min(m, n), "s");
        checkStrictlyPositive(ldu, "ldu");
        checkStrictlyPositive(ldvt, "ldvt");
        // ldu + u
        if (jobz == TSvdJob.ALL || jobz == TSvdJob.PART || (m < n && jobz == TSvdJob.OVERWRITE)) {
            checkValueAtLeast(ldu, m, "ldu");
            int ucol = (jobz == TSvdJob.PART) ? Math.min(m, n) : m;
            checkMinLen(u, 2 * ldu * ucol, "u");
        }
        // ldvt + vt
        if (jobz == TSvdJob.ALL || (m >= n && jobz == TSvdJob.OVERWRITE)) {
            checkValueAtLeast(ldvt, n, "ldvt");
            checkMinLen(vt, 2 * ldvt * n, "vt");
        } else if (jobz == TSvdJob.PART) {
            checkValueAtLeast(ldvt, Math.min(m, n), "ldvt");
            checkMinLen(vt, 2 * ldvt * n, "vt");
        }

        int info = la.zgesdd(jobz.val(), m, n, a, lda, s, u, ldu, vt, ldvt);
        if (info != 0) {
            if (info < 0) {
                throwIAEPosition(info);
            } else {
                throw new NotConvergedException("Did not converge. Update failed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGESV computes the solution to a real system of linear equations
     *     A * X = B,
     *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     *
     *  The LU decomposition with partial pivoting and row interchanges is
     *  used to factor A as
     *     A = P * L * U,
     *  where P is a permutation matrix, L is unit lower triangular, and U is
     *  upper triangular.  The factored form of A is then used to solve the
     *  system of equations A * X = B.
     *
     *  Arguments
     *  =========
     *
     *  N       (input) INTEGER
     *          The number of linear equations, i.e., the order of the
     *          matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the N-by-N coefficient matrix A.
     *          On exit, the factors L and U from the factorization
     *          A = P*L*U; the unit diagonal elements of L are not stored.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  IPIV    (output) INTEGER array, dimension (N)
     *          The pivot indices that define the permutation matrix P;
     *          row i of the matrix was interchanged with row IPIV(i).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the N-by-NRHS matrix of right hand side matrix B.
     *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param rhsCount
     * @param a
     * @param lda
     * @param indices
     * @param b
     * @param ldb
     */
    public static void dgesv(Lapack la, int n, int rhsCount, double[] a, int lda, int[] indices, double[] b, int ldb) {
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkValueAtLeast(lda, n, "lda");
        checkValueAtLeast(ldb, n, "ldb");
        checkMinLen(indices, n, "indices");
        checkMinLen(a, lda * n, "a");
        checkMinLen(b, ldb * rhsCount, "b");

        intW info = new intW(0);
        la.dgesv(n, rhsCount, a, lda, indices, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "Factor U in the LU decomposition is exactly singular. Solution could not be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     * 
     *  Purpose
     *  =======
     * 
     *  SGESV computes the solution to a real system of linear equations
     *     A * X = B,
     *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     * 
     *  The LU decomposition with partial pivoting and row interchanges is
     *  used to factor A as
     *     A = P * L * U,
     *  where P is a permutation matrix, L is unit lower triangular, and U is
     *  upper triangular.  The factored form of A is then used to solve the
     *  system of equations A * X = B.
     * 
     *  Arguments
     *  =========
     * 
     *  N       (input) INTEGER
     *          The number of linear equations, i.e., the order of the
     *          matrix A.  N >= 0.
     * 
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     * 
     *  A       (input/output) REAL array, dimension (LDA,N)
     *          On entry, the N-by-N coefficient matrix A.
     *          On exit, the factors L and U from the factorization
     *          A = P*L*U; the unit diagonal elements of L are not stored.
     * 
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     * 
     *  IPIV    (output) INTEGER array, dimension (N)
     *          The pivot indices that define the permutation matrix P;
     *          row i of the matrix was interchanged with row IPIV(i).
     * 
     *  B       (input/output) REAL array, dimension (LDB,NRHS)
     *          On entry, the N-by-NRHS matrix of right hand side matrix B.
     *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     * 
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     * 
     *  =====================================================================
     * 
     * </code>
     * </pre>
     *
     * @param n
     * @param rhsCount
     * @param a
     * @param lda
     * @param indices
     * @param b
     * @param ldb
     */
    public static void sgesv(Lapack la, int n, int rhsCount, float[] a, int lda, int[] indices, float[] b, int ldb) {
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkValueAtLeast(lda, n, "lda");
        checkValueAtLeast(ldb, n, "ldb");
        checkMinLen(indices, n, "indices");
        checkMinLen(a, lda * n, "a");
        checkMinLen(b, ldb * rhsCount, "b");

        intW info = new intW(0);
        la.sgesv(n, rhsCount, a, lda, indices, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "Factor U in the LU decomposition is exactly singular. Solution could not be computed.");
            }
        }
    }

    public static void cgesv(Lapack la, int n, int rhsCount, float[] a, int lda, int[] indices, float[] b, int ldb) {
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkValueAtLeast(lda, n, "lda");
        checkValueAtLeast(ldb, n, "ldb");
        checkMinLen(indices, n, "indices");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(b, 2 * ldb * rhsCount, "b");

        int info = la.cgesv(n, rhsCount, a, lda, indices, b, ldb);
        if (info != 0) {
            if (info < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "Factor U in the LU decomposition is exactly singular. Solution could not be computed.");
            }
        }
    }

    public static void zgesv(Lapack la, int n, int rhsCount, double[] a, int lda, int[] indices, double[] b, int ldb) {
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkValueAtLeast(lda, n, "lda");
        checkValueAtLeast(ldb, n, "ldb");
        checkMinLen(indices, n, "indices");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(b, 2 * ldb * rhsCount, "b");

        int info = la.zgesv(n, rhsCount, a, lda, indices, b, ldb);
        if (info != 0) {
            if (info < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "Factor U in the LU decomposition is exactly singular. Solution could not be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGETRF computes an LU factorization of a general M-by-N matrix A
     *  using partial pivoting with row interchanges.
     *
     *  The factorization has the form
     *     A = P * L * U
     *  where P is a permutation matrix, L is lower triangular with unit
     *  diagonal elements (lower trapezoidal if m > n), and U is upper
     *  triangular (upper trapezoidal if m < n).
     *
     *  This is the right-looking Level 3 BLAS version of the algorithm.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the M-by-N matrix to be factored.
     *          On exit, the factors L and U from the factorization
     *          A = P*L*U; the unit diagonal elements of L are not stored.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     *
     *  IPIV    (output) INTEGER array, dimension (min(M,N))
     *          The pivot indices; for 1 <= i <= min(M,N), row i of the
     *          matrix was interchanged with row IPIV(i).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param a
     * @param lda
     * @param indices
     */
    public static void dgetrf(Lapack la, int m, int n, double[] a, int lda, int[] indices) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(indices, Math.min(m, n), "indices");

        intW info = new intW(0);
        la.dgetrf(m, n, a, lda, indices, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("Factor U in the LU decomposition is exactly singular");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  SGETRF computes an LU factorization of a general M-by-N matrix A
     *  using partial pivoting with row interchanges.
     *
     *  The factorization has the form
     *     A = P * L * U
     *  where P is a permutation matrix, L is lower triangular with unit
     *  diagonal elements (lower trapezoidal if m > n), and U is upper
     *  triangular (upper trapezoidal if m < n).
     *
     *  This is the right-looking Level 3 BLAS version of the algorithm.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix A.  M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.  N >= 0.
     *
     *  A       (input/output) REAL array, dimension (LDA,N)
     *          On entry, the M-by-N matrix to be factored.
     *          On exit, the factors L and U from the factorization
     *          A = P*L*U; the unit diagonal elements of L are not stored.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,M).
     *
     *  IPIV    (output) INTEGER array, dimension (min(M,N))
     *          The pivot indices; for 1 <= i <= min(M,N), row i of the
     *          matrix was interchanged with row IPIV(i).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param a
     * @param lda
     * @param indices
     */
    public static void sgetrf(Lapack la, int m, int n, float[] a, int lda, int[] indices) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, lda * n, "a");
        checkMinLen(indices, Math.min(m, n), "indices");

        intW info = new intW(0);
        la.sgetrf(m, n, a, lda, indices, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("Factor U in the LU decomposition is exactly singular");
            }
        }
    }

    public static void cgetrf(Lapack la, int m, int n, float[] a, int lda, int[] indices) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(indices, Math.min(m, n), "indices");

        int info = la.cgetrf(m, n, a, lda, indices);
        if (info != 0) {
            if (info < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("Factor U in the LU decomposition is exactly singular");
            }
        }
    }

    public static void zgetrf(Lapack la, int m, int n, double[] a, int lda, int[] indices) {
        checkStrictlyPositive(m, "m");
        checkStrictlyPositive(n, "n");
        checkValueAtLeast(lda, m, "lda");
        checkMinLen(a, 2 * lda * n, "a");
        checkMinLen(indices, Math.min(m, n), "indices");

        int info = la.zgetrf(m, n, a, lda, indices);
        if (info != 0) {
            if (info < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("Factor U in the LU decomposition is exactly singular");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGETRS solves a system of linear equations
     *     A * X = B  or  A' * X = B
     *  with a general N-by-N matrix A using the LU factorization computed
     *  by DGETRF.
     *
     *  Arguments
     *  =========
     *
     *  TRANS   (input) CHARACTER*1
     *          Specifies the form of the system of equations:
     *          = 'N':  A * X = B  (No transpose)
     *          = 'T':  A'* X = B  (Transpose)
     *          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
     *          The factors L and U from the factorization A = P*L*U
     *          as computed by DGETRF.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  IPIV    (input) INTEGER array, dimension (N)
     *          The pivot indices from DGETRF; for 1<=i<=N, row i of the
     *          matrix was interchanged with row IPIV(i).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the right hand side matrix B.
     *          On exit, the solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param trans
     * @param n
     * @param rhsCount
     * @param a
     * @param lda
     * @param indices
     * @param b
     * @param ldb
     */
    public static void dgetrs(Lapack la, TTrans trans, int n, int rhsCount, double[] a, int lda, int[] indices,
            double[] b, int ldb) {
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkValueAtLeast(lda, n, "lda");
        checkValueAtLeast(ldb, n, "ldb");
        checkMinLen(a, lda * n, "a");
        checkMinLen(b, ldb * rhsCount, "b");
        checkMinLen(indices, n, "indices");

        intW info = new intW(0);
        la.dgetrs(trans.val(), n, rhsCount, a, lda, indices, b, ldb, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DGTSV  solves the equation
     *
     *     A*X = B,
     *
     *  where A is an n by n tridiagonal matrix, by Gaussian elimination with
     *  partial pivoting.
     *
     *  Note that the equation  A'*X = B  may be solved by interchanging the
     *  order of the arguments DU and DL.
     *
     *  Arguments
     *  =========
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
     *          On entry, DL must contain the (n-1) sub-diagonal elements of
     *          A.
     *
     *          On exit, DL is overwritten by the (n-2) elements of the
     *          second super-diagonal of the upper triangular matrix U from
     *          the LU factorization of A, in DL(1), ..., DL(n-2).
     *
     *  D       (input/output) DOUBLE PRECISION array, dimension (N)
     *          On entry, D must contain the diagonal elements of A.
     *
     *          On exit, D is overwritten by the n diagonal elements of U.
     *
     *  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
     *          On entry, DU must contain the (n-1) super-diagonal elements
     *          of A.
     *
     *          On exit, DU is overwritten by the (n-1) elements of the first
     *          super-diagonal of U.
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the N by NRHS matrix of right hand side matrix B.
     *          On exit, if INFO = 0, the N by NRHS solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param rhsCount
     * @param dl
     * @param d
     * @param du
     * @param b
     * @param ldb
     */
    public static void dgtsv(Lapack la, int n, int rhsCount, double[] dl, double[] d, double[] du, double[] b,
            int ldb) {
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(rhsCount, "rhsCount");
        checkValueAtLeast(ldb, n, "ldb");
        checkMinLen(b, ldb * rhsCount, "b");
        checkMinLen(dl, n - 1, "dl");
        checkMinLen(du, n - 1, "du");
        checkMinLen(d, n, "d");

        intW info = new intW(0);
        la.dgtsv(n, rhsCount, dl, d, du, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "Factor U in the LU decomposition is exactly singular. Solution could not be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DLASWP performs a series of row interchanges on the matrix A.
     *  One row interchange is initiated for each of rows K1 through K2 of A.
     *
     *  Arguments
     *  =========
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix A.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the matrix of column dimension N to which the row
     *          interchanges will be applied.
     *          On exit, the permuted matrix.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.
     *
     *  K1      (input) INTEGER
     *          The first element of IPIV for which a row interchange will
     *          be done.
     *
     *  K2      (input) INTEGER
     *          The last element of IPIV for which a row interchange will
     *          be done.
     *
     *  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
     *          The vector of pivot indices.  Only the elements in positions
     *          K1 through K2 of IPIV are accessed.
     *          IPIV(K) = L implies rows K and L are to be interchanged.
     *
     *  INCX    (input) INTEGER
     *          The increment between successive values of IPIV.  If IPIV
     *          is negative, the pivots are applied in reverse order.
     *
     * =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param a
     * @param lda
     * @param pivFirstIdx
     * @param pivLastIdx
     * @param indices
     * @param increment
     */
    public static void dlaswp(Lapack la, int n, double[] a, int lda, int pivFirstIdx, int pivLastIdx, int[] indices,
            int increment) {
        checkStrictlyPositive(n, "n");
        checkStrictlyPositive(lda, "lda");
        checkMinLen(a, lda * n, "a");
        checkStrictlyPositive(pivFirstIdx, "pivFirstIdx");
        checkStrictlyPositive(pivLastIdx, "pivLastIdx");
        checkMinLen(indices, pivLastIdx * Math.abs(increment), "indices");

        la.dlaswp(n, a, lda, pivFirstIdx, pivLastIdx, indices, increment);
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DORGLQ generates an M-by-N real matrix Q with orthonormal rows,
     *  which is defined as the first M rows of a product of K elementary
     *  reflectors of order N
     *
     *        Q  =  H(k) . . . H(2) H(1)
     *
     *  as returned by DGELQF.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix Q. M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix Q. N >= M.
     *
     *  K       (input) INTEGER
     *          The number of elementary reflectors whose product defines the
     *          matrix Q. M >= K >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the i-th row must contain the vector which defines
     *          the elementary reflector H(i), for i = 1,2,...,k, as returned
     *          by DGELQF in the first k rows of its array argument A.
     *          On exit, the M-by-N matrix Q.
     *
     *  LDA     (input) INTEGER
     *          The first dimension of the array A. LDA >= max(1,M).
     *
     *  TAU     (input) DOUBLE PRECISION array, dimension (K)
     *          TAU(i) must contain the scalar factor of the elementary
     *          reflector H(i), as returned by DGELQF.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param k
     * @param a
     * @param lda
     * @param tau
     */
    public static void dorglq(Lapack la, int m, int n, int k, double[] a, int lda, double[] tau) {
        checkStrictlyPositive(m, "m");
        checkValueAtLeast(m, k, "m");
        checkValueAtLeast(n, m, "n");
        checkStrictlyPositive(k, "k");
        checkValueAtLeast(lda, Math.max(1, m), "lda");
        checkMinLen(a, m * n, "a");
        checkMinLen(tau, k, "tau");

        intW info = new intW(0);
        double[] work = new double[1];
        la.dorglq(m, n, k, new double[0], lda, new double[0], work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dorglq(m, n, k, a, lda, tau, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
     *  which is defined as the first N columns of a product of K elementary
     *  reflectors of order M
     *
     *        Q  =  H(1) H(2) . . . H(k)
     *
     *  as returned by DGEQRF.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix Q. M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix Q. M >= N >= 0.
     *
     *  K       (input) INTEGER
     *          The number of elementary reflectors whose product defines the
     *          matrix Q. N >= K >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the i-th column must contain the vector which
     *          defines the elementary reflector H(i), for i = 1,2,...,k, as
     *          returned by DGEQRF in the first k columns of its array
     *          argument A.
     *          On exit, the M-by-N matrix Q.
     *
     *  LDA     (input) INTEGER
     *          The first dimension of the array A. LDA >= max(1,M).
     *
     *  TAU     (input) DOUBLE PRECISION array, dimension (K)
     *          TAU(i) must contain the scalar factor of the elementary
     *          reflector H(i), as returned by DGEQRF.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param k
     * @param a
     * @param lda
     * @param tau
     */
    public static void dorgqr(Lapack la, int m, int n, int k, double[] a, int lda, double[] tau) {
        checkStrictlyPositive(k, "k");
        checkValueAtLeast(n, k, "n");
        checkValueAtLeast(m, n, "m");
        checkValueAtLeast(lda, Math.max(1, m), "lda");
        checkMinLen(a, m * n, "a");
        checkMinLen(tau, k, "tau");

        intW info = new intW(0);
        double[] work = new double[1];
        la.dorgqr(m, n, k, new double[0], lda, new double[0], work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dorgqr(m, n, k, a, lda, tau, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * Complex counterpart of {@link #dorgqr}.
     */
    public static void zungqr(Lapack la, int m, int n, int k, double[] a, int lda, double[] tau) {
        checkStrictlyPositive(k, "k");
        checkValueAtLeast(n, k, "n");
        checkValueAtLeast(m, n, "m");
        checkValueAtLeast(lda, Math.max(1, m), "lda");
        checkMinLen(a, 2 * m * n, "a");
        checkMinLen(tau, 2 * k, "tau");

        int info = la.zungqr(m, n, k, a, lda, tau);
        if (info != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  SORGQR generates an M-by-N real matrix Q with orthonormal columns,
     *  which is defined as the first N columns of a product of K elementary
     *  reflectors of order M
     *
     *        Q  =  H(1) H(2) . . . H(k)
     *
     *  as returned by SGEQRF.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix Q. M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix Q. M >= N >= 0.
     *
     *  K       (input) INTEGER
     *          The number of elementary reflectors whose product defines the
     *          matrix Q. N >= K >= 0.
     *
     *  A       (input/output) REAL array, dimension (LDA,N)
     *          On entry, the i-th column must contain the vector which
     *          defines the elementary reflector H(i), for i = 1,2,...,k, as
     *          returned by SGEQRF in the first k columns of its array
     *          argument A.
     *          On exit, the M-by-N matrix Q.
     *
     *  LDA     (input) INTEGER
     *          The first dimension of the array A. LDA >= max(1,M).
     *
     *  TAU     (input) REAL array, dimension (K)
     *          TAU(i) must contain the scalar factor of the elementary
     *          reflector H(i), as returned by SGEQRF.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param k
     * @param a
     * @param lda
     * @param tau
     */
    public static void sorgqr(Lapack la, int m, int n, int k, float[] a, int lda, float[] tau) {
        checkStrictlyPositive(k, "k");
        checkValueAtLeast(n, k, "n");
        checkValueAtLeast(m, n, "m");
        checkValueAtLeast(lda, Math.max(1, m), "lda");
        checkMinLen(a, m * n, "a");
        checkMinLen(tau, k, "tau");

        intW info = new intW(0);
        float[] work = new float[1];
        la.sorgqr(m, n, k, new float[0], lda, new float[0], work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new float[(int) work[0]];
        la.sorgqr(m, n, k, a, lda, tau, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * Complex counterpart of {@link #sorgqr}.
     */
    public static void cungqr(Lapack la, int m, int n, int k, float[] a, int lda, float[] tau) {
        checkStrictlyPositive(k, "k");
        checkValueAtLeast(n, k, "n");
        checkValueAtLeast(m, n, "m");
        checkValueAtLeast(lda, Math.max(1, m), "lda");
        checkMinLen(a, 2 * m * n, "a");
        checkMinLen(tau, 2 * k, "tau");

        int info = la.cungqr(m, n, k, a, lda, tau);
        if (info != 0) {
            throwIAEPosition(info);
        }
    }

    // TODO: finish sanity checks for the remaining arguments

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DORGRQ generates an M-by-N real matrix Q with orthonormal rows,
     *  which is defined as the last M rows of a product of K elementary
     *  reflectors of order N
     *
     *        Q  =  H(1) H(2) . . . H(k)
     *
     *  as returned by DGERQF.
     *
     *  Arguments
     *  =========
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix Q. M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix Q. N >= M.
     *
     *  K       (input) INTEGER
     *          The number of elementary reflectors whose product defines the
     *          matrix Q. M >= K >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the (m-k+i)-th row must contain the vector which
     *          defines the elementary reflector H(i), for i = 1,2,...,k, as
     *          returned by DGERQF in the last k rows of its array argument
     *          A.
     *          On exit, the M-by-N matrix Q.
     *
     *  LDA     (input) INTEGER
     *          The first dimension of the array A. LDA >= max(1,M).
     *
     *  TAU     (input) DOUBLE PRECISION array, dimension (K)
     *          TAU(i) must contain the scalar factor of the elementary
     *          reflector H(i), as returned by DGERQF.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param m
     * @param n
     * @param k
     * @param a
     * @param lda
     * @param tau
     */
    public static void dorgrq(Lapack la, int m, int n, int k, double[] a, int lda, double[] tau) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        double[] work = new double[1];
        la.dorgrq(m, n, k, new double[0], lda, new double[0], work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dorgrq(m, n, k, a, lda, tau, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DORMRZ overwrites the general real M-by-N matrix C with
     *
     *                  SIDE = 'L'     SIDE = 'R'
     *  TRANS = 'N':      Q * C          C * Q
     *  TRANS = 'T':      Q**T * C       C * Q**T
     *
     *  where Q is a real orthogonal matrix defined as the product of k
     *  elementary reflectors
     *
     *        Q = H(1) H(2) . . . H(k)
     *
     *  as returned by DTZRZF. Q is of order M if SIDE = 'L' and of order N
     *  if SIDE = 'R'.
     *
     *  Arguments
     *  =========
     *
     *  SIDE    (input) CHARACTER*1
     *          = 'L': apply Q or Q**T from the Left;
     *          = 'R': apply Q or Q**T from the Right.
     *
     *  TRANS   (input) CHARACTER*1
     *          = 'N':  No transpose, apply Q;
     *          = 'T':  Transpose, apply Q**T.
     *
     *  M       (input) INTEGER
     *          The number of rows of the matrix C. M >= 0.
     *
     *  N       (input) INTEGER
     *          The number of columns of the matrix C. N >= 0.
     *
     *  K       (input) INTEGER
     *          The number of elementary reflectors whose product defines
     *          the matrix Q.
     *          If SIDE = 'L', M >= K >= 0;
     *          if SIDE = 'R', N >= K >= 0.
     *
     *  L       (input) INTEGER
     *          The number of columns of the matrix A containing
     *          the meaningful part of the Householder reflectors.
     *          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
     *
     *  A       (input) DOUBLE PRECISION array, dimension
     *                               (LDA,M) if SIDE = 'L',
     *                               (LDA,N) if SIDE = 'R'
     *          The i-th row must contain the vector which defines the
     *          elementary reflector H(i), for i = 1,2,...,k, as returned by
     *          DTZRZF in the last k rows of its array argument A.
     *          A is modified by the routine but restored on exit.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A. LDA >= max(1,K).
     *
     *  TAU     (input) DOUBLE PRECISION array, dimension (K)
     *          TAU(i) must contain the scalar factor of the elementary
     *          reflector H(i), as returned by DTZRZF.
     *
     *  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
     *          On entry, the M-by-N matrix C.
     *          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
     *
     *  LDC     (input) INTEGER
     *          The leading dimension of the array C. LDC >= max(1,M).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param side
     * @param trans
     * @param m
     * @param n
     * @param k
     * @param l
     * @param a
     * @param lda
     * @param tau
     * @param c
     * @param ldc
     */
    public static void dormrz(Lapack la, TSide side, TTrans trans, int m, int n, int k, int l, double[] a, int lda,
            double[] tau, double[] c, int ldc) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        double[] work = new double[1];
        la.dormrz(side.val(), trans.val(), m, n, k, l, new double[0], lda, new double[0], new double[0], ldc, work, -1,
                info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dormrz(side.val(), trans.val(), m, n, k, l, a, lda, tau, c, ldc, work, work.length, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPBCON estimates the reciprocal of the condition number (in the
     *  1-norm) of a real symmetric positive definite band matrix using the
     *  Cholesky factorization A = U**T*U or A = L*L**T computed by DPBTRF.
     *
     *  An estimate is obtained for norm(inv(A)), and the reciprocal of the
     *  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangular factor stored in AB;
     *          = 'L':  Lower triangular factor stored in AB.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  KD      (input) INTEGER
     *          The number of superdiagonals of the matrix A if UPLO = 'U',
     *          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
     *
     *  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
     *          The triangular factor U or L from the Cholesky factorization
     *          A = U**T*U or A = L*L**T of the band matrix A, stored in the
     *          first KD+1 rows of the array.  The j-th column of U or L is
     *          stored in the j-th column of the array AB as follows:
     *          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
     *          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
     *
     *  LDAB    (input) INTEGER
     *          The leading dimension of the array AB.  LDAB >= KD+1.
     *
     *  ANORM   (input) DOUBLE PRECISION
     *          The 1-norm (or infinity-norm) of the symmetric band matrix A.
     *
     *  RCOND   (output) DOUBLE PRECISION
     *          The reciprocal of the condition number of the matrix A,
     *          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
     *          estimate of the 1-norm of inv(A) computed in this routine.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param diagCount
     * @param ab
     * @param normA
     */
    public static double dpbcon(Lapack la, TUpLo uplo, int n, int diagCount, double[] ab, double normA) {
        intW info = new intW(0);
        doubleW rcond = new doubleW(0.0);
        // TODO: remaining params sanity check
        double[] work = new double[3 * n];
        int[] iwork = new int[n];
        int ldab = Math.max(1, diagCount + 1);
        la.dpbcon(uplo.val(), n, diagCount, ab, ldab, normA, rcond, work, iwork, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        return rcond.val;
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPBSV computes the solution to a real system of linear equations
     *     A * X = B,
     *  where A is an N-by-N symmetric positive definite band matrix and X
     *  and B are N-by-NRHS matrices.
     *
     *  The Cholesky decomposition is used to factor A as
     *     A = U**T * U,  if UPLO = 'U', or
     *     A = L * L**T,  if UPLO = 'L',
     *  where U is an upper triangular band matrix, and L is a lower
     *  triangular band matrix, with the same number of superdiagonals or
     *  subdiagonals as A.  The factored form of A is then used to solve the
     *  system of equations A * X = B.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The number of linear equations, i.e., the order of the
     *          matrix A.  N >= 0.
     *
     *  KD      (input) INTEGER
     *          The number of superdiagonals of the matrix A if UPLO = 'U',
     *          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
     *          On entry, the upper or lower triangle of the symmetric band
     *          matrix A, stored in the first KD+1 rows of the array.  The
     *          j-th column of A is stored in the j-th column of the array AB
     *          as follows:
     *          if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j;
     *          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD).
     *          See below for further details.
     *
     *          On exit, if INFO = 0, the triangular factor U or L from the
     *          Cholesky factorization A = U**T*U or A = L*L**T of the band
     *          matrix A, in the same storage format as A.
     *
     *  LDAB    (input) INTEGER
     *          The leading dimension of the array AB.  LDAB >= KD+1.
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the N-by-NRHS right hand side matrix B.
     *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the leading minor of order i of A is not
     *                positive definite, so the factorization could not be
     *                completed, and the solution has not been computed.
     *
     *  Further Details
     *  ===============
     *
     *  The band storage scheme is illustrated by the following example, when
     *  N = 6, KD = 2, and UPLO = 'U':
     *
     *  On entry:                       On exit:
     *
     *      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
     *      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
     *     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
     *
     *  Similarly, if UPLO = 'L' the format of A is as follows:
     *
     *  On entry:                       On exit:
     *
     *     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
     *     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
     *     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
     *
     *  Array elements marked * are not used by the routine.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param diagCount
     * @param rhsCount
     * @param ab
     * @param b
     * @param ldb
     */
    public static void dpbsv(Lapack la, TUpLo uplo, int n, int diagCount, int rhsCount, double[] ab, double[] b,
            int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        int ldab = Math.max(1, diagCount + 1);
        la.dpbsv(uplo.val(), n, diagCount, rhsCount, ab, ldab, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "A leading minor of A is not positive-definite. Solution could not be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPBTRF computes the Cholesky factorization of a real symmetric
     *  positive definite band matrix A.
     *
     *  The factorization has the form
     *     A = U**T * U,  if UPLO = 'U', or
     *     A = L  * L**T,  if UPLO = 'L',
     *  where U is an upper triangular matrix and L is lower triangular.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  KD      (input) INTEGER
     *          The number of superdiagonals of the matrix A if UPLO = 'U',
     *          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
     *
     *  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
     *          On entry, the upper or lower triangle of the symmetric band
     *          matrix A, stored in the first KD+1 rows of the array.  The
     *          j-th column of A is stored in the j-th column of the array AB
     *          as follows:
     *          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
     *          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
     *
     *          On exit, if INFO = 0, the triangular factor U or L from the
     *          Cholesky factorization A = U**T*U or A = L*L**T of the band
     *          matrix A, in the same storage format as A.
     *
     *  LDAB    (input) INTEGER
     *          The leading dimension of the array AB.  LDAB >= KD+1.
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the leading minor of order i is not
     *                positive definite, and the factorization could not be
     *                completed.
     *
     *  Further Details
     *  ===============
     *
     *  The band storage scheme is illustrated by the following example, when
     *  N = 6, KD = 2, and UPLO = 'U':
     *
     *  On entry:                       On exit:
     *
     *      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
     *      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
     *     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
     *
     *  Similarly, if UPLO = 'L' the format of A is as follows:
     *
     *  On entry:                       On exit:
     *
     *     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
     *     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
     *     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
     *
     *  Array elements marked * are not used by the routine.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param diagCount
     * @param ab
     */
    public static void dpbtrf(Lapack la, TUpLo uplo, int n, int diagCount, double[] ab) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        int ldab = Math.max(1, diagCount + 1);
        la.dpbtrf(uplo.val(), n, diagCount, ab, ldab, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("The leading minor of order " + info.val
                        + " is not positive-definite. The factorization could not be completed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPBTRS solves a system of linear equations A*X = B with a symmetric
     *  positive definite band matrix A using the Cholesky factorization
     *  A = U**T*U or A = L*L**T computed by DPBTRF.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangular factor stored in AB;
     *          = 'L':  Lower triangular factor stored in AB.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  KD      (input) INTEGER
     *          The number of superdiagonals of the matrix A if UPLO = 'U',
     *          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
     *          The triangular factor U or L from the Cholesky factorization
     *          A = U**T*U or A = L*L**T of the band matrix A, stored in the
     *          first KD+1 rows of the array.  The j-th column of U or L is
     *          stored in the j-th column of the array AB as follows:
     *          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
     *          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
     *
     *  LDAB    (input) INTEGER
     *          The leading dimension of the array AB.  LDAB >= KD+1.
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the right hand side matrix B.
     *          On exit, the solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param diagCount
     * @param rhsCount
     * @param ab
     * @param b
     * @param ldb
     */
    public static void dpbtrs(Lapack la, TUpLo uplo, int n, int diagCount, int rhsCount, double[] ab, double[] b,
            int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        int ldab = Math.max(1, diagCount + 1);
        la.dpbtrs(uplo.val(), n, diagCount, rhsCount, ab, ldab, b, ldb, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPOCON estimates the reciprocal of the condition number (in the
     *  1-norm) of a real symmetric positive definite matrix using the
     *  Cholesky factorization A = U**T*U or A = L*L**T computed by DPOTRF.
     *
     *  An estimate is obtained for norm(inv(A)), and the reciprocal of the
     *  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
     *          The triangular factor U or L from the Cholesky factorization
     *          A = U**T*U or A = L*L**T, as computed by DPOTRF.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  ANORM   (input) DOUBLE PRECISION
     *          The 1-norm (or infinity-norm) of the symmetric matrix A.
     *
     *  RCOND   (output) DOUBLE PRECISION
     *          The reciprocal of the condition number of the matrix A,
     *          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
     *          estimate of the 1-norm of inv(A) computed in this routine.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param a
     * @param lda
     * @param normA
     */
    public static double dpocon(Lapack la, TUpLo uplo, int n, double[] a, int lda, double normA) {
        intW info = new intW(0);
        doubleW rcond = new doubleW(0.0);
        // TODO: remaining params sanity check
        double[] work = new double[3 * n];
        int[] iwork = new int[n];
        la.dpocon(uplo.val(), n, a, lda, normA, rcond, work, iwork, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        return rcond.val;
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPOSV computes the solution to a real system of linear equations
     *     A * X = B,
     *  where A is an N-by-N symmetric positive definite matrix and X and B
     *  are N-by-NRHS matrices.
     *
     *  The Cholesky decomposition is used to factor A as
     *     A = U**T* U,  if UPLO = 'U', or
     *     A = L * L**T,  if UPLO = 'L',
     *  where U is an upper triangular matrix and L is a lower triangular
     *  matrix.  The factored form of A is then used to solve the system of
     *  equations A * X = B.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The number of linear equations, i.e., the order of the
     *          matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
     *          N-by-N upper triangular part of A contains the upper
     *          triangular part of the matrix A, and the strictly lower
     *          triangular part of A is not referenced.  If UPLO = 'L', the
     *          leading N-by-N lower triangular part of A contains the lower
     *          triangular part of the matrix A, and the strictly upper
     *          triangular part of A is not referenced.
     *
     *          On exit, if INFO = 0, the factor U or L from the Cholesky
     *          factorization A = U**T*U or A = L*L**T.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the N-by-NRHS right hand side matrix B.
     *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the leading minor of order i of A is not
     *                positive definite, so the factorization could not be
     *                completed, and the solution has not been computed.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param rhsCount
     * @param a
     * @param lda
     * @param b
     * @param ldb
     */
    public static void dposv(Lapack la, TUpLo uplo, int n, int rhsCount, double[] a, int lda, double[] b, int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        la.dposv(uplo.val(), n, rhsCount, a, lda, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("The leading minor of order " + info.val
                        + " of A is not positive-definite. The factorization could not be completed "
                        + " and the solution has not been computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPOTRF computes the Cholesky factorization of a real symmetric
     *  positive definite matrix A.
     *
     *  The factorization has the form
     *     A = U**T * U,  if UPLO = 'U', or
     *     A = L  * L**T,  if UPLO = 'L',
     *  where U is an upper triangular matrix and L is lower triangular.
     *
     *  This is the block version of the algorithm, calling Level 3 BLAS.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
     *          N-by-N upper triangular part of A contains the upper
     *          triangular part of the matrix A, and the strictly lower
     *          triangular part of A is not referenced.  If UPLO = 'L', the
     *          leading N-by-N lower triangular part of A contains the lower
     *          triangular part of the matrix A, and the strictly upper
     *          triangular part of A is not referenced.
     *
     *          On exit, if INFO = 0, the factor U or L from the Cholesky
     *          factorization A = U**T*U or A = L*L**T.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the leading minor of order i is not
     *                positive definite, and the factorization could not be
     *                completed.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param a
     * @param lda
     */
    public static void dpotrf(Lapack la, TUpLo uplo, int n, double[] a, int lda) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        la.dpotrf(uplo.val(), n, a, lda, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("The leading minor of order " + info.val
                        + " is not positive-definite. The factorization could not be completed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPOTRS solves a system of linear equations A*X = B with a symmetric
     *  positive definite matrix A using the Cholesky factorization
     *  A = U**T*U or A = L*L**T computed by DPOTRF.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
     *          The triangular factor U or L from the Cholesky factorization
     *          A = U**T*U or A = L*L**T, as computed by DPOTRF.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the right hand side matrix B.
     *          On exit, the solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param rhsCount
     * @param a
     * @param lda
     * @param b
     * @param ldb
     */
    public static void dpotrs(Lapack la, TUpLo uplo, int n, int rhsCount, double[] a, int lda, double[] b, int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        la.dpotrs(uplo.val(), n, rhsCount, a, lda, b, ldb, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPPCON estimates the reciprocal of the condition number (in the
     *  1-norm) of a real symmetric positive definite packed matrix using
     *  the Cholesky factorization A = U**T*U or A = L*L**T computed by
     *  DPPTRF.
     *
     *  An estimate is obtained for norm(inv(A)), and the reciprocal of the
     *  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
     *          The triangular factor U or L from the Cholesky factorization
     *          A = U**T*U or A = L*L**T, packed columnwise in a linear
     *          array.  The j-th column of U or L is stored in the array AP
     *          as follows:
     *          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
     *          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
     *
     *  ANORM   (input) DOUBLE PRECISION
     *          The 1-norm (or infinity-norm) of the symmetric matrix A.
     *
     *  RCOND   (output) DOUBLE PRECISION
     *          The reciprocal of the condition number of the matrix A,
     *          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
     *          estimate of the 1-norm of inv(A) computed in this routine.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param ap
     * @param normA
     */
    public static double dppcon(Lapack la, TUpLo uplo, int n, double[] ap, double normA) {
        intW info = new intW(0);
        doubleW rcond = new doubleW(0.0);
        // TODO: remaining params sanity check
        double[] work = new double[3 * n];
        int[] iwork = new int[n];
        la.dppcon(uplo.val(), n, ap, normA, rcond, work, iwork, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        return rcond.val;
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPPSV computes the solution to a real system of linear equations
     *     A * X = B,
     *  where A is an N-by-N symmetric positive definite matrix stored in
     *  packed format and X and B are N-by-NRHS matrices.
     *
     *  The Cholesky decomposition is used to factor A as
     *     A = U**T* U,  if UPLO = 'U', or
     *     A = L * L**T,  if UPLO = 'L',
     *  where U is an upper triangular matrix and L is a lower triangular
     *  matrix.  The factored form of A is then used to solve the system of
     *  equations A * X = B.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The number of linear equations, i.e., the order of the
     *          matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
     *          On entry, the upper or lower triangle of the symmetric matrix
     *          A, packed columnwise in a linear array.  The j-th column of A
     *          is stored in the array AP as follows:
     *          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
     *          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
     *          See below for further details.
     *
     *          On exit, if INFO = 0, the factor U or L from the Cholesky
     *          factorization A = U**T*U or A = L*L**T, in the same storage
     *          format as A.
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the N-by-NRHS right hand side matrix B.
     *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the leading minor of order i of A is not
     *                positive definite, so the factorization could not be
     *                completed, and the solution has not been computed.
     *
     *  Further Details
     *  ===============
     *
     *  The packed storage scheme is illustrated by the following example
     *  when N = 4, UPLO = 'U':
     *
     *  Two-dimensional storage of the symmetric matrix A:
     *
     *     a11 a12 a13 a14
     *         a22 a23 a24
     *             a33 a34     (aij = conjg(aji))
     *                 a44
     *
     *  Packed storage of the upper triangle of A:
     *
     *  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param rhsCount
     * @param ap
     * @param b
     * @param ldb
     */
    public static void dppsv(Lapack la, TUpLo uplo, int n, int rhsCount, double[] ap, double[] b, int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        la.dppsv(uplo.val(), n, rhsCount, ap, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("The leading minor of order " + info.val
                        + " of A is not positive-definite. The factorization could not be completed "
                        + " and the solution has not been computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPPTRF computes the Cholesky factorization of a real symmetric
     *  positive definite matrix A stored in packed format.
     *
     *  The factorization has the form
     *     A = U**T * U,  if UPLO = 'U', or
     *     A = L  * L**T,  if UPLO = 'L',
     *  where U is an upper triangular matrix and L is lower triangular.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
     *          On entry, the upper or lower triangle of the symmetric matrix
     *          A, packed columnwise in a linear array.  The j-th column of A
     *          is stored in the array AP as follows:
     *          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
     *          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
     *          See below for further details.
     *
     *          On exit, if INFO = 0, the triangular factor U or L from the
     *          Cholesky factorization A = U**T*U or A = L*L**T, in the same
     *          storage format as A.
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the leading minor of order i is not
     *                positive definite, and the factorization could not be
     *                completed.
     *
     *  Further Details
     *  ======= =======
     *
     *  The packed storage scheme is illustrated by the following example
     *  when N = 4, UPLO = 'U':
     *
     *  Two-dimensional storage of the symmetric matrix A:
     *
     *     a11 a12 a13 a14
     *         a22 a23 a24
     *             a33 a34     (aij = aji)
     *                 a44
     *
     *  Packed storage of the upper triangle of A:
     *
     *  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param ap
     */
    public static void dpptrf(Lapack la, TUpLo uplo, int n, double[] ap) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        la.dpptrf(uplo.val(), n, ap, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("The leading minor of order " + info.val
                        + " is not positive-definite. The factorization could not be completed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPPTRS solves a system of linear equations A*X = B with a symmetric
     *  positive definite matrix A in packed storage using the Cholesky
     *  factorization A = U**T*U or A = L*L**T computed by DPPTRF.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
     *          The triangular factor U or L from the Cholesky factorization
     *          A = U**T*U or A = L*L**T, packed columnwise in a linear
     *          array.  The j-th column of U or L is stored in the array AP
     *          as follows:
     *          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
     *          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the right hand side matrix B.
     *          On exit, the solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param rhsCount
     * @param ap
     * @param b
     * @param ldb
     */
    public static void dpptrs(Lapack la, TUpLo uplo, int n, int rhsCount, double[] ap, double[] b, int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        la.dpptrs(uplo.val(), n, rhsCount, ap, b, ldb, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DPTSV computes the solution to a real system of linear equations
     *  A*X = B, where A is an N-by-N symmetric positive definite tridiagonal
     *  matrix, and X and B are N-by-NRHS matrices.
     *
     *  A is factored as A = L*D*L**T, and the factored form of A is then
     *  used to solve the system of equations.
     *
     *  Arguments
     *  =========
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  D       (input/output) DOUBLE PRECISION array, dimension (N)
     *          On entry, the n diagonal elements of the tridiagonal matrix
     *          A.  On exit, the n diagonal elements of the diagonal matrix
     *          D from the factorization A = L*D*L**T.
     *
     *  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
     *          On entry, the (n-1) subdiagonal elements of the tridiagonal
     *          matrix A.  On exit, the (n-1) subdiagonal elements of the
     *          unit bidiagonal factor L from the L*D*L**T factorization of
     *          A.  (E can also be regarded as the superdiagonal of the unit
     *          bidiagonal factor U from the U**T*D*U factorization of A.)
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the N-by-NRHS right hand side matrix B.
     *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the leading minor of order i is not
     *                positive definite, and the solution has not been
     *                computed.  The factorization has not been completed
     *                unless i = N.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param n
     * @param rhsCount
     * @param d
     * @param e
     * @param b
     * @param ldb
     */
    public static void dptsv(Lapack la, int n, int rhsCount, double[] d, double[] e, double[] b, int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        la.dptsv(n, rhsCount, d, e, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("The leading minor of order " + info.val
                        + " is not positive-definite. Solution could not be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSBEVD computes all the eigenvalues and, optionally, eigenvectors of
     *  a real symmetric band matrix A. If eigenvectors are desired, it uses
     *  a divide and conquer algorithm.
     *
     *  The divide and conquer algorithm makes very mild assumptions about
     *  floating point arithmetic. It will work on machines with a guard
     *  digit in add/subtract, or on those binary machines without guard
     *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
     *  without guard digits, but we know of none.
     *
     *  Arguments
     *  =========
     *
     *  JOBZ    (input) CHARACTER*1
     *          = 'N':  Compute eigenvalues only;
     *          = 'V':  Compute eigenvalues and eigenvectors.
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  KD      (input) INTEGER
     *          The number of superdiagonals of the matrix A if UPLO = 'U',
     *          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
     *
     *  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)
     *          On entry, the upper or lower triangle of the symmetric band
     *          matrix A, stored in the first KD+1 rows of the array.  The
     *          j-th column of A is stored in the j-th column of the array AB
     *          as follows:
     *          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
     *          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
     *
     *          On exit, AB is overwritten by values generated during the
     *          reduction to tridiagonal form.  If UPLO = 'U', the first
     *          superdiagonal and the diagonal of the tridiagonal matrix T
     *          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
     *          the diagonal and first subdiagonal of T are returned in the
     *          first two rows of AB.
     *
     *  LDAB    (input) INTEGER
     *          The leading dimension of the array AB.  LDAB >= KD + 1.
     *
     *  W       (output) DOUBLE PRECISION array, dimension (N)
     *          If INFO = 0, the eigenvalues in ascending order.
     *
     *  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
     *          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
     *          eigenvectors of the matrix A, with the i-th column of Z
     *          holding the eigenvector associated with W(i).
     *          If JOBZ = 'N', then Z is not referenced.
     *
     *  LDZ     (input) INTEGER
     *          The leading dimension of the array Z.  LDZ >= 1, and if
     *          JOBZ = 'V', LDZ >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the algorithm failed to converge; i
     *                off-diagonal elements of an intermediate tridiagonal
     *                form did not converge to zero.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param jobz
     * @param uplo
     * @param n
     * @param diagCount
     * @param ab
     * @param w
     * @param z
     * @param ldz
     */
    public static void dsbevd(Lapack la, TEigJob jobz, TUpLo uplo, int n, int diagCount, double[] ab, double[] w,
            double[] z, int ldz) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        double[] work = new double[1];
        int[] iwork = new int[1];
        int ldab = Math.max(1, diagCount + 1);
        la.dsbevd(jobz.val(), uplo.val(), n, diagCount, new double[0], ldab, new double[0], new double[0], ldz, work,
                -1, iwork, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        iwork = new int[iwork[0]];
        la.dsbevd(jobz.val(), uplo.val(), n, diagCount, ab, ldab, w, z, ldz, work, work.length, iwork, iwork.length,
                info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new NotConvergedException("Failed to converge");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSPEVD computes all the eigenvalues and, optionally, eigenvectors
     *  of a real symmetric matrix A in packed storage. If eigenvectors are
     *  desired, it uses a divide and conquer algorithm.
     *
     *  The divide and conquer algorithm makes very mild assumptions about
     *  floating point arithmetic. It will work on machines with a guard
     *  digit in add/subtract, or on those binary machines without guard
     *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
     *  without guard digits, but we know of none.
     *
     *  Arguments
     *  =========
     *
     *  JOBZ    (input) CHARACTER*1
     *          = 'N':  Compute eigenvalues only;
     *          = 'V':  Compute eigenvalues and eigenvectors.
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
     *          On entry, the upper or lower triangle of the symmetric matrix
     *          A, packed columnwise in a linear array.  The j-th column of A
     *          is stored in the array AP as follows:
     *          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
     *          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
     *
     *          On exit, AP is overwritten by values generated during the
     *          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
     *          and first superdiagonal of the tridiagonal matrix T overwrite
     *          the corresponding elements of A, and if UPLO = 'L', the
     *          diagonal and first subdiagonal of T overwrite the
     *          corresponding elements of A.
     *
     *  W       (output) DOUBLE PRECISION array, dimension (N)
     *          If INFO = 0, the eigenvalues in ascending order.
     *
     *  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
     *          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
     *          eigenvectors of the matrix A, with the i-th column of Z
     *          holding the eigenvector associated with W(i).
     *          If JOBZ = 'N', then Z is not referenced.
     *
     *  LDZ     (input) INTEGER
     *          The leading dimension of the array Z.  LDZ >= 1, and if
     *          JOBZ = 'V', LDZ >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value.
     *          > 0:  if INFO = i, the algorithm failed to converge; i
     *                off-diagonal elements of an intermediate tridiagonal
     *                form did not converge to zero.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param jobz
     * @param uplo
     * @param n
     * @param ap
     * @param w
     * @param z
     * @param ldz
     */
    public static void dspevd(Lapack la, TEigJob jobz, TUpLo uplo, int n, double[] ap, double[] w, double[] z,
            int ldz) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        double[] work = new double[1];
        int[] iwork = new int[1];
        la.dspevd(jobz.val(), uplo.val(), n, new double[0], new double[0], new double[0], ldz, work, -1, iwork, -1,
                info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        iwork = new int[iwork[0]];
        la.dspevd(jobz.val(), uplo.val(), n, ap, w, z, ldz, work, work.length, iwork, iwork.length, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new NotConvergedException("Failed to converge");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSPSV computes the solution to a real system of linear equations
     *     A * X = B,
     *  where A is an N-by-N symmetric matrix stored in packed format and X
     *  and B are N-by-NRHS matrices.
     *
     *  The diagonal pivoting method is used to factor A as
     *     A = U * D * U**T,  if UPLO = 'U', or
     *     A = L * D * L**T,  if UPLO = 'L',
     *  where U (or L) is a product of permutation and unit upper (lower)
     *  triangular matrices, D is symmetric and block diagonal with 1-by-1
     *  and 2-by-2 diagonal blocks.  The factored form of A is then used to
     *  solve the system of equations A * X = B.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The number of linear equations, i.e., the order of the
     *          matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
     *          On entry, the upper or lower triangle of the symmetric matrix
     *          A, packed columnwise in a linear array.  The j-th column of A
     *          is stored in the array AP as follows:
     *          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
     *          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
     *          See below for further details.
     *
     *          On exit, the block diagonal matrix D and the multipliers used
     *          to obtain the factor U or L from the factorization
     *          A = U*D*U**T or A = L*D*L**T as computed by DSPTRF, stored as
     *          a packed triangular matrix in the same storage format as A.
     *
     *  IPIV    (output) INTEGER array, dimension (N)
     *          Details of the interchanges and the block structure of D, as
     *          determined by DSPTRF.  If IPIV(k) > 0, then rows and columns
     *          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
     *          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
     *          then rows and columns k-1 and -IPIV(k) were interchanged and
     *          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
     *          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
     *          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
     *          diagonal block.
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the N-by-NRHS right hand side matrix B.
     *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
     *                has been completed, but the block diagonal matrix D is
     *                exactly singular, so the solution could not be
     *                computed.
     *
     *  Further Details
     *  ===============
     *
     *  The packed storage scheme is illustrated by the following example
     *  when N = 4, UPLO = 'U':
     *
     *  Two-dimensional storage of the symmetric matrix A:
     *
     *     a11 a12 a13 a14
     *         a22 a23 a24
     *             a33 a34     (aij = aji)
     *                 a44
     *
     *  Packed storage of the upper triangle of A:
     *
     *  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param rhsCount
     * @param ap
     * @param indices
     * @param b
     * @param ldb
     */
    public static void dspsv(Lapack la, TUpLo uplo, int n, int rhsCount, double[] ap, int[] indices, double[] b,
            int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        la.dspsv(uplo.val(), n, rhsCount, ap, indices, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "The block diagonal matrix D is exactly singular. Solution could not be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSTEVR computes selected eigenvalues and, optionally, eigenvectors
     *  of a real symmetric tridiagonal matrix T.  Eigenvalues and
     *  eigenvectors can be selected by specifying either a range of values
     *  or a range of indices for the desired eigenvalues.
     *
     *  Whenever possible, DSTEVR calls DSTEMR to compute the
     *  eigenspectrum using Relatively Robust Representations.  DSTEMR
     *  computes eigenvalues by the dqds algorithm, while orthogonal
     *  eigenvectors are computed from various "good" L D L^T representations
     *  (also known as Relatively Robust Representations). Gram-Schmidt
     *  orthogonalization is avoided as far as possible. More specifically,
     *  the various steps of the algorithm are as follows. For the i-th
     *  unreduced block of T,
     *     (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
     *          is a relatively robust representation,
     *     (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high
     *         relative accuracy by the dqds algorithm,
     *     (c) If there is a cluster of close eigenvalues, "choose" sigma_i
     *         close to the cluster, and go to step (a),
     *     (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,
     *         compute the corresponding eigenvector by forming a
     *         rank-revealing twisted factorization.
     *  The desired accuracy of the output can be specified by the input
     *  parameter ABSTOL.
     *
     *  For more details, see "A new O(n^2) algorithm for the symmetric
     *  tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
     *  Computer Science Division Technical Report No. UCB//CSD-97-971,
     *  UC Berkeley, May 1997.
     *
     *
     *  Note 1 : DSTEVR calls DSTEMR when the full spectrum is requested
     *  on machines which conform to the ieee-754 floating point standard.
     *  DSTEVR calls DSTEBZ and DSTEIN on non-ieee machines and
     *  when partial spectrum requests are made.
     *
     *  Normal execution of DSTEMR may create NaNs and infinities and
     *  hence may abort due to a floating point exception in environments
     *  which do not handle NaNs and infinities in the ieee standard default
     *  manner.
     *
     *  Arguments
     *  =========
     *
     *  JOBZ    (input) CHARACTER*1
     *          = 'N':  Compute eigenvalues only;
     *          = 'V':  Compute eigenvalues and eigenvectors.
     *
     *  RANGE   (input) CHARACTER*1
     *          = 'A': all eigenvalues will be found.
     *          = 'V': all eigenvalues in the half-open interval (VL,VU]
     *                 will be found.
     *          = 'I': the IL-th through IU-th eigenvalues will be found.
     ********** For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
     ********** DSTEIN are called
     *
     *  N       (input) INTEGER
     *          The order of the matrix.  N >= 0.
     *
     *  D       (input/output) DOUBLE PRECISION array, dimension (N)
     *          On entry, the n diagonal elements of the tridiagonal matrix
     *          A.
     *          On exit, D may be multiplied by a constant factor chosen
     *          to avoid over/underflow in computing the eigenvalues.
     *
     *  E       (input/output) DOUBLE PRECISION array, dimension (max(1,N-1))
     *          On entry, the (n-1) subdiagonal elements of the tridiagonal
     *          matrix A in elements 1 to N-1 of E.
     *          On exit, E may be multiplied by a constant factor chosen
     *          to avoid over/underflow in computing the eigenvalues.
     *
     *  VL      (input) DOUBLE PRECISION
     *  VU      (input) DOUBLE PRECISION
     *          If RANGE='V', the lower and upper bounds of the interval to
     *          be searched for eigenvalues. VL < VU.
     *          Not referenced if RANGE = 'A' or 'I'.
     *
     *  IL      (input) INTEGER
     *  IU      (input) INTEGER
     *          If RANGE='I', the indices (in ascending order) of the
     *          smallest and largest eigenvalues to be returned.
     *          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
     *          Not referenced if RANGE = 'A' or 'V'.
     *
     *  ABSTOL  (input) DOUBLE PRECISION
     *          The absolute error tolerance for the eigenvalues.
     *          An approximate eigenvalue is accepted as converged
     *          when it is determined to lie in an interval [a,b]
     *          of width less than or equal to
     *
     *                  ABSTOL + EPS *   max( |a|,|b| ) ,
     *
     *          where EPS is the machine precision.  If ABSTOL is less than
     *          or equal to zero, then  EPS*|T|  will be used in its place,
     *          where |T| is the 1-norm of the tridiagonal matrix obtained
     *          by reducing A to tridiagonal form.
     *
     *          See "Computing Small Singular Values of Bidiagonal Matrices
     *          with Guaranteed High Relative Accuracy," by Demmel and
     *          Kahan, LAPACK Working Note #3.
     *
     *          If high relative accuracy is important, set ABSTOL to
     *          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that
     *          eigenvalues are computed to high relative accuracy when
     *          possible in future releases.  The current code does not
     *          make any guarantees about high relative accuracy, but
     *          future releases will. See J. Barlow and J. Demmel,
     *          "Computing Accurate Eigensystems of Scaled Diagonally
     *          Dominant Matrices", LAPACK Working Note #7, for a discussion
     *          of which matrices define their eigenvalues to high relative
     *          accuracy.
     *
     *  M       (output) INTEGER
     *          The total number of eigenvalues found.  0 <= M <= N.
     *          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
     *
     *  W       (output) DOUBLE PRECISION array, dimension (N)
     *          The first M elements contain the selected eigenvalues in
     *          ascending order.
     *
     *  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )
     *          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
     *          contain the orthonormal eigenvectors of the matrix A
     *          corresponding to the selected eigenvalues, with the i-th
     *          column of Z holding the eigenvector associated with W(i).
     *          Note: the user must ensure that at least max(1,M) columns are
     *          supplied in the array Z; if RANGE = 'V', the exact value of M
     *          is not known in advance and an upper bound must be used.
     *
     *  LDZ     (input) INTEGER
     *          The leading dimension of the array Z.  LDZ >= 1, and if
     *          JOBZ = 'V', LDZ >= max(1,N).
     *
     *  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )
     *          The support of the eigenvectors in Z, i.e., the indices
     *          indicating the nonzero elements in Z. The i-th eigenvector
     *          is nonzero only in elements ISUPPZ( 2*i-1 ) through
     *          ISUPPZ( 2*i ).
     ********** Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  Internal error
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param jobz
     * @param range
     * @param n
     * @param d
     * @param e
     * @param vLower
     * @param vUpper
     * @param iLower
     * @param iUpper
     * @param abstol
     * @param w
     * @param z
     * @param ldz
     * @param supportZ
     */
    public static int dstevr(Lapack la, TEigJob jobz, TRange range, int n, double[] d, double[] e, double vLower,
            double vUpper, int iLower, int iUpper, double abstol, double[] w, double[] z, int ldz, int[] supportZ) {
        intW info = new intW(0);
        intW m = new intW(0);
        // TODO: remaining params sanity check
        double[] work = new double[1];
        int[] iwork = new int[1];
        la.dstevr(jobz.val(), range.val(), n, new double[0], new double[0], vLower, vUpper, iLower, iUpper, abstol, m,
                new double[0], new double[0], ldz, new int[0], work, -1, iwork, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        iwork = new int[iwork[0]];
        la.dstevr(jobz.val(), range.val(), n, d, e, vLower, vUpper, iLower, iUpper, abstol, m, w, z, ldz, supportZ,
                work, work.length, iwork, iwork.length, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new IllegalStateException("Internal error");
            }
        }
        return m.val;
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSYEVR computes selected eigenvalues and, optionally, eigenvectors
     *  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
     *  selected by specifying either a range of values or a range of
     *  indices for the desired eigenvalues.
     *
     *  DSYEVR first reduces the matrix A to tridiagonal form T with a call
     *  to DSYTRD.  Then, whenever possible, DSYEVR calls DSTEMR to compute
     *  the eigenspectrum using Relatively Robust Representations.  DSTEMR
     *  computes eigenvalues by the dqds algorithm, while orthogonal
     *  eigenvectors are computed from various "good" L D L^T representations
     *  (also known as Relatively Robust Representations). Gram-Schmidt
     *  orthogonalization is avoided as far as possible. More specifically,
     *  the various steps of the algorithm are as follows.
     *
     *  For each unreduced block (submatrix) of T,
     *     (a) Compute T - sigma I  = L D L^T, so that L and D
     *         define all the wanted eigenvalues to high relative accuracy.
     *         This means that small relative changes in the entries of D and
     *         cause only small relative changes in the eigenvalues and
     *         eigenvectors. The standard (unfactored) representation of the
     *         tridiagonal matrix T does not have this property in general.
     *     (b) Compute the eigenvalues to suitable accuracy.
     *         If the eigenvectors are desired, the algorithm attains full
     *         accuracy of the computed eigenvalues only right before
     *         the corresponding vectors have to be computed, see steps c) an
     *     (c) For each cluster of close eigenvalues, select a new
     *         shift close to the cluster, find a new factorization, and refi
     *         the shifted eigenvalues to suitable accuracy.
     *     (d) For each eigenvalue with a large enough relative separation co
     *         the corresponding eigenvector by forming a rank revealing twis
     *         factorization. Go back to (c) for any clusters that remain.
     *
     *  The desired accuracy of the output can be specified by the input
     *  parameter ABSTOL.
     *
     *  Note 1 : DSYEVR calls DSTEMR when the full spectrum is requested
     *  on machines which conform to the ieee-754 floating point standard.
     *  DSYEVR calls DSTEBZ and SSTEIN on non-ieee machines and
     *  when partial spectrum requests are made.
     *
     *  Normal execution of DSTEMR may create NaNs and infinities and
     *  hence may abort due to a floating point exception in environments
     *  which do not handle NaNs and infinities in the ieee standard default
     *  manner.
     *
     *  Arguments
     *  =========
     *
     *  JOBZ    (input) CHARACTER*1
     *          = 'N':  Compute eigenvalues only;
     *          = 'V':  Compute eigenvalues and eigenvectors.
     *
     *  RANGE   (input) CHARACTER*1
     *          = 'A': all eigenvalues will be found.
     *          = 'V': all eigenvalues in the half-open interval (VL,VU]
     *                 will be found.
     *          = 'I': the IL-th through IU-th eigenvalues will be found.
     ********** For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
     ********** DSTEIN are called
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
     *          On entry, the symmetric matrix A.  If UPLO = 'U', the
     *          leading N-by-N upper triangular part of A contains the
     *          upper triangular part of the matrix A.  If UPLO = 'L',
     *          the leading N-by-N lower triangular part of A contains
     *          the lower triangular part of the matrix A.
     *          On exit, the lower triangle (if UPLO='L') or the upper
     *          triangle (if UPLO='U') of A, including the diagonal, is
     *          destroyed.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  VL      (input) DOUBLE PRECISION
     *  VU      (input) DOUBLE PRECISION
     *          If RANGE='V', the lower and upper bounds of the interval to
     *          be searched for eigenvalues. VL < VU.
     *          Not referenced if RANGE = 'A' or 'I'.
     *
     *  IL      (input) INTEGER
     *  IU      (input) INTEGER
     *          If RANGE='I', the indices (in ascending order) of the
     *          smallest and largest eigenvalues to be returned.
     *          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
     *          Not referenced if RANGE = 'A' or 'V'.
     *
     *  ABSTOL  (input) DOUBLE PRECISION
     *          The absolute error tolerance for the eigenvalues.
     *          An approximate eigenvalue is accepted as converged
     *          when it is determined to lie in an interval [a,b]
     *          of width less than or equal to
     *
     *                  ABSTOL + EPS *   max( |a|,|b| ) ,
     *
     *          where EPS is the machine precision.  If ABSTOL is less than
     *          or equal to zero, then  EPS*|T|  will be used in its place,
     *          where |T| is the 1-norm of the tridiagonal matrix obtained
     *          by reducing A to tridiagonal form.
     *
     *          See "Computing Small Singular Values of Bidiagonal Matrices
     *          with Guaranteed High Relative Accuracy," by Demmel and
     *          Kahan, LAPACK Working Note #3.
     *
     *          If high relative accuracy is important, set ABSTOL to
     *          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that
     *          eigenvalues are computed to high relative accuracy when
     *          possible in future releases.  The current code does not
     *          make any guarantees about high relative accuracy, but
     *          future releases will. See J. Barlow and J. Demmel,
     *          "Computing Accurate Eigensystems of Scaled Diagonally
     *          Dominant Matrices", LAPACK Working Note #7, for a discussion
     *          of which matrices define their eigenvalues to high relative
     *          accuracy.
     *
     *  M       (output) INTEGER
     *          The total number of eigenvalues found.  0 <= M <= N.
     *          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
     *
     *  W       (output) DOUBLE PRECISION array, dimension (N)
     *          The first M elements contain the selected eigenvalues in
     *          ascending order.
     *
     *  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
     *          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
     *          contain the orthonormal eigenvectors of the matrix A
     *          corresponding to the selected eigenvalues, with the i-th
     *          column of Z holding the eigenvector associated with W(i).
     *          If JOBZ = 'N', then Z is not referenced.
     *          Note: the user must ensure that at least max(1,M) columns are
     *          supplied in the array Z; if RANGE = 'V', the exact value of M
     *          is not known in advance and an upper bound must be used.
     *          Supplying N columns is always safe.
     *
     *  LDZ     (input) INTEGER
     *          The leading dimension of the array Z.  LDZ >= 1, and if
     *          JOBZ = 'V', LDZ >= max(1,N).
     *
     *  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )
     *          The support of the eigenvectors in Z, i.e., the indices
     *          indicating the nonzero elements in Z. The i-th eigenvector
     *          is nonzero only in elements ISUPPZ( 2*i-1 ) through
     *          ISUPPZ( 2*i ).
     ********** Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  Internal error
     *
     * =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param jobz
     * @param range
     * @param uplo
     * @param n
     * @param a
     * @param lda
     * @param vLower
     * @param vUpper
     * @param iLower
     * @param iUpper
     * @param abstol
     * @param w
     * @param z
     * @param ldz
     * @param supportZ
     */
    public static int dsyevr(Lapack la, TEigJob jobz, TRange range, TUpLo uplo, int n, double[] a, int lda,
            double vLower, double vUpper, int iLower, int iUpper, double abstol, double[] w, double[] z, int ldz,
            int[] supportZ) {
        intW info = new intW(0);
        intW m = new intW(0);
        // TODO: remaining params sanity check
        double[] work = new double[1];
        int[] iwork = new int[1];
        la.dsyevr(jobz.val(), range.val(), uplo.val(), n, new double[0], lda, vLower, vUpper, iLower, iUpper, abstol, m,
                new double[0], new double[0], ldz, new int[0], work, -1, iwork, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        iwork = new int[iwork[0]];
        la.dsyevr(jobz.val(), range.val(), uplo.val(), n, a, lda, vLower, vUpper, iLower, iUpper, abstol, m, w, z, ldz,
                supportZ, work, work.length, iwork, iwork.length, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new IllegalStateException("Internal error");
            }
        }
        return m.val;
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSYGVD computes all the eigenvalues, and optionally, the eigenvectors
     *  of a real generalized symmetric-definite eigenproblem, of the form
     *  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
     *  B are assumed to be symmetric and B is also positive definite.
     *  If eigenvectors are desired, it uses a divide and conquer algorithm.
     *
     *  The divide and conquer algorithm makes very mild assumptions about
     *  floating point arithmetic. It will work on machines with a guard
     *  digit in add/subtract, or on those binary machines without guard
     *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
     *  without guard digits, but we know of none.
     *
     *  Arguments
     *  =========
     *
     *  ITYPE   (input) INTEGER
     *          Specifies the problem type to be solved:
     *          = 1:  A*x = (lambda)*B*x
     *          = 2:  A*B*x = (lambda)*x
     *          = 3:  B*A*x = (lambda)*x
     *
     *  JOBZ    (input) CHARACTER*1
     *          = 'N':  Compute eigenvalues only;
     *          = 'V':  Compute eigenvalues and eigenvectors.
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangles of A and B are stored;
     *          = 'L':  Lower triangles of A and B are stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrices A and B.  N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
     *          On entry, the symmetric matrix A.  If UPLO = 'U', the
     *          leading N-by-N upper triangular part of A contains the
     *          upper triangular part of the matrix A.  If UPLO = 'L',
     *          the leading N-by-N lower triangular part of A contains
     *          the lower triangular part of the matrix A.
     *
     *          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
     *          matrix Z of eigenvectors.  The eigenvectors are normalized
     *          as follows:
     *          if ITYPE = 1 or 2, Z**T*B*Z = I;
     *          if ITYPE = 3, Z**T*inv(B)*Z = I.
     *          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
     *          or the lower triangle (if UPLO='L') of A, including the
     *          diagonal, is destroyed.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
     *          On entry, the symmetric matrix B.  If UPLO = 'U', the
     *          leading N-by-N upper triangular part of B contains the
     *          upper triangular part of the matrix B.  If UPLO = 'L',
     *          the leading N-by-N lower triangular part of B contains
     *          the lower triangular part of the matrix B.
     *
     *          On exit, if INFO <= N, the part of B containing the matrix is
     *          overwritten by the triangular factor U or L from the Cholesky
     *          factorization B = U**T*U or B = L*L**T.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  W       (output) DOUBLE PRECISION array, dimension (N)
     *          If INFO = 0, the eigenvalues in ascending order.
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  DPOTRF or DSYEVD returned an error code:
     *             <= N:  if INFO = i and JOBZ = 'N', then the algorithm
     *                    failed to converge; i off-diagonal elements of an
     *                    intermediate tridiagonal form did not converge to
     *                    zero;
     *                    if INFO = i and JOBZ = 'V', then the algorithm
     *                    failed to compute an eigenvalue while working on
     *                    the submatrix lying in rows and columns INFO/(N+1)
     *                    through mod(INFO,N+1);
     *             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
     *                    minor of order i of B is not positive definite.
     *                    The factorization of B could not be completed and
     *                    no eigenvalues or eigenvectors were computed.
     *
     * </code>
     * </pre>
     *
     * @param type
     * @param jobz
     * @param uplo
     * @param n
     * @param a
     * @param lda
     * @param b
     * @param ldb
     * @param w
     */
    public static void dsygvd(Lapack la, int type, TEigJob jobz, TUpLo uplo, int n, double[] a, int lda, double[] b,
            int ldb, double[] w) {
        if (type < 1 || type > 3) {
            throw new IllegalArgumentException("type : " + type + " must be from closed interval [1..3]");
        }
        // TODO: remaining params sanity check
        intW info = new intW(0);
        double[] work = new double[1];
        int[] iwork = new int[1];
        la.dsygvd(type, jobz.val(), uplo.val(), n, new double[0], lda, new double[0], ldb, new double[0], work, -1,
                iwork, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        iwork = new int[iwork[0]];
        la.dsygvd(type, jobz.val(), uplo.val(), n, a, lda, b, ldb, w, work, work.length, iwork, iwork.length, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                if (info.val > n) {
                    throw new ComputationTruncatedException("The leading minor of order " + info.val
                            + " of B is not positive-semidifinte. Factorization could not be completed.");
                } else {
                    throw new NotConvergedException("Failed to converge");
                }
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DSYSV computes the solution to a real system of linear equations
     *     A * X = B,
     *  where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     *  matrices.
     *
     *  The diagonal pivoting method is used to factor A as
     *     A = U * D * U**T,  if UPLO = 'U', or
     *     A = L * D * L**T,  if UPLO = 'L',
     *  where U (or L) is a product of permutation and unit upper (lower)
     *  triangular matrices, and D is symmetric and block diagonal with
     *  1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     *  used to solve the system of equations A * X = B.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The number of linear equations, i.e., the order of the
     *          matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
     *          N-by-N upper triangular part of A contains the upper
     *          triangular part of the matrix A, and the strictly lower
     *          triangular part of A is not referenced.  If UPLO = 'L', the
     *          leading N-by-N lower triangular part of A contains the lower
     *          triangular part of the matrix A, and the strictly upper
     *          triangular part of A is not referenced.
     *
     *          On exit, if INFO = 0, the block diagonal matrix D and the
     *          multipliers used to obtain the factor U or L from the
     *          factorization A = U*D*U**T or A = L*D*L**T as computed by
     *          DSYTRF.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  IPIV    (output) INTEGER array, dimension (N)
     *          Details of the interchanges and the block structure of D, as
     *          determined by DSYTRF.  If IPIV(k) > 0, then rows and columns
     *          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
     *          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
     *          then rows and columns k-1 and -IPIV(k) were interchanged and
     *          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
     *          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
     *          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
     *          diagonal block.
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the N-by-NRHS right hand side matrix B.
     *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0: successful exit
     *          < 0: if INFO = -i, the i-th argument had an illegal value
     *          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
     *               has been completed, but the block diagonal matrix D is
     *               exactly singular, so the solution could not be computed.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param n
     * @param rhsCount
     * @param a
     * @param lda
     * @param indices
     * @param b
     * @param ldb
     */
    public static void dsysv(Lapack la, TUpLo uplo, int n, int rhsCount, double[] a, int lda, int[] indices, double[] b,
            int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        double[] work = new double[1];
        la.dsysv(uplo.val(), n, rhsCount, new double[0], lda, new int[0], new double[0], ldb, work, -1, info);
        if (info.val != 0) {
            throwIAEPosition(info);
        }
        work = new double[(int) work[0]];
        la.dsysv(uplo.val(), n, rhsCount, a, lda, indices, b, ldb, work, work.length, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException(
                        "The block-diagonal matrix D is exactly singular. The solution couldn't be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DTBTRS solves a triangular system of the form
     *
     *     A * X = B  or  A**T * X = B,
     *
     *  where A is a triangular band matrix of order N, and B is an
     *  N-by NRHS matrix.  A check is made to verify that A is nonsingular.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  A is upper triangular;
     *          = 'L':  A is lower triangular.
     *
     *  TRANS   (input) CHARACTER*1
     *          Specifies the form the system of equations:
     *          = 'N':  A * X = B  (No transpose)
     *          = 'T':  A**T * X = B  (Transpose)
     *          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
     *
     *  DIAG    (input) CHARACTER*1
     *          = 'N':  A is non-unit triangular;
     *          = 'U':  A is unit triangular.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  KD      (input) INTEGER
     *          The number of superdiagonals or subdiagonals of the
     *          triangular band matrix A.  KD >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
     *          The upper or lower triangular band matrix A, stored in the
     *          first kd+1 rows of AB.  The j-th column of A is stored
     *          in the j-th column of the array AB as follows:
     *          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
     *          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
     *          If DIAG = 'U', the diagonal elements of A are not referenced
     *          and are assumed to be 1.
     *
     *  LDAB    (input) INTEGER
     *          The leading dimension of the array AB.  LDAB >= KD+1.
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the right hand side matrix B.
     *          On exit, if INFO = 0, the solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the i-th diagonal element of A is zero,
     *                indicating that the matrix is singular and the
     *                solutions X have not been computed.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param trans
     * @param diag
     * @param n
     * @param diagCount
     * @param rhsCount
     * @param ab
     * @param b
     * @param ldb
     */
    public static void dtbtrs(Lapack la, TUpLo uplo, TTrans trans, TDiag diag, int n, int diagCount, int rhsCount,
            double[] ab, double[] b, int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        int ldab = Math.max(1, diagCount + 1);
        la.dtbtrs(uplo.val(), trans.val(), diag.val(), n, diagCount, rhsCount, ab, ldab, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("Matrix A is singular. The solution couldn't be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DTPTRS solves a triangular system of the form
     *
     *     A * X = B  or  A**T * X = B,
     *
     *  where A is a triangular matrix of order N stored in packed format,
     *  and B is an N-by-NRHS matrix.  A check is made to verify that A is
     *  nonsingular.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  A is upper triangular;
     *          = 'L':  A is lower triangular.
     *
     *  TRANS   (input) CHARACTER*1
     *          Specifies the form of the system of equations:
     *          = 'N':  A * X = B  (No transpose)
     *          = 'T':  A**T * X = B  (Transpose)
     *          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
     *
     *  DIAG    (input) CHARACTER*1
     *          = 'N':  A is non-unit triangular;
     *          = 'U':  A is unit triangular.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
     *          The upper or lower triangular matrix A, packed columnwise in
     *          a linear array.  The j-th column of A is stored in the array
     *          AP as follows:
     *          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
     *          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the right hand side matrix B.
     *          On exit, if INFO = 0, the solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the i-th diagonal element of A is zero,
     *                indicating that the matrix is singular and the
     *                solutions X have not been computed.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param trans
     * @param diag
     * @param n
     * @param rhsCount
     * @param ap
     * @param b
     * @param ldb
     */
    public static void dtptrs(Lapack la, TUpLo uplo, TTrans trans, TDiag diag, int n, int rhsCount, double[] ap,
            double[] b, int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        la.dtptrs(uplo.val(), trans.val(), diag.val(), n, rhsCount, ap, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("Matrix A is singular. The solution couldn't be computed.");
            }
        }
    }

    /**
     * <pre>
     * <code>
     *
     *  Purpose
     *  =======
     *
     *  DTRTRS solves a triangular system of the form
     *
     *     A * X = B  or  A**T * X = B,
     *
     *  where A is a triangular matrix of order N, and B is an N-by-NRHS
     *  matrix.  A check is made to verify that A is nonsingular.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  A is upper triangular;
     *          = 'L':  A is lower triangular.
     *
     *  TRANS   (input) CHARACTER*1
     *          Specifies the form of the system of equations:
     *          = 'N':  A * X = B  (No transpose)
     *          = 'T':  A**T * X = B  (Transpose)
     *          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
     *
     *  DIAG    (input) CHARACTER*1
     *          = 'N':  A is non-unit triangular;
     *          = 'U':  A is unit triangular.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
     *          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
     *          upper triangular part of the array A contains the upper
     *          triangular matrix, and the strictly lower triangular part of
     *          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
     *          triangular part of the array A contains the lower triangular
     *          matrix, and the strictly upper triangular part of A is not
     *          referenced.  If DIAG = 'U', the diagonal elements of A are
     *          also not referenced and are assumed to be 1.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the right hand side matrix B.
     *          On exit, if INFO = 0, the solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0: if INFO = -i, the i-th argument had an illegal value
     *          > 0: if INFO = i, the i-th diagonal element of A is zero,
     *               indicating that the matrix is singular and the solutions
     *               X have not been computed.
     *
     *  =====================================================================
     *
     * </code>
     * </pre>
     *
     * @param uplo
     * @param trans
     * @param diag
     * @param n
     * @param rhsCount
     * @param a
     * @param lda
     * @param b
     * @param ldb
     */
    public static void dtrtrs(Lapack la, TUpLo uplo, TTrans trans, TDiag diag, int n, int rhsCount, double[] a, int lda,
            double[] b, int ldb) {
        intW info = new intW(0);
        // TODO: remaining params sanity check
        la.dtrtrs(uplo.val(), trans.val(), diag.val(), n, rhsCount, a, lda, b, ldb, info);
        if (info.val != 0) {
            if (info.val < 0) {
                throwIAEPosition(info);
            } else {
                throw new ComputationTruncatedException("Matrix A is singular. The solution couldn't be computed.");
            }
        }
    }

    // helper methods for argument checks

    private static void checkNonNegative(int value, String name) {
        if (value < 0) {
            throw new IllegalArgumentException("Parameter " + name + " must be non-negative (value = " + value + ")");
        }
    }

    private static void checkStrictlyPositive(int value, String name) {
        if (value <= 0) {
            throw new IllegalArgumentException(
                    "Parameter " + name + " must be strictly positive (value = " + value + ")");
        }
    }

    private static void checkValueAtLeast(int value, int minVal, String name) {
        if (value < Math.max(1, minVal)) {
            throw new IllegalArgumentException(
                    "Parameter " + name + " must be at least " + Math.max(1, minVal) + " (value = " + value + ")");
        }
    }

    private static void checkMinLen(double[] array, int minLen, String name) {
        if (array.length < minLen) {
            throw new IllegalArgumentException("Length of array '" + name + "' argument must be at least " + minLen
                    + " (length = " + array.length + ")");
        }
    }

    private static void checkMinLen(float[] array, int minLen, String name) {
        if (array.length < minLen) {
            throw new IllegalArgumentException("Length of array '" + name + "' argument must be at least " + minLen
                    + " (length = " + array.length + ")");
        }
    }

    private static void checkMinLen(int[] array, int minLen, String name) {
        if (array.length < minLen) {
            throw new IllegalArgumentException("Length of array '" + name + "' argument must be at least " + minLen
                    + " (length = " + array.length + ")");
        }
    }

    private static void throwIAEPosition(intW info) {
        throw new IllegalArgumentException("Illegal argument at position " + -info.val);
    }

    private static void throwIAEPosition(int info) {
        throw new IllegalArgumentException("Illegal argument at position " + -info);
    }

    protected PlainLapack() {
    }
}
