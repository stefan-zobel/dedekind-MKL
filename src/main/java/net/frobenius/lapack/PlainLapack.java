/*
 * Copyright 2019 Stefan Zobel
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
        // ldvt
        if (jobz == TSvdJob.ALL || (m >= n && jobz == TSvdJob.OVERWRITE)) {
            checkValueAtLeast(ldvt, n, "ldvt");
        } else if (jobz == TSvdJob.PART) {
            checkValueAtLeast(ldvt, Math.min(m, n), "ldvt");
        }
        // vt
        checkMinLen(vt, ldvt * n, "vt");

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
     * @param k1
     * @param k2
     * @param indices
     * @param incx
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

    private static void checkMinLen(int[] array, int minLen, String name) {
        throw new IllegalArgumentException("Length of array '" + name + "' argument must be at least " + minLen
                + " (length = " + array.length + ")");
    }

    private static void throwIAEPosition(intW info) {
        throw new IllegalArgumentException("Illegal argument at position " + -info.val);
    }

    protected PlainLapack() {
    }
}
