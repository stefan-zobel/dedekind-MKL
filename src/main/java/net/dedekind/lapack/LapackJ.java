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
package net.dedekind.lapack;

import org.netlib.lapack.*;
import org.netlib.util.doubleW;
import org.netlib.util.intW;

/**
 * {@code Netlib} implementation.
 */
public class LapackJ extends Lapack {

    @Override
    public final void dgbcon(String norm, int n, int kl, int ku, double[] ab, int abOffset, int ldab, int[] ipiv,
            int ipivOffset, double anorm, doubleW rcond, double[] work, int workOffset, int[] iwork, int iworkOffset,
            intW info) {
        Dgbcon.dgbcon(norm, n, kl, ku, ab, abOffset, ldab, ipiv, ipivOffset, anorm, rcond, work, workOffset, iwork,
                iworkOffset, info);
    }

    @Override
    public final void dgbsv(int n, int kl, int ku, int nrhs, double[] ab, int abOffset, int ldab, int[] ipiv,
            int ipivOffset, double[] b, int bOffset, int ldb, intW info) {
        Dgbsv.dgbsv(n, kl, ku, nrhs, ab, abOffset, ldab, ipiv, ipivOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void dgbtrf(int m, int n, int kl, int ku, double[] ab, int abOffset, int ldab, int[] ipiv,
            int ipivOffset, intW info) {
        Dgbtrf.dgbtrf(m, n, kl, ku, ab, abOffset, ldab, ipiv, ipivOffset, info);
    }

    @Override
    public final void dgbtrs(String trans, int n, int kl, int ku, int nrhs, double[] ab, int abOffset, int ldab,
            int[] ipiv, int ipivOffset, double[] b, int bOffset, int ldb, intW info) {
        Dgbtrs.dgbtrs(trans, n, kl, ku, nrhs, ab, abOffset, ldab, ipiv, ipivOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void dgecon(String norm, int n, double[] a, int aOffset, int lda, double anorm, doubleW rcond,
            double[] work, int workOffset, int[] iwork, int iworkOffset, intW info) {
        Dgecon.dgecon(norm, n, a, aOffset, lda, anorm, rcond, work, workOffset, iwork, iworkOffset, info);
    }

    @Override
    public final void dgeev(String jobvl, String jobvr, int n, double[] a, int aOffset, int lda, double[] wr,
            int wrOffset, double[] wi, int wiOffset, double[] vl, int vlOffset, int ldvl, double[] vr, int vrOffset,
            int ldvr, double[] work, int workOffset, int lwork, intW info) {
        Dgeev.dgeev(jobvl, jobvr, n, a, aOffset, lda, wr, wrOffset, wi, wiOffset, vl, vlOffset, ldvl, vr, vrOffset,
                ldvr, work, workOffset, lwork, info);
    }

    @Override
    public final void dgelqf(int m, int n, double[] a, int aOffset, int lda, double[] tau, int tauOffset, double[] work,
            int workOffset, int lwork, intW info) {
        Dgelqf.dgelqf(m, n, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork, info);
    }

    @Override
    public final void dgels(String trans, int m, int n, int nrhs, double[] a, int aOffset, int lda, double[] b,
            int bOffset, int ldb, double[] work, int workOffset, int lwork, intW info) {
        Dgels.dgels(trans, m, n, nrhs, a, aOffset, lda, b, bOffset, ldb, work, workOffset, lwork, info);
    }

    @Override
    public final void dgeqlf(int m, int n, double[] a, int aOffset, int lda, double[] tau, int tauOffset, double[] work,
            int workOffset, int lwork, intW info) {
        Dgeqlf.dgeqlf(m, n, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork, info);
    }

    @Override
    public final void dgeqp3(int m, int n, double[] a, int aOffset, int lda, int[] jpvt, int jpvtOffset, double[] tau,
            int tauOffset, double[] work, int workOffset, int lwork, intW info) {
        Dgeqp3.dgeqp3(m, n, a, aOffset, lda, jpvt, jpvtOffset, tau, tauOffset, work, workOffset, lwork, info);
    }

    @Override
    public final void dgeqrf(int m, int n, double[] a, int aOffset, int lda, double[] tau, int tauOffset, double[] work,
            int workOffset, int lwork, intW info) {
        Dgeqrf.dgeqrf(m, n, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork, info);
    }

    @Override
    public final void dgerqf(int m, int n, double[] a, int aOffset, int lda, double[] tau, int tauOffset, double[] work,
            int workOffset, int lwork, intW info) {
        Dgerqf.dgerqf(m, n, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork, info);
    }

    @Override
    public final void dgesdd(String jobz, int m, int n, double[] a, int aOffset, int lda, double[] s, int sOffset,
            double[] u, int uOffset, int ldu, double[] vt, int vtOffset, int ldvt, double[] work, int workOffset,
            int lwork, int[] iwork, int iworkOffset, intW info) {
        Dgesdd.dgesdd(jobz, m, n, a, aOffset, lda, s, sOffset, u, uOffset, ldu, vt, vtOffset, ldvt, work, workOffset,
                lwork, iwork, iworkOffset, info);
    }

    @Override
    public final void dgesv(int n, int nrhs, double[] a, int aOffset, int lda, int[] ipiv, int ipivOffset, double[] b,
            int bOffset, int ldb, intW info) {
        Dgesv.dgesv(n, nrhs, a, aOffset, lda, ipiv, ipivOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void dgetrf(int m, int n, double[] a, int aOffset, int lda, int[] ipiv, int ipivOffset, intW info) {
        Dgetrf.dgetrf(m, n, a, aOffset, lda, ipiv, ipivOffset, info);
    }

    @Override
    public final void dgetrs(String trans, int n, int nrhs, double[] a, int aOffset, int lda, int[] ipiv,
            int ipivOffset, double[] b, int bOffset, int ldb, intW info) {
        Dgetrs.dgetrs(trans, n, nrhs, a, aOffset, lda, ipiv, ipivOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void dgtsv(int n, int nrhs, double[] dl, int dlOffset, double[] d, int dOffset, double[] du,
            int duOffset, double[] b, int bOffset, int ldb, intW info) {
        Dgtsv.dgtsv(n, nrhs, dl, dlOffset, d, dOffset, du, duOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void dlaswp(int n, double[] a, int aOffset, int lda, int k1, int k2, int[] ipiv, int ipivOffset,
            int incx) {
        Dlaswp.dlaswp(n, a, aOffset, lda, k1, k2, ipiv, ipivOffset, incx);
    }

    @Override
    public final void dorglq(int m, int n, int k, double[] a, int aOffset, int lda, double[] tau, int tauOffset,
            double[] work, int workOffset, int lwork, intW info) {
        Dorglq.dorglq(m, n, k, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork, info);
    }

    @Override
    public final void dorgql(int m, int n, int k, double[] a, int aOffset, int lda, double[] tau, int tauOffset,
            double[] work, int workOffset, int lwork, intW info) {
        Dorgql.dorgql(m, n, k, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork, info);
    }

    @Override
    public final void dorgqr(int m, int n, int k, double[] a, int aOffset, int lda, double[] tau, int tauOffset,
            double[] work, int workOffset, int lwork, intW info) {
        Dorgqr.dorgqr(m, n, k, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork, info);
    }

    @Override
    public final void dorgrq(int m, int n, int k, double[] a, int aOffset, int lda, double[] tau, int tauOffset,
            double[] work, int workOffset, int lwork, intW info) {
        Dorgrq.dorgrq(m, n, k, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork, info);
    }

    @Override
    public final void dormrz(String side, String trans, int m, int n, int k, int l, double[] a, int aOffset, int lda,
            double[] tau, int tauOffset, double[] c, int cOffset, int ldc, double[] work, int workOffset, int lwork,
            intW info) {
        Dormrz.dormrz(side, trans, m, n, k, l, a, aOffset, lda, tau, tauOffset, c, cOffset, ldc, work, workOffset,
                lwork, info);
    }

    @Override
    public final void dpbcon(String uplo, int n, int kd, double[] ab, int abOffset, int ldab, double anorm,
            doubleW rcond, double[] work, int workOffset, int[] iwork, int iworkOffset, intW info) {
        Dpbcon.dpbcon(uplo, n, kd, ab, abOffset, ldab, anorm, rcond, work, workOffset, iwork, iworkOffset, info);
    }

    @Override
    public final void dpbsv(String uplo, int n, int kd, int nrhs, double[] ab, int abOffset, int ldab, double[] b,
            int bOffset, int ldb, intW info) {
        Dpbsv.dpbsv(uplo, n, kd, nrhs, ab, abOffset, ldab, b, bOffset, ldb, info);
    }

    @Override
    public final void dpbtrf(String uplo, int n, int kd, double[] ab, int abOffset, int ldab, intW info) {
        Dpbtrf.dpbtrf(uplo, n, kd, ab, abOffset, ldab, info);
    }

    @Override
    public final void dpbtrs(String uplo, int n, int kd, int nrhs, double[] ab, int abOffset, int ldab, double[] b,
            int bOffset, int ldb, intW info) {
        Dpbtrs.dpbtrs(uplo, n, kd, nrhs, ab, abOffset, ldab, b, bOffset, ldb, info);
    }

    @Override
    public final void dpocon(String uplo, int n, double[] a, int aOffset, int lda, double anorm, doubleW rcond,
            double[] work, int workOffset, int[] iwork, int iworkOffset, intW info) {
        Dpocon.dpocon(uplo, n, a, aOffset, lda, anorm, rcond, work, workOffset, iwork, iworkOffset, info);
    }

    @Override
    public final void dposv(String uplo, int n, int nrhs, double[] a, int aOffset, int lda, double[] b, int bOffset,
            int ldb, intW info) {
        Dposv.dposv(uplo, n, nrhs, a, aOffset, lda, b, bOffset, ldb, info);
    }

    @Override
    public final void dpotrf(String uplo, int n, double[] a, int aOffset, int lda, intW info) {
        Dpotrf.dpotrf(uplo, n, a, aOffset, lda, info);
    }

    @Override
    public final void dpotrs(String uplo, int n, int nrhs, double[] a, int aOffset, int lda, double[] b, int bOffset,
            int ldb, intW info) {
        Dpotrs.dpotrs(uplo, n, nrhs, a, aOffset, lda, b, bOffset, ldb, info);
    }

    @Override
    public final void dppcon(String uplo, int n, double[] ap, int apOffset, double anorm, doubleW rcond, double[] work,
            int workOffset, int[] iwork, int iworkOffset, intW info) {
        Dppcon.dppcon(uplo, n, ap, apOffset, anorm, rcond, work, workOffset, iwork, iworkOffset, info);
    }

    @Override
    public final void dppsv(String uplo, int n, int nrhs, double[] ap, int apOffset, double[] b, int bOffset, int ldb,
            intW info) {
        Dppsv.dppsv(uplo, n, nrhs, ap, apOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void dpptrf(String uplo, int n, double[] ap, int apOffset, intW info) {
        Dpptrf.dpptrf(uplo, n, ap, apOffset, info);
    }

    @Override
    public final void dpptrs(String uplo, int n, int nrhs, double[] ap, int apOffset, double[] b, int bOffset, int ldb,
            intW info) {
        Dpptrs.dpptrs(uplo, n, nrhs, ap, apOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void dptsv(int n, int nrhs, double[] d, int dOffset, double[] e, int eOffset, double[] b, int bOffset,
            int ldb, intW info) {
        Dptsv.dptsv(n, nrhs, d, dOffset, e, eOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void dsbevd(String jobz, String uplo, int n, int kd, double[] ab, int abOffset, int ldab, double[] w,
            int wOffset, double[] z, int zOffset, int ldz, double[] work, int workOffset, int lwork, int[] iwork,
            int iworkOffset, int liwork, intW info) {
        Dsbevd.dsbevd(jobz, uplo, n, kd, ab, abOffset, ldab, w, wOffset, z, zOffset, ldz, work, workOffset, lwork,
                iwork, iworkOffset, liwork, info);
    }

    @Override
    public final void dspevd(String jobz, String uplo, int n, double[] ap, int apOffset, double[] w, int wOffset,
            double[] z, int zOffset, int ldz, double[] work, int workOffset, int lwork, int[] iwork, int iworkOffset,
            int liwork, intW info) {
        Dspevd.dspevd(jobz, uplo, n, ap, apOffset, w, wOffset, z, zOffset, ldz, work, workOffset, lwork, iwork,
                iworkOffset, liwork, info);
    }

    @Override
    public final void dspsv(String uplo, int n, int nrhs, double[] ap, int apOffset, int[] ipiv, int ipivOffset,
            double[] b, int bOffset, int ldb, intW info) {
        Dspsv.dspsv(uplo, n, nrhs, ap, apOffset, ipiv, ipivOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void dstevr(String jobz, String range, int n, double[] d, int dOffset, double[] e, int eOffset,
            double vl, double vu, int il, int iu, double abstol, intW m, double[] w, int wOffset, double[] z,
            int zOffset, int ldz, int[] isuppz, int isuppzOffset, double[] work, int workOffset, int lwork, int[] iwork,
            int iworkOffset, int liwork, intW info) {
        Dstevr.dstevr(jobz, range, n, d, dOffset, e, eOffset, vl, vu, il, iu, abstol, m, w, wOffset, z, zOffset, ldz,
                isuppz, isuppzOffset, work, workOffset, lwork, iwork, iworkOffset, liwork, info);
    }

    @Override
    public final void dsyevr(String jobz, String range, String uplo, int n, double[] a, int aOffset, int lda, double vl,
            double vu, int il, int iu, double abstol, intW m, double[] w, int wOffset, double[] z, int zOffset, int ldz,
            int[] isuppz, int isuppzOffset, double[] work, int workOffset, int lwork, int[] iwork, int iworkOffset,
            int liwork, intW info) {
        Dsyevr.dsyevr(jobz, range, uplo, n, a, aOffset, lda, vl, vu, il, iu, abstol, m, w, wOffset, z, zOffset, ldz,
                isuppz, isuppzOffset, work, workOffset, lwork, iwork, iworkOffset, liwork, info);
    }

    @Override
    public final void dsygvd(int itype, String jobz, String uplo, int n, double[] a, int aOffset, int lda, double[] b,
            int bOffset, int ldb, double[] w, int wOffset, double[] work, int workOffset, int lwork, int[] iwork,
            int iworkOffset, int liwork, intW info) {
        Dsygvd.dsygvd(itype, jobz, uplo, n, a, aOffset, lda, b, bOffset, ldb, w, wOffset, work, workOffset, lwork,
                iwork, iworkOffset, liwork, info);
    }

    @Override
    public final void dsysv(String uplo, int n, int nrhs, double[] a, int aOffset, int lda, int[] ipiv, int ipivOffset,
            double[] b, int bOffset, int ldb, double[] work, int workOffset, int lwork, intW info) {
        Dsysv.dsysv(uplo, n, nrhs, a, aOffset, lda, ipiv, ipivOffset, b, bOffset, ldb, work, workOffset, lwork, info);
    }

    @Override
    public final void dtbtrs(String uplo, String trans, String diag, int n, int kd, int nrhs, double[] ab, int abOffset,
            int ldab, double[] b, int bOffset, int ldb, intW info) {
        Dtbtrs.dtbtrs(uplo, trans, diag, n, kd, nrhs, ab, abOffset, ldab, b, bOffset, ldb, info);
    }

    @Override
    public final void dtptrs(String uplo, String trans, String diag, int n, int nrhs, double[] ap, int apOffset,
            double[] b, int bOffset, int ldb, intW info) {
        Dtptrs.dtptrs(uplo, trans, diag, n, nrhs, ap, apOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void dtrtrs(String uplo, String trans, String diag, int n, int nrhs, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, intW info) {
        Dtrtrs.dtrtrs(uplo, trans, diag, n, nrhs, a, aOffset, lda, b, bOffset, ldb, info);
    }

    // miscellaneous float routines

    @Override
    public final void sgeev(String jobvl, String jobvr, int n, float[] a, int aOffset, int lda, float[] wr,
            int wrOffset, float[] wi, int wiOffset, float[] vl, int vlOffset, int ldvl, float[] vr, int vrOffset,
            int ldvr, float[] work, int workOffset, int lwork, intW info) {
        Sgeev.sgeev(jobvl, jobvr, n, a, aOffset, lda, wr, wrOffset, wi, wiOffset, vl, vlOffset, ldvl, vr, vrOffset,
                ldvr, work, workOffset, lwork, info);
    }

    @Override
    public final void sgels(String trans, int m, int n, int nrhs, float[] a, int aOffset, int lda, float[] b,
            int bOffset, int ldb, float[] work, int workOffset, int lwork, intW info) {
        Sgels.sgels(trans, m, n, nrhs, a, aOffset, lda, b, bOffset, ldb, work, workOffset, lwork, info);
    }

    @Override
    public final void sgesdd(String jobz, int m, int n, float[] a, int aOffset, int lda, float[] s, int sOffset,
            float[] u, int uOffset, int ldu, float[] vt, int vtOffset, int ldvt, float[] work, int workOffset,
            int lwork, int[] iwork, int iworkOffset, intW info) {
        Sgesdd.sgesdd(jobz, m, n, a, aOffset, lda, s, sOffset, u, uOffset, ldu, vt, vtOffset, ldvt, work, workOffset,
                lwork, iwork, iworkOffset, info);
    }

    @Override
    public final void sgesv(int n, int nrhs, float[] a, int aOffset, int lda, int[] ipiv, int ipivOffset, float[] b,
            int bOffset, int ldb, intW info) {
        Sgesv.sgesv(n, nrhs, a, aOffset, lda, ipiv, ipivOffset, b, bOffset, ldb, info);
    }

    @Override
    public final void sgeqrf(int m, int n, float[] a, int aOffset, int lda, float[] tau, int tauOffset, float[] work,
            int workOffset, int lwork, intW info) {
        Sgeqrf.sgeqrf(m, n, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork, info);
    }

    @Override
    public final void sorgqr(int m, int n, int k, float[] a, int aOffset, int lda, float[] tau, int tauOffset,
            float[] work, int workOffset, int lwork, intW info) {
        Sorgqr.sorgqr(m, n, k, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork, info);
    }

    // miscellaneous complex routines

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int cgeev(String jobvl, String jobvr, int n, float[] a, int lda, float[] w, float[] vl, int ldvl,
            float[] vr, int ldvr) {
        throw new UnsupportedOperationException("cgeev only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int zgeev(String jobvl, String jobvr, int n, double[] a, int lda, double[] w, double[] vl, int ldvl,
            double[] vr, int ldvr) {
        throw new UnsupportedOperationException("zgeev only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int cgesdd(String jobz, int m, int n, float[] a, int lda, float[] s, float[] u, int ldu, float[] vt,
            int ldvt) {
        throw new UnsupportedOperationException("cgesdd only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int zgesdd(String jobz, int m, int n, double[] a, int lda, double[] s, double[] u, int ldu,
            double[] vt, int ldvt) {
        throw new UnsupportedOperationException("zgesdd only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int cgels(String trans, int m, int n, int nrhs, float[] a, int lda, float[] b, int ldb) {
        throw new UnsupportedOperationException("cgels only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int zgels(String trans, int m, int n, int nrhs, double[] a, int lda, double[] b, int ldb) {
        throw new UnsupportedOperationException("zgels only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int cgesv(int n, int nrhs, float[] a, int lda, int[] ipiv, float[] b, int ldb) {
        throw new UnsupportedOperationException("cgesv only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int zgesv(int n, int nrhs, double[] a, int lda, int[] ipiv, double[] b, int ldb) {
        throw new UnsupportedOperationException("zgesv only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int cgeqrf(int m, int n, float[] a, int lda, float[] tau) {
        throw new UnsupportedOperationException("cgeqrf only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int zgeqrf(int m, int n, double[] a, int lda, double[] tau) {
        throw new UnsupportedOperationException("zgeqrf only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int cungqr(int m, int n, int k, float[] a, int lda, float[] tau) {
        throw new UnsupportedOperationException("cungqr only supported for LapackN");
    }

    /**
     * Not supported in the {@code Netlib} implementation. Use {@link LapackN}
     * instead.
     * 
     * @throws UnsupportedOperationException
     *             always
     */
    @Override
    public final int zungqr(int m, int n, int k, double[] a, int lda, double[] tau) {
        throw new UnsupportedOperationException("zungqr only supported for LapackN");
    }

    protected LapackJ() {
    }
}
