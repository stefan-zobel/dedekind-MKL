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
package net.dedekind.lapack;

import java.util.Objects;

import org.netlib.util.doubleW;
import org.netlib.util.intW;

import net.dedekind.Order;

/**
 * {@code Intel MKL} (Math Kernel Library) native implementation.
 */
public class LapackN extends Lapack {

    // TODO
    private static final boolean USE_CRITICAL = true;

    @Override
    public final void dgbcon(String norm, int n, int kl, int ku, double[] ab, int abOffset, int ldab, int[] ipiv,
            int ipivOffset, double anorm, doubleW rcondDW, double[] work, int workOffset, int[] iwork, int iworkOffset,
            intW info) {
        Objects.requireNonNull(ab, "ab");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(rcondDW, "rcond");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = dgbcon_n(Order.COL.code(), norm(norm), n, kl, ku, ab, abOffset, ldab, ipiv, ipivOffset, anorm,
                rcondDW, work, workOffset, iwork, iworkOffset, USE_CRITICAL);
    }

    @Override
    public final void dgbsv(int n, int kl, int ku, int nrhs, double[] ab, int abOffset, int ldab, int[] ipiv,
            int ipivOffset, double[] b, int bOffset, int ldb, intW info) {
        Objects.requireNonNull(ab, "ab");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dgbsv_n(Order.COL.code(), n, kl, ku, nrhs, ab, abOffset, ldab, ipiv, ipivOffset, b, bOffset, ldb,
                USE_CRITICAL);
    }

    @Override
    public final void dgbtrf(int m, int n, int kl, int ku, double[] ab, int abOffset, int ldab, int[] ipiv,
            int ipivOffset, intW info) {
        Objects.requireNonNull(ab, "ab");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(info, "info");
        info.val = dgbtrf_n(Order.COL.code(), m, n, kl, ku, ab, abOffset, ldab, ipiv, ipivOffset, USE_CRITICAL);
    }

    @Override
    public final void dgbtrs(String trans, int n, int kl, int ku, int nrhs, double[] ab, int abOffset, int ldab,
            int[] ipiv, int ipivOffset, double[] b, int bOffset, int ldb, intW info) {
        Objects.requireNonNull(ab, "ab");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dgbtrs_n(Order.COL.code(), trans(trans), n, kl, ku, nrhs, ab, abOffset, ldab, ipiv, ipivOffset, b,
                bOffset, ldb, USE_CRITICAL);
    }

    @Override
    public final void dgecon(String norm, int n, double[] a, int aOffset, int lda, double anorm, doubleW rcondDW,
            double[] work, int workOffset, int[] iwork, int iworkOffset, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(rcondDW, "rcond");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = dgecon_n(Order.COL.code(), norm(norm), n, a, aOffset, lda, anorm, rcondDW, work, workOffset, iwork,
                iworkOffset, USE_CRITICAL);
    }

    @Override
    public final void dgeev(String jobvl, String jobvr, int n, double[] a, int aOffset, int lda, double[] wr,
            int wrOffset, double[] wi, int wiOffset, double[] vl, int vlOffset, int ldvl, double[] vr, int vrOffset,
            int ldvr, double[] work, int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(wr, "wr");
        Objects.requireNonNull(wi, "wi");
        Objects.requireNonNull(vl, "vl");
        Objects.requireNonNull(vr, "vr");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dgeev_n(Order.COL.code(), eigJob(jobvl), eigJob(jobvr), n, a, aOffset, lda, wr, wrOffset, wi,
                wiOffset, vl, vlOffset, ldvl, vr, vrOffset, ldvr, work, workOffset, lwork, USE_CRITICAL);
    }

    @Override
    public final void dgelqf(int m, int n, double[] a, int aOffset, int lda, double[] tau, int tauOffset, double[] work,
            int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(tau, "tau");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dgelqf_n(Order.COL.code(), m, n, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork,
                USE_CRITICAL);
    }

    @Override
    public final void dgels(String trans, int m, int n, int nrhs, double[] a, int aOffset, int lda, double[] b,
            int bOffset, int ldb, double[] work, int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dgels_n(Order.COL.code(), trans(trans), m, n, nrhs, a, aOffset, lda, b, bOffset, ldb, work,
                workOffset, lwork, USE_CRITICAL);
    }

    @Override
    public final void dgeqlf(int m, int n, double[] a, int aOffset, int lda, double[] tau, int tauOffset, double[] work,
            int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(tau, "tau");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dgeqlf_n(Order.COL.code(), m, n, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork,
                USE_CRITICAL);
    }

    @Override
    public final void dgeqp3(int m, int n, double[] a, int aOffset, int lda, int[] jpvt, int jpvtOffset, double[] tau,
            int tauOffset, double[] work, int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(jpvt, "jpvt");
        Objects.requireNonNull(tau, "tau");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dgeqp3_n(Order.COL.code(), m, n, a, aOffset, lda, jpvt, jpvtOffset, tau, tauOffset, work, workOffset,
                lwork, USE_CRITICAL);
    }

    @Override
    public final void dgeqrf(int m, int n, double[] a, int aOffset, int lda, double[] tau, int tauOffset, double[] work,
            int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(tau, "tau");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dgeqrf_n(Order.COL.code(), m, n, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork,
                USE_CRITICAL);
    }

    @Override
    public final void dgerqf(int m, int n, double[] a, int aOffset, int lda, double[] tau, int tauOffset, double[] work,
            int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(tau, "tau");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dgerqf_n(Order.COL.code(), m, n, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork,
                USE_CRITICAL);
    }

    @Override
    public final void dgesdd(String jobz, int m, int n, double[] a, int aOffset, int lda, double[] s, int sOffset,
            double[] u, int uOffset, int ldu, double[] vt, int vtOffset, int ldvt, double[] work, int workOffset,
            int lwork, int[] iwork, int iworkOffset, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(s, "s");
        Objects.requireNonNull(u, "u");
        Objects.requireNonNull(vt, "vt");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = dgesdd_n(Order.COL.code(), svdJob(jobz), m, n, a, aOffset, lda, s, sOffset, u, uOffset, ldu, vt,
                vtOffset, ldvt, work, workOffset, lwork, iwork, iworkOffset, USE_CRITICAL);
    }

    @Override
    public final void dgesv(int n, int nrhs, double[] a, int aOffset, int lda, int[] ipiv, int ipivOffset, double[] b,
            int bOffset, int ldb, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dgesv_n(Order.COL.code(), n, nrhs, a, aOffset, lda, ipiv, ipivOffset, b, bOffset, ldb, USE_CRITICAL);
    }

    @Override
    public final void dgetrf(int m, int n, double[] a, int aOffset, int lda, int[] ipiv, int ipivOffset, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(info, "info");
        info.val = dgetrf_n(Order.COL.code(), m, n, a, aOffset, lda, ipiv, ipivOffset, USE_CRITICAL);
    }

    @Override
    public final void dgetrs(String trans, int n, int nrhs, double[] a, int aOffset, int lda, int[] ipiv,
            int ipivOffset, double[] b, int bOffset, int ldb, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dgetrs_n(Order.COL.code(), trans(trans), n, nrhs, a, aOffset, lda, ipiv, ipivOffset, b, bOffset, ldb,
                USE_CRITICAL);
    }

    @Override
    public final void dgtsv(int n, int nrhs, double[] dl, int dlOffset, double[] d, int dOffset, double[] du,
            int duOffset, double[] b, int bOffset, int ldb, intW info) {
        Objects.requireNonNull(dl, "dl");
        Objects.requireNonNull(d, "d");
        Objects.requireNonNull(du, "du");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dgtsv_n(Order.COL.code(), n, nrhs, dl, dlOffset, d, dOffset, du, duOffset, b, bOffset, ldb,
                USE_CRITICAL);
    }

    @Override
    public final void dlaswp(int n, double[] a, int aOffset, int lda, int k1, int k2, int[] ipiv, int ipivOffset,
            int incx) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(ipiv, "ipiv");
        dlaswp_n(Order.COL.code(), n, a, aOffset, lda, k1, k2, ipiv, ipivOffset, incx, USE_CRITICAL);
    }

    @Override
    public final void dorglq(int m, int n, int k, double[] a, int aOffset, int lda, double[] tau, int tauOffset,
            double[] work, int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(tau, "tau");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dorglq_n(Order.COL.code(), m, n, k, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork,
                USE_CRITICAL);
    }

    @Override
    public final void dorgql(int m, int n, int k, double[] a, int aOffset, int lda, double[] tau, int tauOffset,
            double[] work, int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(tau, "tau");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dorgql_n(Order.COL.code(), m, n, k, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork,
                USE_CRITICAL);
    }

    @Override
    public final void dorgqr(int m, int n, int k, double[] a, int aOffset, int lda, double[] tau, int tauOffset,
            double[] work, int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(tau, "tau");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dorgqr_n(Order.COL.code(), m, n, k, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork,
                USE_CRITICAL);
    }

    @Override
    public final void dorgrq(int m, int n, int k, double[] a, int aOffset, int lda, double[] tau, int tauOffset,
            double[] work, int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(tau, "tau");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dorgrq_n(Order.COL.code(), m, n, k, a, aOffset, lda, tau, tauOffset, work, workOffset, lwork,
                USE_CRITICAL);
    }

    @Override
    public final void dormrz(String side, String trans, int m, int n, int k, int l, double[] a, int aOffset, int lda,
            double[] tau, int tauOffset, double[] c, int cOffset, int ldc, double[] work, int workOffset, int lwork,
            intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(tau, "tau");
        Objects.requireNonNull(c, "c");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dormrz_n(Order.COL.code(), side(side), trans(trans), m, n, k, l, a, aOffset, lda, tau, tauOffset, c,
                cOffset, ldc, work, workOffset, lwork, USE_CRITICAL);
    }

    @Override
    public final void dpbcon(String uplo, int n, int kd, double[] ab, int abOffset, int ldab, double anorm,
            doubleW rcondDW, double[] work, int workOffset, int[] iwork, int iworkOffset, intW info) {
        Objects.requireNonNull(ab, "ab");
        Objects.requireNonNull(rcondDW, "rcond");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dpbcon_n(Order.COL.code(), uplo(uplo), n, kd, ab, abOffset, ldab, anorm, rcondDW, work, workOffset,
                iwork, iworkOffset, USE_CRITICAL);
    }

    @Override
    public final void dpbsv(String uplo, int n, int kd, int nrhs, double[] ab, int abOffset, int ldab, double[] b,
            int bOffset, int ldb, intW info) {
        Objects.requireNonNull(ab, "ab");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dpbsv_n(Order.COL.code(), uplo(uplo), n, kd, nrhs, ab, abOffset, ldab, b, bOffset, ldb,
                USE_CRITICAL);
    }

    @Override
    public final void dpbtrf(String uplo, int n, int kd, double[] ab, int abOffset, int ldab, intW info) {
        Objects.requireNonNull(ab, "ab");
        Objects.requireNonNull(info, "info");
        info.val = dpbtrf_n(Order.COL.code(), uplo(uplo), n, kd, ab, abOffset, ldab, USE_CRITICAL);
    }

    @Override
    public final void dpbtrs(String uplo, int n, int kd, int nrhs, double[] ab, int abOffset, int ldab, double[] b,
            int bOffset, int ldb, intW info) {
        Objects.requireNonNull(ab, "ab");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dpbtrs_n(Order.COL.code(), uplo(uplo), n, kd, nrhs, ab, abOffset, ldab, b, bOffset, ldb,
                USE_CRITICAL);
    }

    @Override
    public final void dpocon(String uplo, int n, double[] a, int aOffset, int lda, double anorm, doubleW rcondDW,
            double[] work, int workOffset, int[] iwork, int iworkOffset, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(rcondDW, "rcond");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = dpocon_n(Order.COL.code(), uplo(uplo), n, a, aOffset, lda, anorm, rcondDW, work, workOffset, iwork,
                iworkOffset, USE_CRITICAL);
    }

    @Override
    public final void dposv(String uplo, int n, int nrhs, double[] a, int aOffset, int lda, double[] b, int bOffset,
            int ldb, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dposv_n(Order.COL.code(), uplo(uplo), n, nrhs, a, aOffset, lda, b, bOffset, ldb, USE_CRITICAL);
    }

    @Override
    public final void dpotrf(String uplo, int n, double[] a, int aOffset, int lda, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(info, "info");
        info.val = dpotrf_n(Order.COL.code(), uplo(uplo), n, a, aOffset, lda, USE_CRITICAL);
    }

    @Override
    public final void dpotrs(String uplo, int n, int nrhs, double[] a, int aOffset, int lda, double[] b, int bOffset,
            int ldb, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dpotrs_n(Order.COL.code(), uplo(uplo), n, nrhs, a, aOffset, lda, b, bOffset, ldb, USE_CRITICAL);
    }

    @Override
    public final void dppcon(String uplo, int n, double[] ap, int apOffset, double anorm, doubleW rcondDW,
            double[] work, int workOffset, int[] iwork, int iworkOffset, intW info) {
        Objects.requireNonNull(ap, "ap");
        Objects.requireNonNull(rcondDW, "rcond");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = dppcon_n(Order.COL.code(), uplo(uplo), n, ap, apOffset, anorm, rcondDW, work, workOffset, iwork,
                iworkOffset, USE_CRITICAL);
    }

    @Override
    public final void dppsv(String uplo, int n, int nrhs, double[] ap, int apOffset, double[] b, int bOffset, int ldb,
            intW info) {
        Objects.requireNonNull(ap, "ap");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dppsv_n(Order.COL.code(), uplo(uplo), n, nrhs, ap, apOffset, b, bOffset, ldb, USE_CRITICAL);
    }

    @Override
    public final void dpptrf(String uplo, int n, double[] ap, int apOffset, intW info) {
        Objects.requireNonNull(ap, "ap");
        Objects.requireNonNull(info, "info");
        info.val = dpptrf_n(Order.COL.code(), uplo(uplo), n, ap, apOffset, USE_CRITICAL);
    }

    @Override
    public final void dpptrs(String uplo, int n, int nrhs, double[] ap, int apOffset, double[] b, int bOffset, int ldb,
            intW info) {
        Objects.requireNonNull(ap, "ap");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dpptrs_n(Order.COL.code(), uplo(uplo), n, nrhs, ap, apOffset, b, bOffset, ldb, USE_CRITICAL);
    }

    @Override
    public final void dptsv(int n, int nrhs, double[] d, int dOffset, double[] e, int eOffset, double[] b, int bOffset,
            int ldb, intW info) {
        Objects.requireNonNull(d, "d");
        Objects.requireNonNull(e, "e");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dptsv_n(Order.COL.code(), n, nrhs, d, dOffset, e, eOffset, b, bOffset, ldb, USE_CRITICAL);
    }

    @Override
    public final void dsbevd(String jobz, String uplo, int n, int kd, double[] ab, int abOffset, int ldab, double[] w,
            int wOffset, double[] z, int zOffset, int ldz, double[] work, int workOffset, int lwork, int[] iwork,
            int iworkOffset, int liwork, intW info) {
        Objects.requireNonNull(ab, "ab");
        Objects.requireNonNull(w, "w");
        Objects.requireNonNull(z, "z");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = dsbevd_n(Order.COL.code(), eigJob(jobz), uplo(uplo), n, kd, ab, abOffset, ldab, w, wOffset, z,
                zOffset, ldz, work, workOffset, lwork, iwork, iworkOffset, liwork, USE_CRITICAL);
    }

    @Override
    public final void dspevd(String jobz, String uplo, int n, double[] ap, int apOffset, double[] w, int wOffset,
            double[] z, int zOffset, int ldz, double[] work, int workOffset, int lwork, int[] iwork, int iworkOffset,
            int liwork, intW info) {
        Objects.requireNonNull(ap, "ap");
        Objects.requireNonNull(w, "w");
        Objects.requireNonNull(z, "z");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = dspevd_n(Order.COL.code(), eigJob(jobz), uplo(uplo), n, ap, apOffset, w, wOffset, z, zOffset, ldz,
                work, workOffset, lwork, iwork, iworkOffset, liwork, USE_CRITICAL);
    }

    @Override
    public final void dspsv(String uplo, int n, int nrhs, double[] ap, int apOffset, int[] ipiv, int ipivOffset,
            double[] b, int bOffset, int ldb, intW info) {
        Objects.requireNonNull(ap, "ap");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dspsv_n(Order.COL.code(), uplo(uplo), n, nrhs, ap, apOffset, ipiv, ipivOffset, b, bOffset, ldb,
                USE_CRITICAL);
    }

    @Override
    public final void dstevr(String jobz, String range, int n, double[] d, int dOffset, double[] e, int eOffset,
            double vl, double vu, int il, int iu, double abstol, intW mIW, double[] w, int wOffset, double[] z,
            int zOffset, int ldz, int[] isuppz, int isuppzOffset, double[] work, int workOffset, int lwork, int[] iwork,
            int iworkOffset, int liwork, intW info) {
        Objects.requireNonNull(d, "d");
        Objects.requireNonNull(e, "e");
        Objects.requireNonNull(mIW, "m");
        Objects.requireNonNull(w, "w");
        Objects.requireNonNull(z, "z");
        Objects.requireNonNull(isuppz, "isuppz");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = dstevr_n(Order.COL.code(), eigJob(jobz), range(range), n, d, dOffset, e, eOffset, vl, vu, il, iu,
                abstol, mIW, w, wOffset, z, zOffset, ldz, isuppz, isuppzOffset, work, workOffset, lwork, iwork,
                iworkOffset, liwork, USE_CRITICAL);
    }

    @Override
    public final void dsyevr(String jobz, String range, String uplo, int n, double[] a, int aOffset, int lda, double vl,
            double vu, int il, int iu, double abstol, intW mIW, double[] w, int wOffset, double[] z, int zOffset,
            int ldz, int[] isuppz, int isuppzOffset, double[] work, int workOffset, int lwork, int[] iwork,
            int iworkOffset, int liwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(mIW, "m");
        Objects.requireNonNull(w, "w");
        Objects.requireNonNull(z, "z");
        Objects.requireNonNull(isuppz, "isuppz");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = dsyevr_n(Order.COL.code(), eigJob(jobz), range(range), uplo(uplo), n, a, aOffset, lda, vl, vu, il,
                iu, abstol, mIW, w, wOffset, z, zOffset, ldz, isuppz, isuppzOffset, work, workOffset, lwork, iwork,
                iworkOffset, liwork, USE_CRITICAL);
    }

    @Override
    public final void dsygvd(int itype, String jobz, String uplo, int n, double[] a, int aOffset, int lda, double[] b,
            int bOffset, int ldb, double[] w, int wOffset, double[] work, int workOffset, int lwork, int[] iwork,
            int iworkOffset, int liwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(w, "w");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = dsygvd_n(Order.COL.code(), itype, eigJob(jobz), uplo(uplo), n, a, aOffset, lda, b, bOffset, ldb, w,
                wOffset, work, workOffset, lwork, iwork, iworkOffset, liwork, USE_CRITICAL);
    }

    @Override
    public final void dsysv(String uplo, int n, int nrhs, double[] a, int aOffset, int lda, int[] ipiv, int ipivOffset,
            double[] b, int bOffset, int ldb, double[] work, int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = dsysv_n(Order.COL.code(), uplo(uplo), n, nrhs, a, aOffset, lda, ipiv, ipivOffset, b, bOffset, ldb,
                work, workOffset, lwork, USE_CRITICAL);
    }

    @Override
    public final void dtbtrs(String uplo, String trans, String diag, int n, int kd, int nrhs, double[] ab, int abOffset,
            int ldab, double[] b, int bOffset, int ldb, intW info) {
        Objects.requireNonNull(ab, "ab");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dtbtrs_n(Order.COL.code(), uplo(uplo), trans(trans), diag(diag), n, kd, nrhs, ab, abOffset, ldab, b,
                bOffset, ldb, USE_CRITICAL);
    }

    @Override
    public final void dtptrs(String uplo, String trans, String diag, int n, int nrhs, double[] ap, int apOffset,
            double[] b, int bOffset, int ldb, intW info) {
        Objects.requireNonNull(ap, "ap");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dtptrs_n(Order.COL.code(), uplo(uplo), trans(trans), diag(diag), n, nrhs, ap, apOffset, b, bOffset,
                ldb, USE_CRITICAL);
    }

    @Override
    public final void dtrtrs(String uplo, String trans, String diag, int n, int nrhs, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = dtrtrs_n(Order.COL.code(), uplo(uplo), trans(trans), diag(diag), n, nrhs, a, aOffset, lda, b,
                bOffset, ldb, USE_CRITICAL);
    }

    // miscellaneous float routines

    @Override
    public final void sgeev(String jobvl, String jobvr, int n, float[] a, int aOffset, int lda, float[] wr,
            int wrOffset, float[] wi, int wiOffset, float[] vl, int vlOffset, int ldvl, float[] vr, int vrOffset,
            int ldvr, float[] work, int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(wr, "wr");
        Objects.requireNonNull(wi, "wi");
        Objects.requireNonNull(vl, "vl");
        Objects.requireNonNull(vr, "vr");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = sgeev_n(Order.COL.code(), eigJob(jobvl), eigJob(jobvr), n, a, aOffset, lda, wr, wrOffset, wi,
                wiOffset, vl, vlOffset, ldvl, vr, vrOffset, ldvr, work, workOffset, lwork, USE_CRITICAL);
    }

    @Override
    public final void sgels(String trans, int m, int n, int nrhs, float[] a, int aOffset, int lda, float[] b,
            int bOffset, int ldb, float[] work, int workOffset, int lwork, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(info, "info");
        info.val = sgels_n(Order.COL.code(), trans(trans), m, n, nrhs, a, aOffset, lda, b, bOffset, ldb, work,
                workOffset, lwork, USE_CRITICAL);
    }

    @Override
    public final void sgesdd(String jobz, int m, int n, float[] a, int aOffset, int lda, float[] s, int sOffset,
            float[] u, int uOffset, int ldu, float[] vt, int vtOffset, int ldvt, float[] work, int workOffset,
            int lwork, int[] iwork, int iworkOffset, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(s, "s");
        Objects.requireNonNull(u, "u");
        Objects.requireNonNull(vt, "vt");
        Objects.requireNonNull(work, "work");
        Objects.requireNonNull(iwork, "iwork");
        Objects.requireNonNull(info, "info");
        info.val = sgesdd_n(Order.COL.code(), svdJob(jobz), m, n, a, aOffset, lda, s, sOffset, u, uOffset, ldu, vt,
                vtOffset, ldvt, work, workOffset, lwork, iwork, iworkOffset, USE_CRITICAL);
    }

    @Override
    public final void sgesv(int n, int nrhs, float[] a, int aOffset, int lda, int[] ipiv, int ipivOffset, float[] b,
            int bOffset, int ldb, intW info) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(ipiv, "ipiv");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(info, "info");
        info.val = sgesv_n(Order.COL.code(), n, nrhs, a, aOffset, lda, ipiv, ipivOffset, b, bOffset, ldb, USE_CRITICAL);
    }

    // native methods

    private static native int dgbcon_n(int order, byte norm, int n, int kl, int ku, double[] ab, int abOffset, int ldab,
            int[] ipiv, int ipivOffset, double anorm, doubleW rcondDW, double[] work, int workOffset, int[] iwork,
            int iworkOffset, boolean useCriticalRegion);

    private static native int dgbsv_n(int order, int n, int kl, int ku, int nrhs, double[] ab, int abOffset, int ldab,
            int[] ipiv, int ipivOffset, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dgbtrf_n(int order, int m, int n, int kl, int ku, double[] ab, int abOffset, int ldab,
            int[] ipiv, int ipivOffset, boolean useCriticalRegion);

    private static native int dgbtrs_n(int order, byte trans, int n, int kl, int ku, int nrhs, double[] ab,
            int abOffset, int ldab, int[] ipiv, int ipivOffset, double[] b, int bOffset, int ldb,
            boolean useCriticalRegion);

    private static native int dgecon_n(int order, byte norm, int n, double[] a, int aOffset, int lda, double anorm,
            doubleW rcondDW, double[] work, int workOffset, int[] iwork, int iworkOffset, boolean useCriticalRegion);

    private static native int dgeev_n(int order, byte jobvl, byte jobvr, int n, double[] a, int aOffset, int lda,
            double[] wr, int wrOffset, double[] wi, int wiOffset, double[] vl, int vlOffset, int ldvl, double[] vr,
            int vrOffset, int ldvr, double[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dgelqf_n(int order, int m, int n, double[] a, int aOffset, int lda, double[] tau,
            int tauOffset, double[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dgels_n(int order, byte trans, int m, int n, int nrhs, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, double[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dgeqlf_n(int order, int m, int n, double[] a, int aOffset, int lda, double[] tau,
            int tauOffset, double[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dgeqp3_n(int order, int m, int n, double[] a, int aOffset, int lda, int[] jpvt,
            int jpvtOffset, double[] tau, int tauOffset, double[] work, int workOffset, int lwork,
            boolean useCriticalRegion);

    private static native int dgeqrf_n(int order, int m, int n, double[] a, int aOffset, int lda, double[] tau,
            int tauOffset, double[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dgerqf_n(int order, int m, int n, double[] a, int aOffset, int lda, double[] tau,
            int tauOffset, double[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dgesdd_n(int order, byte jobz, int m, int n, double[] a, int aOffset, int lda, double[] s,
            int sOffset, double[] u, int uOffset, int ldu, double[] vt, int vtOffset, int ldvt, double[] work,
            int workOffset, int lwork, int[] iwork, int iworkOffset, boolean useCriticalRegion);

    private static native int dgesv_n(int order, int n, int nrhs, double[] a, int aOffset, int lda, int[] ipiv,
            int ipivOffset, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dgetrf_n(int order, int m, int n, double[] a, int aOffset, int lda, int[] ipiv,
            int ipivOffset, boolean useCriticalRegion);

    private static native int dgetrs_n(int order, byte trans, int n, int nrhs, double[] a, int aOffset, int lda,
            int[] ipiv, int ipivOffset, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dgtsv_n(int order, int n, int nrhs, double[] dl, int dlOffset, double[] d, int dOffset,
            double[] du, int duOffset, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native void dlaswp_n(int order, int n, double[] a, int aOffset, int lda, int k1, int k2, int[] ipiv,
            int ipivOffset, int incx, boolean useCriticalRegion);

    private static native int dorglq_n(int order, int m, int n, int k, double[] a, int aOffset, int lda, double[] tau,
            int tauOffset, double[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dorgql_n(int order, int m, int n, int k, double[] a, int aOffset, int lda, double[] tau,
            int tauOffset, double[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dorgqr_n(int order, int m, int n, int k, double[] a, int aOffset, int lda, double[] tau,
            int tauOffset, double[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dorgrq_n(int order, int m, int n, int k, double[] a, int aOffset, int lda, double[] tau,
            int tauOffset, double[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dormrz_n(int order, byte side, byte trans, int m, int n, int k, int l, double[] a,
            int aOffset, int lda, double[] tau, int tauOffset, double[] c, int cOffset, int ldc, double[] work,
            int workOffset, int lwork, boolean useCriticalRegion);

    private static native int dpbcon_n(int order, byte uplo, int n, int kd, double[] ab, int abOffset, int ldab,
            double anorm, doubleW rcondDW, double[] work, int workOffset, int[] iwork, int iworkOffset,
            boolean useCriticalRegion);

    private static native int dpbsv_n(int order, byte uplo, int n, int kd, int nrhs, double[] ab, int abOffset,
            int ldab, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dpbtrf_n(int order, byte uplo, int n, int kd, double[] ab, int abOffset, int ldab,
            boolean useCriticalRegion);

    private static native int dpbtrs_n(int order, byte uplo, int n, int kd, int nrhs, double[] ab, int abOffset,
            int ldab, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dpocon_n(int order, byte uplo, int n, double[] a, int aOffset, int lda, double anorm,
            doubleW rcondDW, double[] work, int workOffset, int[] iwork, int iworkOffset, boolean useCriticalRegion);

    private static native int dposv_n(int order, byte uplo, int n, int nrhs, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dpotrf_n(int order, byte uplo, int n, double[] a, int aOffset, int lda,
            boolean useCriticalRegion);

    private static native int dpotrs_n(int order, byte uplo, int n, int nrhs, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dppcon_n(int order, byte uplo, int n, double[] ap, int apOffset, double anorm,
            doubleW rcondDW, double[] work, int workOffset, int[] iwork, int iworkOffset, boolean useCriticalRegion);

    private static native int dppsv_n(int order, byte uplo, int n, int nrhs, double[] ap, int apOffset, double[] b,
            int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dpptrf_n(int order, byte uplo, int n, double[] ap, int apOffset,
            boolean useCriticalRegion);

    private static native int dpptrs_n(int order, byte uplo, int n, int nrhs, double[] ap, int apOffset, double[] b,
            int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dptsv_n(int order, int n, int nrhs, double[] d, int dOffset, double[] e, int eOffset,
            double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dsbevd_n(int order, byte jobz, byte uplo, int n, int kd, double[] ab, int abOffset,
            int ldab, double[] w, int wOffset, double[] z, int zOffset, int ldz, double[] work, int workOffset,
            int lwork, int[] iwork, int iworkOffset, int liwork, boolean useCriticalRegion);

    private static native int dspevd_n(int order, byte jobz, byte uplo, int n, double[] ap, int apOffset, double[] w,
            int wOffset, double[] z, int zOffset, int ldz, double[] work, int workOffset, int lwork, int[] iwork,
            int iworkOffset, int liwork, boolean useCriticalRegion);

    private static native int dspsv_n(int order, byte uplo, int n, int nrhs, double[] ap, int apOffset, int[] ipiv,
            int ipivOffset, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dstevr_n(int order, byte jobz, byte range, int n, double[] d, int dOffset, double[] e,
            int eOffset, double vl, double vu, int il, int iu, double abstol, intW mIW, double[] w, int wOffset,
            double[] z, int zOffset, int ldz, int[] isuppz, int isuppzOffset, double[] work, int workOffset, int lwork,
            int[] iwork, int iworkOffset, int liwork, boolean useCriticalRegion);

    private static native int dsyevr_n(int order, byte jobz, byte range, byte uplo, int n, double[] a, int aOffset,
            int lda, double vl, double vu, int il, int iu, double abstol, intW mIW, double[] w, int wOffset, double[] z,
            int zOffset, int ldz, int[] isuppz, int isuppzOffset, double[] work, int workOffset, int lwork, int[] iwork,
            int iworkOffset, int liwork, boolean useCriticalRegion);

    private static native int dsygvd_n(int order, int itype, byte jobz, byte uplo, int n, double[] a, int aOffset,
            int lda, double[] b, int bOffset, int ldb, double[] w, int wOffset, double[] work, int workOffset,
            int lwork, int[] iwork, int iworkOffset, int liwork, boolean useCriticalRegion);

    private static native int dsysv_n(int order, byte uplo, int n, int nrhs, double[] a, int aOffset, int lda,
            int[] ipiv, int ipivOffset, double[] b, int bOffset, int ldb, double[] work, int workOffset, int lwork,
            boolean useCriticalRegion);

    private static native int dtbtrs_n(int order, byte uplo, byte trans, byte diag, int n, int kd, int nrhs,
            double[] ab, int abOffset, int ldab, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dtptrs_n(int order, byte uplo, byte trans, byte diag, int n, int nrhs, double[] ap,
            int apOffset, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native int dtrtrs_n(int order, byte uplo, byte trans, byte diag, int n, int nrhs, double[] a,
            int aOffset, int lda, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    // miscellaneous float routines

    private static native int sgeev_n(int order, byte jobvl, byte jobvr, int n, float[] a, int aOffset, int lda,
            float[] wr, int wrOffset, float[] wi, int wiOffset, float[] vl, int vlOffset, int ldvl, float[] vr,
            int vrOffset, int ldvr, float[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int sgels_n(int order, byte trans, int m, int n, int nrhs, float[] a, int aOffset, int lda,
            float[] b, int bOffset, int ldb, float[] work, int workOffset, int lwork, boolean useCriticalRegion);

    private static native int sgesdd_n(int order, byte jobz, int m, int n, float[] a, int aOffset, int lda, float[] s,
            int sOffset, float[] u, int uOffset, int ldu, float[] vt, int vtOffset, int ldvt, float[] work,
            int workOffset, int lwork, int[] iwork, int iworkOffset, boolean useCriticalRegion);

    private static native int sgesv_n(int order, int n, int nrhs, float[] a, int aOffset, int lda, int[] ipiv,
            int ipivOffset, float[] b, int bOffset, int ldb, boolean useCriticalRegion);

    static native void initialize_n();

    private static byte norm(String norm) {
        char c = Character.toUpperCase(norm.charAt(0));
        if (c == '1' || c == 'I' || c == 'O') {
            return (byte) c;
        }
        throw new IllegalArgumentException("Invalid norm: " + norm);
    }

    private static byte trans(String trans) {
        char c = Character.toUpperCase(trans.charAt(0));
        if (c == 'N' || c == 'T' || c == 'C') {
            return (byte) c;
        }
        throw new IllegalArgumentException("Invalid transposition: " + trans);
    }

    private static byte eigJob(String job) {
        char c = Character.toUpperCase(job.charAt(0));
        if (c == 'V' || c == 'N') {
            return (byte) c;
        }
        throw new IllegalArgumentException("Invalid eigenvalue job: " + job);
    }

    private static byte range(String range) {
        char c = Character.toUpperCase(range.charAt(0));
        if (c == 'A' || c == 'I' || c == 'V') {
            return (byte) c;
        }
        throw new IllegalArgumentException("Invalid range: " + range);
    }

    private static byte svdJob(String job) {
        char c = Character.toUpperCase(job.charAt(0));
        if (c == 'A' || c == 'S' || c == 'O' || c == 'N') {
            return (byte) c;
        }
        throw new IllegalArgumentException("Invalid svd job: " + job);
    }

    private static byte uplo(String uplo) {
        char c = Character.toUpperCase(uplo.charAt(0));
        if (c == 'U' || c == 'L') {
            return (byte) c;
        }
        throw new IllegalArgumentException("Invalid uplo: " + uplo);
    }

    private static byte side(String side) {
        char c = Character.toUpperCase(side.charAt(0));
        if (c == 'L' || c == 'R') {
            return (byte) c;
        }
        throw new IllegalArgumentException("Invalid side: " + side);
    }

    private static byte diag(String diag) {
        char c = Character.toUpperCase(diag.charAt(0));
        if (c == 'N' || c == 'U') {
            return (byte) c;
        }
        throw new IllegalArgumentException("Invalid diag: " + diag);
    }

    protected LapackN() {
    }
}
