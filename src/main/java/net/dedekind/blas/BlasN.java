/*
 * Copyright 2019, 2023 Stefan Zobel
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

/**
 * {@code Intel MKL} (Math Kernel Library) native implementation.
 */
public class BlasN extends Blas {

    private static final boolean USE_CRITICAL = true;

    @Override
    public final void dgbmv(String trans, int m, int n, int kl, int ku, double alpha, double[] a, int aOffset, int lda,
            double[] x, int xOffset, int incx, double beta, double[] y, int yOffset, int incy) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(x, "x");
        Objects.requireNonNull(y, "y");
        dgbmv_n(Order.COL.code(), Trans.of(trans).code(), m, n, kl, ku, alpha, a, aOffset, lda, x, xOffset, incx, beta,
                y, yOffset, incy, USE_CRITICAL);
    }

    @Override
    public final void dgemm(String transa, String transb, int m, int n, int k, double alpha, double[] a, int aOffset,
            int lda, double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(c, "c");
        dgemm_n(Order.COL.code(), Trans.of(transa).code(), Trans.of(transb).code(), m, n, k, alpha, a, aOffset, lda, b,
                bOffset, ldb, beta, c, cOffset, ldc, USE_CRITICAL);
    }

    @Override
    public void dgemm_multi(String transa, String transb, int m, int n, int k, double alpha, double[] a, int aOffset,
            int lda, double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc, int howMany,
            int incAOff, int incBOff, int incCOff) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(c, "c");
        throw new UnsupportedOperationException("NOT YET IMPLEMENTED");
    }

    @Override
    public final void dgemv(String trans, int m, int n, double alpha, double[] a, int aOffset, int lda, double[] x,
            int xOffset, int incx, double beta, double[] y, int yOffset, int incy) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(x, "x");
        Objects.requireNonNull(y, "y");
        dgemv_n(Order.COL.code(), Trans.of(trans).code(), m, n, alpha, a, aOffset, lda, x, xOffset, incx, beta, y,
                yOffset, incy, USE_CRITICAL);
    }

    @Override
    public final void dger(int m, int n, double alpha, double[] x, int xOffset, int incx, double[] y, int yOffset,
            int incy, double[] a, int aOffset, int lda) {
        Objects.requireNonNull(x, "x");
        Objects.requireNonNull(y, "y");
        Objects.requireNonNull(a, "a");
        dger_n(Order.COL.code(), m, n, alpha, x, xOffset, incx, y, yOffset, incy, a, aOffset, lda, USE_CRITICAL);
    }

    @Override
    public final void dsbmv(String uplo, int n, int k, double alpha, double[] a, int aOffset, int lda, double[] x,
            int xOffset, int incx, double beta, double[] y, int yOffset, int incy) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(x, "x");
        Objects.requireNonNull(y, "y");
        dsbmv_n(Order.COL.code(), Uplo.of(uplo).code(), n, k, alpha, a, aOffset, lda, x, xOffset, incx, beta, y,
                yOffset, incy, USE_CRITICAL);
    }

    @Override
    public final void dspmv(String uplo, int n, double alpha, double[] ap, int apOffset, double[] x, int xOffset,
            int incx, double beta, double[] y, int yOffset, int incy) {
        Objects.requireNonNull(ap, "ap");
        Objects.requireNonNull(x, "x");
        Objects.requireNonNull(y, "y");
        dspmv_n(Order.COL.code(), Uplo.of(uplo).code(), n, alpha, ap, apOffset, x, xOffset, incx, beta, y, yOffset,
                incy, USE_CRITICAL);
    }

    @Override
    public final void dspr(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] ap,
            int apOffset) {
        Objects.requireNonNull(x, "x");
        Objects.requireNonNull(ap, "ap");
        dspr_n(Order.COL.code(), Uplo.of(uplo).code(), n, alpha, x, xOffset, incx, ap, apOffset, USE_CRITICAL);
    }

    @Override
    public final void dspr2(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] y,
            int yOffset, int incy, double[] ap, int apOffset) {
        Objects.requireNonNull(x, "x");
        Objects.requireNonNull(y, "y");
        Objects.requireNonNull(ap, "ap");
        dspr2_n(Order.COL.code(), Uplo.of(uplo).code(), n, alpha, x, xOffset, incx, y, yOffset, incy, ap, apOffset,
                USE_CRITICAL);
    }

    @Override
    public final void dsymm(String side, String uplo, int m, int n, double alpha, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(c, "c");
        dsymm_n(Order.COL.code(), Side.of(side).code(), Uplo.of(uplo).code(), m, n, alpha, a, aOffset, lda, b, bOffset,
                ldb, beta, c, cOffset, ldc, USE_CRITICAL);
    }

    @Override
    public final void dsymv(String uplo, int n, double alpha, double[] a, int aOffset, int lda, double[] x, int xOffset,
            int incx, double beta, double[] y, int yOffset, int incy) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(x, "x");
        Objects.requireNonNull(y, "y");
        dsymv_n(Order.COL.code(), Uplo.of(uplo).code(), n, alpha, a, aOffset, lda, x, xOffset, incx, beta, y, yOffset,
                incy, USE_CRITICAL);
    }

    @Override
    public final void dsyr(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] a, int aOffset,
            int lda) {
        Objects.requireNonNull(x, "x");
        Objects.requireNonNull(a, "a");
        dsyr_n(Order.COL.code(), Uplo.of(uplo).code(), n, alpha, x, xOffset, incx, a, aOffset, lda, USE_CRITICAL);
    }

    @Override
    public final void dsyr2(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] y,
            int yOffset, int incy, double[] a, int aOffset, int lda) {
        Objects.requireNonNull(x, "x");
        Objects.requireNonNull(y, "y");
        Objects.requireNonNull(a, "a");
        dsyr2_n(Order.COL.code(), Uplo.of(uplo).code(), n, alpha, x, xOffset, incx, y, yOffset, incy, a, aOffset, lda,
                USE_CRITICAL);
    }

    @Override
    public final void dsyr2k(String uplo, String trans, int n, int k, double alpha, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(c, "c");
        dsyr2k_n(Order.COL.code(), Uplo.of(uplo).code(), Trans.of(trans).code(), n, k, alpha, a, aOffset, lda, b,
                bOffset, ldb, beta, c, cOffset, ldc, USE_CRITICAL);
    }

    @Override
    public final void dsyrk(String uplo, String trans, int n, int k, double alpha, double[] a, int aOffset, int lda,
            double beta, double[] c, int cOffset, int ldc) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(c, "c");
        dsyrk_n(Order.COL.code(), Uplo.of(uplo).code(), Trans.of(trans).code(), n, k, alpha, a, aOffset, lda, beta, c,
                cOffset, ldc, USE_CRITICAL);
    }

    @Override
    public final void dtbmv(String uplo, String trans, String diag, int n, int k, double[] a, int aOffset, int lda,
            double[] x, int xOffset, int incx) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(x, "x");
        dtbmv_n(Order.COL.code(), Uplo.of(uplo).code(), Trans.of(trans).code(), Diag.of(diag).code(), n, k, a, aOffset,
                lda, x, xOffset, incx, USE_CRITICAL);
    }

    @Override
    public final void dtpmv(String uplo, String trans, String diag, int n, double[] ap, int apOffset, double[] x,
            int xOffset, int incx) {
        Objects.requireNonNull(ap, "ap");
        Objects.requireNonNull(x, "x");
        dtpmv_n(Order.COL.code(), Uplo.of(uplo).code(), Trans.of(trans).code(), Diag.of(diag).code(), n, ap, apOffset,
                x, xOffset, incx, USE_CRITICAL);
    }

    @Override
    public final void dtrmm(String side, String uplo, String transa, String diag, int m, int n, double alpha,
            double[] a, int aOffset, int lda, double[] b, int bOffset, int ldb) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        dtrmm_n(Order.COL.code(), Side.of(side).code(), Uplo.of(uplo).code(), Trans.of(transa).code(),
                Diag.of(diag).code(), m, n, alpha, a, aOffset, lda, b, bOffset, ldb, USE_CRITICAL);
    }

    @Override
    public final void dtrmv(String uplo, String trans, String diag, int n, double[] a, int aOffset, int lda, double[] x,
            int xOffset, int incx) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(x, "x");
        dtrmv_n(Order.COL.code(), Uplo.of(uplo).code(), Trans.of(trans).code(), Diag.of(diag).code(), n, a, aOffset,
                lda, x, xOffset, incx, USE_CRITICAL);
    }

    @Override
    public final void dtbsv(String uplo, String trans, String diag, int n, int k, double[] a, int aOffset, int lda,
            double[] x, int xOffset, int incx) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(x, "x");
        dtbsv_n(Order.COL.code(), Uplo.of(uplo).code(), Trans.of(trans).code(), Diag.of(diag).code(), n, k, a, aOffset,
                lda, x, xOffset, incx, USE_CRITICAL);
    }

    @Override
    public final void dtpsv(String uplo, String trans, String diag, int n, double[] ap, int apOffset, double[] x,
            int xOffset, int incx) {
        Objects.requireNonNull(ap, "ap");
        Objects.requireNonNull(x, "x");
        dtpsv_n(Order.COL.code(), Uplo.of(uplo).code(), Trans.of(trans).code(), Diag.of(diag).code(), n, ap, apOffset,
                x, xOffset, incx, USE_CRITICAL);
    }

    @Override
    public final void daxpy(int n, double da, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset,
            int incy) {
        Objects.requireNonNull(dx, "dx");
        Objects.requireNonNull(dy, "dy");
        daxpy_n(n, da, dx, dxOffset, incx, dy, dyOffset, incy, USE_CRITICAL);
    }

    @Override
    public final void dcopy(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy) {
        Objects.requireNonNull(dx, "dx");
        Objects.requireNonNull(dy, "dy");
        dcopy_n(n, dx, dxOffset, incx, dy, dyOffset, incy, USE_CRITICAL);
    }

    @Override
    public final void dscal(int n, double da, double[] dx, int dxOffset, int incx) {
        Objects.requireNonNull(da, "da");
        Objects.requireNonNull(dx, "dx");
        dscal_n(n, da, dx, dxOffset, incx, USE_CRITICAL);
    }

    @Override
    public final void dswap(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy) {
        Objects.requireNonNull(dx, "dx");
        Objects.requireNonNull(dy, "dy");
        dswap_n(n, dx, dxOffset, incx, dy, dyOffset, incy, USE_CRITICAL);
    }

    @Override
    public final double ddot(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy) {
        Objects.requireNonNull(dx, "dx");
        Objects.requireNonNull(dy, "dy");
        return ddot_n(n, dx, dxOffset, incx, dy, dyOffset, incy, USE_CRITICAL);
    }

    @Override
    public final void drot(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy, double c,
            double s) {
        Objects.requireNonNull(dx, "dx");
        Objects.requireNonNull(dy, "dy");
        drot_n(n, dx, dxOffset, incx, dy, dyOffset, incy, c, s, USE_CRITICAL);
    }

    @Override
    public final int idamax(int n, double[] dx, int dxOffset, int incx) {
        Objects.requireNonNull(dx, "dx");
        return idamax_n(n, dx, dxOffset, incx, USE_CRITICAL);
    }

    // miscellaneous float routines

    @Override
    public final void sgemm(String transa, String transb, int m, int n, int k, float alpha, float[] a, int aOffset,
            int lda, float[] b, int bOffset, int ldb, float beta, float[] c, int cOffset, int ldc) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(c, "c");
        sgemm_n(Order.COL.code(), Trans.of(transa).code(), Trans.of(transb).code(), m, n, k, alpha, a, aOffset, lda, b,
                bOffset, ldb, beta, c, cOffset, ldc, USE_CRITICAL);
    }

    @Override
    public void sgemm_multi(String transa, String transb, int m, int n, int k, float alpha, float[] a, int aOffset,
            int lda, float[] b, int bOffset, int ldb, float beta, float[] c, int cOffset, int ldc, int howMany,
            int incAOff, int incBOff, int incCOff) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(c, "c");
        throw new UnsupportedOperationException("NOT YET IMPLEMENTED");
    }

    // miscellaneous complex routines

    @Override
    public final void cgemm(Trans transa, Trans transb, int m, int n, int k, float alphar, float alphai, float[] a,
            int lda, float[] b, int ldb, float betar, float betai, float[] c, int ldc) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(c, "c");
        Objects.requireNonNull(transa, "transa");
        Objects.requireNonNull(transb, "transb");
        cgemm_n(Order.COL.code(), transa.code(), transb.code(), m, n, k, alphar, alphai, a, lda, b, ldb, betar, betai,
                c, ldc, USE_CRITICAL);
    }

    @Override
    public final void zgemm(Trans transa, Trans transb, int m, int n, int k, double alphar, double alphai, double[] a,
            int lda, double[] b, int ldb, double betar, double betai, double[] c, int ldc) {
        Objects.requireNonNull(a, "a");
        Objects.requireNonNull(b, "b");
        Objects.requireNonNull(c, "c");
        Objects.requireNonNull(transa, "transa");
        Objects.requireNonNull(transb, "transb");
        zgemm_n(Order.COL.code(), transa.code(), transb.code(), m, n, k, alphar, alphai, a, lda, b, ldb, betar, betai,
                c, ldc, USE_CRITICAL);
    }

    // native methods

    private static native void dgbmv_n(int order, int trans, int m, int n, int kl, int ku, double alpha, double[] a,
            int aOffset, int lda, double[] x, int xOffset, int incx, double beta, double[] y, int yOffset, int incy,
            boolean useCriticalRegion);

    private static native void dgemm_n(int order, int transa, int transb, int m, int n, int k, double alpha, double[] a,
            int aOffset, int lda, double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc,
            boolean useCriticalRegion);

    private static native void dgemm_multi_n(int order, int transa, int transb, int m, int n, int k, double alpha,
            double[] a, int aOffset, int lda, double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset,
            int ldc, boolean useCriticalRegion, int howMany, int incAOff, int incBOff, int incCOff);

    private static native void dgemv_n(int order, int trans, int m, int n, double alpha, double[] a, int aOffset,
            int lda, double[] x, int xOffset, int incx, double beta, double[] y, int yOffset, int incy,
            boolean useCriticalRegion);

    private static native void dger_n(int order, int m, int n, double alpha, double[] x, int xOffset, int incx,
            double[] y, int yOffset, int incy, double[] a, int aOffset, int lda, boolean useCriticalRegion);

    private static native void dsbmv_n(int order, int uplo, int n, int k, double alpha, double[] a, int aOffset,
            int lda, double[] x, int xOffset, int incx, double beta, double[] y, int yOffset, int incy,
            boolean useCriticalRegion);

    private static native void dspmv_n(int order, int uplo, int n, double alpha, double[] ap, int apOffset, double[] x,
            int xOffset, int incx, double beta, double[] y, int yOffset, int incy, boolean useCriticalRegion);

    private static native void dspr_n(int order, int uplo, int n, double alpha, double[] x, int xOffset, int incx,
            double[] ap, int apOffset, boolean useCriticalRegion);

    private static native void dspr2_n(int order, int uplo, int n, double alpha, double[] x, int xOffset, int incx,
            double[] y, int yOffset, int incy, double[] ap, int apOffset, boolean useCriticalRegion);

    private static native void dsymm_n(int order, int side, int uplo, int m, int n, double alpha, double[] a,
            int aOffset, int lda, double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc,
            boolean useCriticalRegion);

    private static native void dsymv_n(int order, int uplo, int n, double alpha, double[] a, int aOffset, int lda,
            double[] x, int xOffset, int incx, double beta, double[] y, int yOffset, int incy,
            boolean useCriticalRegion);

    private static native void dsyr_n(int order, int uplo, int n, double alpha, double[] x, int xOffset, int incx,
            double[] a, int aOffset, int lda, boolean useCriticalRegion);

    private static native void dsyr2_n(int order, int uplo, int n, double alpha, double[] x, int xOffset, int incx,
            double[] y, int yOffset, int incy, double[] a, int aOffset, int lda, boolean useCriticalRegion);

    private static native void dsyr2k_n(int order, int uplo, int trans, int n, int k, double alpha, double[] a,
            int aOffset, int lda, double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc,
            boolean useCriticalRegion);

    private static native void dsyrk_n(int order, int uplo, int trans, int n, int k, double alpha, double[] a,
            int aOffset, int lda, double beta, double[] c, int cOffset, int ldc, boolean useCriticalRegion);

    private static native void dtbmv_n(int order, int uplo, int trans, int diag, int n, int k, double[] a, int aOffset,
            int lda, double[] x, int xOffset, int incx, boolean useCriticalRegion);

    private static native void dtpmv_n(int order, int uplo, int trans, int diag, int n, double[] ap, int apOffset,
            double[] x, int xOffset, int incx, boolean useCriticalRegion);

    private static native void dtrmm_n(int order, int side, int uplo, int transa, int diag, int m, int n, double alpha,
            double[] a, int aOffset, int lda, double[] b, int bOffset, int ldb, boolean useCriticalRegion);

    private static native void dtrmv_n(int order, int uplo, int trans, int diag, int n, double[] a, int aOffset,
            int lda, double[] x, int xOffset, int incx, boolean useCriticalRegion);

    private static native void dtbsv_n(int order, int uplo, int trans, int diag, int n, int k, double[] a, int aOffset,
            int lda, double[] x, int xOffset, int incx, boolean useCriticalRegion);

    private static native void dtpsv_n(int order, int uplo, int trans, int diag, int n, double[] ap, int apOffset,
            double[] x, int xOffset, int incx, boolean useCriticalRegion);

    private static native void daxpy_n(int n, double da, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset,
            int incy, boolean useCriticalRegion);

    private static native void dcopy_n(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy,
            boolean useCriticalRegion);

    private static native void dscal_n(int n, double da, double[] dx, int dxOffset, int incx,
            boolean useCriticalRegion);

    private static native void dswap_n(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy,
            boolean useCriticalRegion);

    private static native double ddot_n(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy,
            boolean useCriticalRegion);

    private static native void drot_n(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy,
            double c, double s, boolean useCriticalRegion);

    private static native int idamax_n(int n, double[] dx, int dxOffset, int incx, boolean useCriticalRegion);

    static native void redirect_xerbla_n();

    // miscellaneous float routines

    private static native void sgemm_n(int order, int transa, int transb, int m, int n, int k, float alpha, float[] a,
            int aOffset, int lda, float[] b, int bOffset, int ldb, float beta, float[] c, int cOffset, int ldc,
            boolean useCriticalRegion);

    private static native void sgemm_multi_n(int order, int transa, int transb, int m, int n, int k, float alpha,
            float[] a, int aOffset, int lda, float[] b, int bOffset, int ldb, float beta, float[] c, int cOffset,
            int ldc, boolean useCriticalRegion, int howMany, int incAOff, int incBOff, int incCOff);

    // miscellaneous complex routines

    private static native void cgemm_n(int order, int transa, int transb, int m, int n, int k, float alphar,
            float alphai, float[] a, int lda, float[] b, int ldb, float betar, float betai, float[] c, int ldc,
            boolean useCriticalRegion);

    private static native void zgemm_n(int order, int transa, int transb, int m, int n, int k, double alphar,
            double alphai, double[] a, int lda, double[] b, int ldb, double betar, double betai, double[] c, int ldc,
            boolean useCriticalRegion);

    protected BlasN() {
    }
}
