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
package net.dedekind.blas;

import org.netlib.blas.*;

/**
 * {@code Netlib} implementation.
 */
public class BlasJ extends Blas {

    @Override
    public final void dgbmv(String trans, int m, int n, int kl, int ku, double alpha, double[] a, int aOffset, int lda,
            double[] x, int xOffset, int incx, double beta, double[] y, int yOffset, int incy) {
        Dgbmv.dgbmv(trans, m, n, kl, ku, alpha, a, aOffset, lda, x, xOffset, incx, beta, y, yOffset, incy);
    }

    @Override
    public final void dgemm(String transa, String transb, int m, int n, int k, double alpha, double[] a, int aOffset,
            int lda, double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc) {
        Dgemm.dgemm(transa, transb, m, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, c, cOffset, ldc);
    }

    @Override
    public final void dgemv(String trans, int m, int n, double alpha, double[] a, int aOffset, int lda, double[] x,
            int xOffset, int incx, double beta, double[] y, int yOffset, int incy) {
        Dgemv.dgemv(trans, m, n, alpha, a, aOffset, lda, x, xOffset, incx, beta, y, yOffset, incy);
    }

    @Override
    public final void dger(int m, int n, double alpha, double[] x, int xOffset, int incx, double[] y, int yOffset,
            int incy, double[] a, int aOffset, int lda) {
        Dger.dger(m, n, alpha, x, xOffset, incx, y, yOffset, incy, a, aOffset, lda);
    }

    @Override
    public final void dsbmv(String uplo, int n, int k, double alpha, double[] a, int aOffset, int lda, double[] x,
            int xOffset, int incx, double beta, double[] y, int yOffset, int incy) {
        Dsbmv.dsbmv(uplo, n, k, alpha, a, aOffset, lda, x, xOffset, incx, beta, y, yOffset, incy);
    }

    @Override
    public final void dspmv(String uplo, int n, double alpha, double[] ap, int apOffset, double[] x, int xOffset,
            int incx, double beta, double[] y, int yOffset, int incy) {
        Dspmv.dspmv(uplo, n, alpha, ap, apOffset, x, xOffset, incx, beta, y, yOffset, incy);
    }

    @Override
    public final void dspr(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] ap,
            int apOffset) {
        Dspr.dspr(uplo, n, alpha, x, xOffset, incx, ap, apOffset);
    }

    @Override
    public final void dspr2(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] y,
            int yOffset, int incy, double[] ap, int apOffset) {
        Dspr2.dspr2(uplo, n, alpha, x, xOffset, incx, y, yOffset, incy, ap, apOffset);
    }

    @Override
    public final void dsymm(String side, String uplo, int m, int n, double alpha, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc) {
        Dsymm.dsymm(side, uplo, m, n, alpha, a, aOffset, lda, b, bOffset, ldb, beta, c, cOffset, ldc);
    }

    @Override
    public final void dsymv(String uplo, int n, double alpha, double[] a, int aOffset, int lda, double[] x, int xOffset,
            int incx, double beta, double[] y, int yOffset, int incy) {
        Dsymv.dsymv(uplo, n, alpha, a, aOffset, lda, x, xOffset, incx, beta, y, yOffset, incy);
    }

    @Override
    public final void dsyr(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] a, int aOffset,
            int lda) {
        Dsyr.dsyr(uplo, n, alpha, x, xOffset, incx, a, aOffset, lda);
    }

    @Override
    public final void dsyr2(String uplo, int n, double alpha, double[] x, int xOffset, int incx, double[] y,
            int yOffset, int incy, double[] a, int aOffset, int lda) {
        Dsyr2.dsyr2(uplo, n, alpha, x, xOffset, incx, y, yOffset, incy, a, aOffset, lda);
    }

    @Override
    public final void dsyr2k(String uplo, String trans, int n, int k, double alpha, double[] a, int aOffset, int lda,
            double[] b, int bOffset, int ldb, double beta, double[] c, int cOffset, int ldc) {
        Dsyr2k.dsyr2k(uplo, trans, n, k, alpha, a, aOffset, lda, b, bOffset, ldb, beta, c, cOffset, ldc);
    }

    @Override
    public final void dsyrk(String uplo, String trans, int n, int k, double alpha, double[] a, int aOffset, int lda,
            double beta, double[] c, int cOffset, int ldc) {
        Dsyrk.dsyrk(uplo, trans, n, k, alpha, a, aOffset, lda, beta, c, cOffset, ldc);
    }

    @Override
    public final void dtbmv(String uplo, String trans, String diag, int n, int k, double[] a, int aOffset, int lda,
            double[] x, int xOffset, int incx) {
        Dtbmv.dtbmv(uplo, trans, diag, n, k, a, aOffset, lda, x, xOffset, incx);
    }

    @Override
    public final void dtpmv(String uplo, String trans, String diag, int n, double[] ap, int apOffset, double[] x,
            int xOffset, int incx) {
        Dtpmv.dtpmv(uplo, trans, diag, n, ap, apOffset, x, xOffset, incx);
    }

    @Override
    public final void dtrmm(String side, String uplo, String transa, String diag, int m, int n, double alpha,
            double[] a, int aOffset, int lda, double[] b, int bOffset, int ldb) {
        Dtrmm.dtrmm(side, uplo, transa, diag, m, n, alpha, a, aOffset, lda, b, bOffset, ldb);
    }

    @Override
    public final void dtrmv(String uplo, String trans, String diag, int n, double[] a, int aOffset, int lda, double[] x,
            int xOffset, int incx) {
        Dtrmv.dtrmv(uplo, trans, diag, n, a, aOffset, lda, x, xOffset, incx);
    }

    @Override
    public void dtbsv(String uplo, String trans, String diag, int n, int k, double[] a, int aOffset, int lda,
            double[] x, int xOffset, int incx) {
        Dtbsv.dtbsv(uplo, trans, diag, n, k, a, aOffset, lda, x, xOffset, incx);

    }

    @Override
    public void dtpsv(String uplo, String trans, String diag, int n, double[] ap, int apOffset, double[] x, int xOffset,
            int incx) {
        Dtpsv.dtpsv(uplo, trans, diag, n, ap, apOffset, x, xOffset, incx);
    }

    @Override
    public void daxpy(int n, double da, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy) {
        Daxpy.daxpy(n, da, dx, dxOffset, incx, dy, dyOffset, incy);
    }

    @Override
    public void dcopy(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy) {
        Dcopy.dcopy(n, dx, dxOffset, incx, dy, dyOffset, incy);
    }

    @Override
    public void dscal(int n, double da, double[] dx, int dxOffset, int incx) {
        Dscal.dscal(n, da, dx, dxOffset, incx);
    }

    @Override
    public void dswap(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy) {
        Dswap.dswap(n, dx, dxOffset, incx, dy, dyOffset, incy);
    }

    @Override
    public double ddot(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy) {
        return Ddot.ddot(n, dx, dxOffset, incx, dy, dyOffset, incy);
    }

    @Override
    public void drot(int n, double[] dx, int dxOffset, int incx, double[] dy, int dyOffset, int incy, double c,
            double s) {
        Drot.drot(n, dx, dxOffset, incx, dy, dyOffset, incy, c, s);
    }

    @Override
    public int idamax(int n, double[] dx, int dxOffset, int incx) {
        return Idamax.idamax(n, dx, dxOffset, incx);
    }

    protected BlasJ() {
    }
}
