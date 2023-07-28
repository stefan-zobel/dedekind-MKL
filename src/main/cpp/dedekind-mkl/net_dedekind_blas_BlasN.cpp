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

#ifndef _JAVASOFT_JNI_H_
#include <jni.h>
#endif /* _JAVASOFT_JNI_H_ */

#ifndef DEDEKIND_INCLUDED_
#include "dedekind-mkl.h"
#endif /* DEDEKIND_INCLUDED_ */

#ifndef DOUBLEARRAY_INCLUDED_
#include "DoubleArray.h"
#endif /* DOUBLEARRAY_INCLUDED_ */

#ifndef FLOATARRAY_INCLUDED_
#include "FloatArray.h"
#endif /* FLOATARRAY_INCLUDED_ */

#ifndef COMPLEX_FLOATARRAY_INCLUDED_
#include "ComplexFloatArray.h"
#endif /* COMPLEX_FLOATARRAY_INCLUDED_ */

#ifndef COMPLEX_DOUBLEARRAY_INCLUDED_
#include "ComplexDoubleArray.h"
#endif /* COMPLEX_DOUBLEARRAY_INCLUDED_ */

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */

#ifndef JEXCEPTIONUTILS_INCLUDED_
#include "JExceptionUtils.h"
#endif /* JEXCEPTIONUTILS_INCLUDED_ */

#ifndef THREADAUTOATTACH_INCLUDED_
#include "ThreadAutoAttach.h"
#endif /* THREADAUTOATTACH_INCLUDED_ */

#ifndef _MKL_H_
#include <mkl.h>
#endif /* _MKL_H_ */


#ifdef __cplusplus
extern "C" {
#endif

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dgbmv_n
 * Signature: (IIIIIID[DII[DIID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dgbmv_1n(JNIEnv* env, jclass,
  jint order,
  jint trans,
  jint m,
  jint n,
  jint kl,
  jint ku,
  jdouble alpha,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jdouble beta,
  jdoubleArray y,
  jint yOffset,
  jint incy,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);
        DoubleArray ya = DoubleArray(env, y, yOffset, useCrit);

        cblas_dgbmv(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_TRANSPOSE>(trans),
            m, n, kl, ku, alpha, aa.ptr(), lda, xa.ptr(), incx, beta, ya.ptr(), incy);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgbmv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgbmv_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dgemm_n
 * Signature: (IIIIIID[DII[DIID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dgemm_1n(JNIEnv* env, jclass,
  jint order,
  jint transa,
  jint transb,
  jint m,
  jint n,
  jint k,
  jdouble alpha,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jdouble beta,
  jdoubleArray c,
  jint cOffset,
  jint ldc,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);
        DoubleArray ca = DoubleArray(env, c, cOffset, useCrit);

        cblas_dgemm(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_TRANSPOSE>(transa),
            static_cast<CBLAS_TRANSPOSE>(transb), m, n, k, alpha, aa.ptr(), lda, ba.ptr(),
            ldb, beta, ca.ptr(), ldc);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgemm_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgemm_n: caught unknown exception");
    }
}
/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dgemm_multi_n
 * Signature: (IIIIIID[DII[DIID[DIIZIIII)V
 */
JNIEXPORT void JNICALL Java_net_dedekind_blas_BlasN_dgemm_1multi_1n(JNIEnv* env, jclass,
  jint order,
  jint transa,
  jint transb,
  jint m,
  jint n,
  jint k,
  jdouble alpha,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jdouble beta,
  jdoubleArray c,
  jint cOffset,
  jint ldc,
  jboolean useCrit,
  jint howMany,
  jint incAOff,
  jint incBOff,
  jint incCOff) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);
        DoubleArray ca = DoubleArray(env, c, cOffset, useCrit);

        double* pa = aa.ptr();
        double* pb = ba.ptr();
        double* pc = ca.ptr();

        for (int i = 0; i < howMany; ++i) {
            cblas_dgemm(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_TRANSPOSE>(transa),
                static_cast<CBLAS_TRANSPOSE>(transb), m, n, k, alpha, pa, lda, pb, ldb, beta,
                pc, ldc);

            pa += incAOff;
            pb += incBOff;
            pc += incCOff;
        }

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgemm_multi_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgemm_multi_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dgemv_n
 * Signature: (IIIID[DII[DIID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dgemv_1n(JNIEnv* env, jclass,
  jint order,
  jint trans,
  jint m,
  jint n,
  jdouble alpha,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jdouble beta,
  jdoubleArray y,
  jint yOffset,
  jint incy,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);
        DoubleArray ya = DoubleArray(env, y, yOffset, useCrit);

        cblas_dgemv(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_TRANSPOSE>(trans),
            m, n, alpha, aa.ptr(), lda, xa.ptr(), incx, beta, ya.ptr(), incy);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgemv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgemv_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dger_n
 * Signature: (IIID[DII[DII[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dger_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jdouble alpha,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jdoubleArray y,
  jint yOffset,
  jint incy,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jboolean useCrit) {
    try {
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);
        DoubleArray ya = DoubleArray(env, y, yOffset, useCrit);
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);

        cblas_dger(static_cast<CBLAS_LAYOUT>(order), m, n, alpha, xa.ptr(), incx,
            ya.ptr(), incy, aa.ptr(), lda);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dger_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dger_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dsbmv_n
 * Signature: (IIIID[DII[DIID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dsbmv_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint n,
  jint k,
  jdouble alpha,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jdouble beta,
  jdoubleArray y,
  jint yOffset,
  jint incy,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);
        DoubleArray ya = DoubleArray(env, y, yOffset, useCrit);

        cblas_dsbmv(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            n, k, alpha, aa.ptr(), lda, xa.ptr(), incx, beta, ya.ptr(), incy);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsbmv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsbmv_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dspmv_n
 * Signature: (IIID[DI[DIID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dspmv_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint n,
  jdouble alpha,
  jdoubleArray ap,
  jint apOffset,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jdouble beta,
  jdoubleArray y,
  jint yOffset,
  jint incy,
  jboolean useCrit) {
    try {
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);
        DoubleArray ya = DoubleArray(env, y, yOffset, useCrit);

        cblas_dspmv(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            n, alpha, apa.ptr(), xa.ptr(), incx, beta, ya.ptr(), incy);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dspmv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dspmv_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dspr_n
 * Signature: (IIID[DII[DIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dspr_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint n,
  jdouble alpha,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jdoubleArray ap,
  jint apOffset,
  jboolean useCrit) {
    try {
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);

        cblas_dspr(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            n, alpha, xa.ptr(), incx, apa.ptr());

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dspr_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dspr_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dspr2_n
 * Signature: (IIID[DII[DII[DIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dspr2_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint n,
  jdouble alpha,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jdoubleArray y,
  jint yOffset,
  jint incy,
  jdoubleArray ap,
  jint apOffset,
  jboolean useCrit) {
    try {
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);
        DoubleArray ya = DoubleArray(env, y, yOffset, useCrit);
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);

        cblas_dspr2(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            n, alpha, xa.ptr(), incx, ya.ptr(), incy, apa.ptr());

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dspr2_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dspr2_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dsymm_n
 * Signature: (IIIIID[DII[DIID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dsymm_1n(JNIEnv* env, jclass,
  jint order,
  jint side,
  jint uplo,
  jint m,
  jint n,
  jdouble alpha,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jdouble beta,
  jdoubleArray c,
  jint cOffset,
  jint ldc,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);
        DoubleArray ca = DoubleArray(env, c, cOffset, useCrit);

        cblas_dsymm(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_SIDE>(side),
            static_cast<CBLAS_UPLO>(uplo), m, n, alpha, aa.ptr(), lda, ba.ptr(), ldb,
            beta, ca.ptr(), ldc);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsymm_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsymm_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dsymv_n
 * Signature: (IIID[DII[DIID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dsymv_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint n,
  jdouble alpha,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jdouble beta,
  jdoubleArray y,
  jint yOffset,
  jint incy,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);
        DoubleArray ya = DoubleArray(env, y, yOffset, useCrit);

        cblas_dsymv(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            n, alpha, aa.ptr(), lda, xa.ptr(), incx, beta, ya.ptr(), incy);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsymv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsymv_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dsyr_n
 * Signature: (IIID[DII[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dsyr_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint n,
  jdouble alpha,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jboolean useCrit) {
    try {
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);

        cblas_dsyr(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            n, alpha, xa.ptr(), incx, aa.ptr(), lda);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsyr_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsyr_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dsyr2_n
 * Signature: (IIID[DII[DII[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dsyr2_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint n,
  jdouble alpha,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jdoubleArray y,
  jint yOffset,
  jint incy,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jboolean useCrit) {
    try {
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);
        DoubleArray ya = DoubleArray(env, y, yOffset, useCrit);
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);

        cblas_dsyr2(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            n, alpha, xa.ptr(), incx, ya.ptr(), incy, aa.ptr(), lda);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsyr2_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsyr2_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dsyr2k_n
 * Signature: (IIIIID[DII[DIID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dsyr2k_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint trans,
  jint n,
  jint k,
  jdouble alpha,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jdouble beta,
  jdoubleArray c,
  jint cOffset,
  jint ldc,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);
        DoubleArray ca = DoubleArray(env, c, cOffset, useCrit);

        cblas_dsyr2k(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            static_cast<CBLAS_TRANSPOSE>(trans), n, k, alpha, aa.ptr(), lda,
            ba.ptr(), ldb, beta, ca.ptr(), ldc);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsyr2k_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsyr2k_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dsyrk_n
 * Signature: (IIIIID[DIID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dsyrk_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint trans,
  jint n,
  jint k,
  jdouble alpha,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdouble beta,
  jdoubleArray c,
  jint cOffset,
  jint ldc,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray ca = DoubleArray(env, c, cOffset, useCrit);

        cblas_dsyrk(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            static_cast<CBLAS_TRANSPOSE>(trans), n, k, alpha, aa.ptr(), lda,
            beta, ca.ptr(), ldc);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsyrk_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsyrk_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dtbmv_n
 * Signature: (IIIIII[DII[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dtbmv_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint trans,
  jint diag,
  jint n,
  jint k,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);

        cblas_dtbmv(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            static_cast<CBLAS_TRANSPOSE>(trans), static_cast<CBLAS_DIAG>(diag),
            n, k, aa.ptr(), lda, xa.ptr(), incx);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dtbmv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dtbmv_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dtpmv_n
 * Signature: (IIIII[DI[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dtpmv_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint trans,
  jint diag,
  jint n,
  jdoubleArray ap,
  jint apOffset,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jboolean useCrit) {
    try {
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);

        cblas_dtpmv(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            static_cast<CBLAS_TRANSPOSE>(trans), static_cast<CBLAS_DIAG>(diag),
            n, apa.ptr(), xa.ptr(), incx);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dtpmv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dtpmv_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dtrmm_n
 * Signature: (IIIIIIID[DII[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dtrmm_1n(JNIEnv* env, jclass,
  jint order,
  jint side,
  jint uplo,
  jint transa,
  jint diag,
  jint m,
  jint n,
  jdouble alpha,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        cblas_dtrmm(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_SIDE>(side),
            static_cast<CBLAS_UPLO>(uplo), static_cast<CBLAS_TRANSPOSE>(transa),
            static_cast<CBLAS_DIAG>(diag), m, n, alpha, aa.ptr(), lda, ba.ptr(), ldb);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dtrmm_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dtrmm_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dtrmv_n
 * Signature: (IIIII[DII[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dtrmv_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint trans,
  jint diag,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);

        cblas_dtrmv(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            static_cast<CBLAS_TRANSPOSE>(trans), static_cast<CBLAS_DIAG>(diag),
            n, aa.ptr(), lda, xa.ptr(), incx);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dtrmv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dtrmv_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dtbsv_n
 * Signature: (IIIIII[DII[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dtbsv_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint trans,
  jint diag,
  jint n,
  jint k,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);

        cblas_dtbsv(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            static_cast<CBLAS_TRANSPOSE>(trans), static_cast<CBLAS_DIAG>(diag),
            n, k, aa.ptr(), lda, xa.ptr(), incx);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dtbsv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dtbsv_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dtpsv_n
 * Signature: (IIIII[DI[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dtpsv_1n(JNIEnv* env, jclass,
  jint order,
  jint uplo,
  jint trans,
  jint diag,
  jint n,
  jdoubleArray ap,
  jint apOffset,
  jdoubleArray x,
  jint xOffset,
  jint incx,
  jboolean useCrit) {
    try {
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);
        DoubleArray xa = DoubleArray(env, x, xOffset, useCrit);

        cblas_dtpsv(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_UPLO>(uplo),
            static_cast<CBLAS_TRANSPOSE>(trans), static_cast<CBLAS_DIAG>(diag),
            n, apa.ptr(), xa.ptr(), incx);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dtpsv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dtpsv_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    daxpy_n
 * Signature: (ID[DII[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_daxpy_1n(JNIEnv* env, jclass,
  jint n,
  jdouble da,
  jdoubleArray dx,
  jint dxOffset,
  jint incx,
  jdoubleArray dy,
  jint dyOffset,
  jint incy,
  jboolean useCrit) {
    try {
        DoubleArray dxa = DoubleArray(env, dx, dxOffset, useCrit);
        DoubleArray dya = DoubleArray(env, dy, dyOffset, useCrit);

        cblas_daxpy(n, da, dxa.ptr(), incx, dya.ptr(), incy);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "daxpy_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "daxpy_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dcopy_n
 * Signature: (I[DII[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dcopy_1n(JNIEnv* env, jclass,
  jint n,
  jdoubleArray dx,
  jint dxOffset,
  jint incx,
  jdoubleArray dy,
  jint dyOffset,
  jint incy,
  jboolean useCrit) {
    try {
        DoubleArray dxa = DoubleArray(env, dx, dxOffset, useCrit);
        DoubleArray dya = DoubleArray(env, dy, dyOffset, useCrit);

        cblas_dcopy(n, dxa.ptr(), incx, dya.ptr(), incy);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dcopy_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dcopy_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dscal_n
 * Signature: (ID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dscal_1n(JNIEnv* env, jclass,
  jint n,
  jdouble da,
  jdoubleArray dx,
  jint dxOffset,
  jint incx,
  jboolean useCrit) {
    try {
        DoubleArray dxa = DoubleArray(env, dx, dxOffset, useCrit);

        cblas_dscal(n, da, dxa.ptr(), incx);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dscal_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dscal_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    dswap_n
 * Signature: (I[DII[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_dswap_1n(JNIEnv* env, jclass,
  jint n,
  jdoubleArray dx,
  jint dxOffset,
  jint incx,
  jdoubleArray dy,
  jint dyOffset,
  jint incy,
  jboolean useCrit) {
    try {
        DoubleArray dxa = DoubleArray(env, dx, dxOffset, useCrit);
        DoubleArray dya = DoubleArray(env, dy, dyOffset, useCrit);

        cblas_dswap(n, dxa.ptr(), incx, dya.ptr(), incy);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dswap_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dswap_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    ddot_n
 * Signature: (I[DII[DIIZ)D
 */
JNIEXPORT jdouble JNICALL
Java_net_dedekind_blas_BlasN_ddot_1n(JNIEnv* env, jclass,
  jint n,
  jdoubleArray dx,
  jint dxOffset,
  jint incx,
  jdoubleArray dy,
  jint dyOffset,
  jint incy,
  jboolean useCrit) {
    try {
        DoubleArray dxa = DoubleArray(env, dx, dxOffset, useCrit);
        DoubleArray dya = DoubleArray(env, dy, dyOffset, useCrit);

        return cblas_ddot(n, dxa.ptr(), incx, dya.ptr(), incy);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "ddot_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "ddot_n: caught unknown exception");
    }
    return 0.0;
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    drot_n
 * Signature: (I[DII[DIIDDZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_drot_1n(JNIEnv* env, jclass,
  jint n,
  jdoubleArray dx,
  jint dxOffset,
  jint incx,
  jdoubleArray dy,
  jint dyOffset,
  jint incy,
  jdouble c,
  jdouble s,
  jboolean useCrit) {
    try {
        DoubleArray dxa = DoubleArray(env, dx, dxOffset, useCrit);
        DoubleArray dya = DoubleArray(env, dy, dyOffset, useCrit);

        cblas_drot(n, dxa.ptr(), incx, dya.ptr(), incy, c, s);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "drot_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "drot_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    idamax_n
 * Signature: (I[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_blas_BlasN_idamax_1n(JNIEnv* env, jclass,
  jint n,
  jdoubleArray dx,
  jint dxOffset,
  jint incx,
  jboolean useCrit) {
    try {
        DoubleArray dxa = DoubleArray(env, dx, dxOffset, useCrit);

        return static_cast<jint>(cblas_idamax(n, dxa.ptr(), incx));

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "idamax_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "idamax_n: caught unknown exception");
    }
    return -1;
}

    // miscellaneous float routines

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    sgemm_n
 * Signature: (IIIIIIF[FII[FIIF[FIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_sgemm_1n(JNIEnv* env, jclass,
  jint order,
  jint transa,
  jint transb,
  jint m,
  jint n,
  jint k,
  jfloat alpha,
  jfloatArray a,
  jint aOffset,
  jint lda,
  jfloatArray b,
  jint bOffset,
  jint ldb,
  jfloat beta,
  jfloatArray c,
  jint cOffset,
  jint ldc,
  jboolean useCrit) {
    try {
        FloatArray aa = FloatArray(env, a, aOffset, useCrit);
        FloatArray ba = FloatArray(env, b, bOffset, useCrit);
        FloatArray ca = FloatArray(env, c, cOffset, useCrit);

        cblas_sgemm(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_TRANSPOSE>(transa),
            static_cast<CBLAS_TRANSPOSE>(transb), m, n, k, alpha, aa.ptr(), lda, ba.ptr(),
            ldb, beta, ca.ptr(), ldc);

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "sgemm_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "sgemm_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    sgemm_multi_n
 * Signature: (IIIIIIF[FII[FIIF[FIIZIIII)V
 */
JNIEXPORT void JNICALL Java_net_dedekind_blas_BlasN_sgemm_1multi_1n(JNIEnv* env, jclass,
  jint order,
  jint transa,
  jint transb,
  jint m,
  jint n,
  jint k,
  jfloat alpha,
  jfloatArray a,
  jint aOffset,
  jint lda,
  jfloatArray b,
  jint bOffset,
  jint ldb,
  jfloat beta,
  jfloatArray c,
  jint cOffset,
  jint ldc,
  jboolean useCrit,
  jint howMany,
  jint incAOff,
  jint incBOff,
  jint incCOff) {
    try {
        FloatArray aa = FloatArray(env, a, aOffset, useCrit);
        FloatArray ba = FloatArray(env, b, bOffset, useCrit);
        FloatArray ca = FloatArray(env, c, cOffset, useCrit);

        float* pa = aa.ptr();
        float* pb = ba.ptr();
        float* pc = ca.ptr();

        for (int i = 0; i < howMany; ++i) {
            cblas_sgemm(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_TRANSPOSE>(transa),
                static_cast<CBLAS_TRANSPOSE>(transb), m, n, k, alpha, pa, lda, pb, ldb, beta,
                pc, ldc);

            pa += incAOff;
            pb += incBOff;
            pc += incCOff;
        }

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "sgemm_multi_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "sgemm_multi_n: caught unknown exception");
    }
}

    // miscellaneous complex routines

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    cgemm_n
 * Signature: (IIIIIIFF[FI[FIFF[FIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_cgemm_1n(JNIEnv* env, jclass,
  jint order,
  jint transa,
  jint transb,
  jint m,
  jint n,
  jint k,
  jfloat alphar,
  jfloat alphai,
  jfloatArray a,
  jint lda,
  jfloatArray b,
  jint ldb,
  jfloat betar,
  jfloat betai,
  jfloatArray c,
  jint ldc,
  jboolean useCrit) {
    try {
        FloatArray aa = FloatArray(env, a, 0, useCrit);
        FloatArray ba = FloatArray(env, b, 0, useCrit);
        FloatArray ca = FloatArray(env, c, 0, useCrit);

        ComplexFloatArray aac = ComplexFloatArray(aa);
        ComplexFloatArray bac = ComplexFloatArray(ba);
        ComplexFloatArray cac = ComplexFloatArray(ca);

        MKL_Complex8* pa = aac.ptr();
        MKL_Complex8* pb = bac.ptr();
        MKL_Complex8* pc = cac.ptr();

        if (pa && pb && pc) {
            MKL_Complex8 alpha = { alphar, alphai };
            MKL_Complex8 beta = { betar, betai };

            cblas_cgemm(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_TRANSPOSE>(transa),
                static_cast<CBLAS_TRANSPOSE>(transb), m, n, k, &alpha, pa, lda, pb, ldb, &beta, pc, ldc);

            long len = cac.complexLength();
            if (len > 0 && cac.hasCopy()) {
                floatCopy(len, ca.ptr(), pc);
            }
        }

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "cgemm_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "cgemm_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    zgemm_n
 * Signature: (IIIIIIDD[DI[DIDD[DIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_zgemm_1n(JNIEnv* env, jclass,
  jint order,
  jint transa,
  jint transb,
  jint m,
  jint n,
  jint k,
  jdouble alphar,
  jdouble alphai,
  jdoubleArray a,
  jint lda,
  jdoubleArray b,
  jint ldb,
  jdouble betar,
  jdouble betai,
  jdoubleArray c,
  jint ldc,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, 0, useCrit);
        DoubleArray ba = DoubleArray(env, b, 0, useCrit);
        DoubleArray ca = DoubleArray(env, c, 0, useCrit);

        ComplexDoubleArray aac = ComplexDoubleArray(aa);
        ComplexDoubleArray bac = ComplexDoubleArray(ba);
        ComplexDoubleArray cac = ComplexDoubleArray(ca);

        MKL_Complex16* pa = aac.ptr();
        MKL_Complex16* pb = bac.ptr();
        MKL_Complex16* pc = cac.ptr();

        if (pa && pb && pc) {
            MKL_Complex16 alpha = { alphar, alphai };
            MKL_Complex16 beta = { betar, betai };

            cblas_zgemm(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_TRANSPOSE>(transa),
                static_cast<CBLAS_TRANSPOSE>(transb), m, n, k, &alpha, pa, lda, pb, ldb, &beta, pc, ldc);

            long len = cac.complexLength();
            if (len > 0 && cac.hasCopy()) {
                doubleCopy(len, ca.ptr(), pc);
            }
        }

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "zgemm_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "zgemm_n: caught unknown exception");
    }
}

    // xerbla

/*
 * Class:     net_dedekind_blas_BlasN
 * Method:    redirect_xerbla_n
 * Signature: ()V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasN_redirect_1xerbla_1n(JNIEnv*, jclass) {
    mkl_set_xerbla(xerbla);
}

#if defined (_WIN64) || defined (_WIN32)
// disable "This function may be unsafe" warnings for older C API functions
#pragma warning( disable: 4996 ) // instead of _CRT_SECURE_NO_WARNINGS, which doesn't work
#endif /* (_WIN64) || (_WIN32) */

void xerbla(const char* name, const int* num, const int /*len*/) {
    try {
        // MKL native thread attach to JVM
        ThreadAutoAttach treadAttach;
        JNIEnv* env = treadAttach.getContext();
        if (!env) {
            // TODO __LOG_WARN ...
            return;
        }
        try {
            int code = *num;
            if (code == 1001) {
                throwJavaRuntimeException(env, "%s %s %s", "Intel MKL ERROR:",
                    "Incompatible optional parameters on entry to", name);
            } else if (code == 1000 || code == 1089) {
                throwJavaRuntimeException(env, "%s %s %s", "Intel MKL INTERNAL ERROR:",
                    "Insufficient workspace available in function", name);
            } else if (code < 0) {
                throwJavaRuntimeException(env, "%s %s %s", "Intel MKL INTERNAL ERROR:",
                    "detected in function", name);
            } else {
                int position = code - 1;
                throwJavaRuntimeException(env, "%s %s %d %s %s", "Intel MKL ERROR:",
                    "Parameter", position, "was incorrect on entry to", name);
            }
        } catch (const JException& ex) {
            throwJavaRuntimeException(env, "%s %s", "xerbla", ex.what());
        } catch (...) {
            throwJavaRuntimeException(env, "%s", "xerbla: caught unknown exception");
        }
    } catch (const JException& /*ex*/) {
        // TODO __LOG_WARN ...
    } catch (...) {
        // TODO __LOG_WARN ...
    }
}

#ifdef __cplusplus
}
#endif // #ifdef __cplusplus
