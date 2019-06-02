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

#ifndef _JAVASOFT_JNI_H_
#include <jni.h>
#endif /* _JAVASOFT_JNI_H_ */

#ifndef DEDEKIND_INCLUDED_
#include "dedekind-mkl.h"
#endif /* DEDEKIND_INCLUDED_ */

#ifndef INTARRAY_INCLUDED_
#include "IntArray.h"
#endif /* INTARRAY_INCLUDED_ */

#ifndef DOUBLEARRAY_INCLUDED_
#include "DoubleArray.h"
#endif /* DOUBLEARRAY_INCLUDED_ */

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */

#ifndef JEXCEPTIONUTILS_INCLUDED_
#include "JExceptionUtils.h"
#endif /* JEXCEPTIONUTILS_INCLUDED_ */

#ifndef _MKL_H_
#include <mkl.h>
#endif /* _MKL_H_ */


static const int NOT_REACHED = -10000;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgbcon_n
 * Signature: (IBIII[DII[IIDLorg/netlib/util/doubleW;[DI[IIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgbcon_1n(JNIEnv* env, jclass,
  jint order,
  jbyte norm,
  jint n,
  jint kl,
  jint ku,
  jdoubleArray ab,
  jint abOffset,
  jint ldab,
  jintArray ipiv,
  jint ipivOffset,
  jdouble anorm,
  jobject rcondDW,
  jdoubleArray work,
  jint workOffset,
  jintArray iwork,
  jint iworkOffset,
  jboolean useCrit) {
    try {
        DoubleArray aba = DoubleArray(env, ab, abOffset, useCrit);
        IntArray ipiva = IntArray(env, ipiv, ipivOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        double rcond = 0.0;
        int info = LAPACKE_dgbcon_work(order, norm, n, kl, ku, aba.ptr(), ldab,
            ipiva.ptr(), anorm, &rcond, worka.ptr(), iworka.ptr());
        setDoubleWValue(env, rcondDW, rcond);
        return info;
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgbcon_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgbcon_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgbsv_n
 * Signature: (IIIII[DII[II[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgbsv_1n(JNIEnv* env, jclass,
  jint order,
  jint n,
  jint kl,
  jint ku,
  jint nrhs,
  jdoubleArray ab,
  jint abOffset,
  jint ldab,
  jintArray ipiv,
  jint ipivOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray aba = DoubleArray(env, ab, abOffset, useCrit);
        IntArray ipiva = IntArray(env, ipiv, ipivOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dgbsv(order, n, kl, ku, nrhs, aba.ptr(), ldab,
            ipiva.ptr(), ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgbsv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgbsv_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgbtrf_n
 * Signature: (IIIII[DII[IIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgbtrf_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jint kl,
  jint ku,
  jdoubleArray ab,
  jint abOffset,
  jint ldab,
  jintArray ipiv,
  jint ipivOffset,
  jboolean useCrit) {
    try {
        DoubleArray aba = DoubleArray(env, ab, abOffset, useCrit);
        IntArray ipiva = IntArray(env, ipiv, ipivOffset, useCrit);

        return LAPACKE_dgbtrf(order, m, n, kl, ku, aba.ptr(), ldab, ipiva.ptr());
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgbtrf_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgbtrf_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgbtrs_n
 * Signature: (IBIIII[DII[II[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgbtrs_1n(JNIEnv* env, jclass,
  jint order,
  jbyte trans,
  jint n,
  jint kl,
  jint ku,
  jint nrhs,
  jdoubleArray ab,
  jint abOffset,
  jint ldab,
  jintArray ipiv,
  jint ipivOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray aba = DoubleArray(env, ab, abOffset, useCrit);
        IntArray ipiva = IntArray(env, ipiv, ipivOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dgbtrs(order, trans, n, kl, ku, nrhs, aba.ptr(),
            ldab, ipiva.ptr(), ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgbtrs_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgbtrs_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgecon_n
 * Signature: (IBI[DIIDLorg/netlib/util/doubleW;[DI[IIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgecon_1n(JNIEnv* env, jclass,
  jint order,
  jbyte norm,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdouble anorm,
  jobject rcondDW,
  jdoubleArray work,
  jint workOffset,
  jintArray iwork,
  jint iworkOffset,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        double rcond = 0.0;
        int info = LAPACKE_dgecon_work(order, norm, n, aa.ptr(), lda,
            anorm, &rcond, worka.ptr(), iworka.ptr());
        setDoubleWValue(env, rcondDW, rcond);
        return info;
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgecon_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgecon_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgeev_n
 * Signature: (IBBI[DII[DI[DI[DII[DII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgeev_1n(JNIEnv* env, jclass,
  jint order,
  jbyte jobvl,
  jbyte jobvr,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray wr,
  jint wrOffset,
  jdoubleArray wi,
  jint wiOffset,
  jdoubleArray vl,
  jint vlOffset,
  jint ldvl,
  jdoubleArray vr,
  jint vrOffset,
  jint ldvr,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray wra = DoubleArray(env, wr, wrOffset, useCrit);
        DoubleArray wia = DoubleArray(env, wi, wiOffset, useCrit);
        DoubleArray vla = DoubleArray(env, vl, vlOffset, useCrit);
        DoubleArray vra = DoubleArray(env, vr, vrOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dgeev_work(order, jobvl, jobvr, n, aa.ptr(), lda,
            wra.ptr(), wia.ptr(), vla.ptr(), ldvl, vra.ptr(), ldvr, worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgeev_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgeev_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgelqf_n
 * Signature: (III[DII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgelqf_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray tau,
  jint tauOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray taua = DoubleArray(env, tau, tauOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dgelqf_work(order, m, n, aa.ptr(), lda, taua.ptr(),
            worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgelqf_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgelqf_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgels_n
 * Signature: (IBIII[DII[DII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgels_1n(JNIEnv* env, jclass,
  jint order,
  jbyte trans,
  jint m,
  jint n,
  jint nrhs,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dgels_work(order, trans, m, n, nrhs, aa.ptr(),
            lda, ba.ptr(), ldb, worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgels_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgels_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgeqlf_n
 * Signature: (III[DII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgeqlf_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray tau,
  jint tauOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray taua = DoubleArray(env, tau, tauOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dgeqlf_work(order, m, n, aa.ptr(), lda,
            taua.ptr(), worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgeqlf_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgeqlf_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgeqp3_n
 * Signature: (III[DII[II[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgeqp3_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jintArray jpvt,
  jint jpvtOffset,
  jdoubleArray tau,
  jint tauOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        IntArray jpvta = IntArray(env, jpvt, jpvtOffset, useCrit);
        DoubleArray taua = DoubleArray(env, tau, tauOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dgeqp3_work(order, m, n, aa.ptr(), lda, jpvta.ptr(),
            taua.ptr(), worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgeqp3_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgeqp3_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgeqrf_n
 * Signature: (III[DII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgeqrf_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray tau,
  jint tauOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray taua = DoubleArray(env, tau, tauOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dgeqrf_work(order, m, n, aa.ptr(), lda, taua.ptr(),
            worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgeqrf_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgeqrf_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgerqf_n
 * Signature: (III[DII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgerqf_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray tau,
  jint tauOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray taua = DoubleArray(env, tau, tauOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dgerqf_work(order, m, n, aa.ptr(), lda, taua.ptr(),
            worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgerqf_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgerqf_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgesdd_n
 * Signature: (IBII[DII[DI[DII[DII[DII[IIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgesdd_1n(JNIEnv* env, jclass,
  jint order,
  jbyte jobz,
  jint m,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray s,
  jint sOffset,
  jdoubleArray u,
  jint uOffset,
  jint ldu,
  jdoubleArray vt,
  jint vtOffset,
  jint ldvt,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jintArray iwork,
  jint iworkOffset,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray sa = DoubleArray(env, s, sOffset, useCrit);
        DoubleArray ua = DoubleArray(env, u, uOffset, useCrit);
        DoubleArray vta = DoubleArray(env, vt, vtOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        return LAPACKE_dgesdd_work(order, jobz, m, n, aa.ptr(), lda, sa.ptr(),
            ua.ptr(), ldu, vta.ptr(), ldvt, worka.ptr(), lwork, iworka.ptr());
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgesdd_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgesdd_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgesv_n
 * Signature: (III[DII[II[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgesv_1n(JNIEnv* env, jclass,
  jint order,
  jint n,
  jint nrhs,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jintArray ipiv,
  jint ipivOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        IntArray ipiva = IntArray(env, ipiv, ipivOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dgesv(order, n, nrhs, aa.ptr(), lda, ipiva.ptr(), ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgesv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgesv_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgetrf_n
 * Signature: (III[DII[IIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgetrf_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jintArray ipiv,
  jint ipivOffset,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        IntArray ipiva = IntArray(env, ipiv, ipivOffset, useCrit);

        return LAPACKE_dgetrf(order, m, n, aa.ptr(), lda, ipiva.ptr());
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgetrf_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgetrf_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgetrs_n
 * Signature: (IBII[DII[II[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgetrs_1n(JNIEnv* env, jclass,
  jint order,
  jbyte trans,
  jint n,
  jint nrhs,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jintArray ipiv,
  jint ipivOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        IntArray ipiva = IntArray(env, ipiv, ipivOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dgetrs(order, trans, n, nrhs, aa.ptr(), lda, ipiva.ptr(),
            ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgetrs_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgetrs_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dgtsv_n
 * Signature: (III[DI[DI[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dgtsv_1n(JNIEnv* env, jclass,
  jint order,
  jint n,
  jint nrhs,
  jdoubleArray dl,
  jint dlOffset,
  jdoubleArray d,
  jint dOffset,
  jdoubleArray du,
  jint duOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray dla = DoubleArray(env, dl, dlOffset, useCrit);
        DoubleArray da = DoubleArray(env, d, dOffset, useCrit);
        DoubleArray dua = DoubleArray(env, du, duOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dgtsv(order, n, nrhs, dla.ptr(), da.ptr(), dua.ptr(),
            ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dgtsv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dgtsv_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dlaswp_n
 * Signature: (II[DIIII[IIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_lapack_LapackN_dlaswp_1n(JNIEnv* env, jclass,
  jint order,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jint k1,
  jint k2,
  jintArray ipiv,
  jint ipivOffset,
  jint incx,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        IntArray ipiva = IntArray(env, ipiv, ipivOffset, useCrit);

        LAPACKE_dlaswp(order, n, aa.ptr(), lda, k1, k2, ipiva.ptr(), incx);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dlaswp_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dlaswp_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dorglq_n
 * Signature: (IIII[DII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dorglq_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jint k,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray tau,
  jint tauOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray taua = DoubleArray(env, tau, tauOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dorglq_work(order, m, n, k, aa.ptr(), lda, taua.ptr(),
            worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dorglq_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dorglq_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dorgql_n
 * Signature: (IIII[DII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dorgql_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jint k,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray tau,
  jint tauOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray taua = DoubleArray(env, tau, tauOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dorgql_work(order, m, n, k, aa.ptr(), lda, taua.ptr(),
            worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dorgql_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dorgql_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dorgqr_n
 * Signature: (IIII[DII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dorgqr_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jint k,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray tau,
  jint tauOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray taua = DoubleArray(env, tau, tauOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dorgqr_work(order, m, n, k, aa.ptr(), lda, taua.ptr(),
            worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dorgqr_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dorgqr_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dorgrq_n
 * Signature: (IIII[DII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dorgrq_1n(JNIEnv* env, jclass,
  jint order,
  jint m,
  jint n,
  jint k,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray tau,
  jint tauOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray taua = DoubleArray(env, tau, tauOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dorgrq_work(order, m, n, k, aa.ptr(), lda, taua.ptr(),
            worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dorgrq_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dorgrq_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dormrz_n
 * Signature: (IBBIIII[DII[DI[DII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dormrz_1n(JNIEnv* env, jclass,
  jint order,
  jbyte side,
  jbyte trans,
  jint m,
  jint n,
  jint k,
  jint l,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray tau,
  jint tauOffset,
  jdoubleArray c,
  jint cOffset,
  jint ldc,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray taua = DoubleArray(env, tau, tauOffset, useCrit);
        DoubleArray ca = DoubleArray(env, c, cOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dormrz_work(order, side, trans, m, n, k, l, aa.ptr(),
            lda, taua.ptr(), ca.ptr(), ldc, worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dormrz_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dormrz_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dpbcon_n
 * Signature: (IBII[DIIDLorg/netlib/util/doubleW;[DI[IIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dpbcon_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jint kd,
  jdoubleArray ab,
  jint abOffset,
  jint ldab,
  jdouble anorm,
  jobject rcondDW,
  jdoubleArray work,
  jint workOffset,
  jintArray iwork,
  jint iworkOffset,
  jboolean useCrit) {
    try {
        DoubleArray aba = DoubleArray(env, ab, abOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        double rcond = 0.0;
        int info = LAPACKE_dpbcon_work(order, uplo, n, kd, aba.ptr(), ldab,
            anorm, &rcond, worka.ptr(), iworka.ptr());
        setDoubleWValue(env, rcondDW, rcond);
        return info;
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dpbcon_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dpbcon_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dpbsv_n
 * Signature: (IBIII[DII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dpbsv_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jint kd,
  jint nrhs,
  jdoubleArray ab,
  jint abOffset,
  jint ldab,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray aba = DoubleArray(env, ab, abOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dpbsv(order, uplo, n, kd, nrhs, aba.ptr(), ldab,
            ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dpbsv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dpbsv_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dpbtrf_n
 * Signature: (IBII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dpbtrf_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jint kd,
  jdoubleArray ab,
  jint abOffset,
  jint ldab,
  jboolean useCrit) {
    try {
        DoubleArray aba = DoubleArray(env, ab, abOffset, useCrit);

        return LAPACKE_dpbtrf(order, uplo, n, kd, aba.ptr(), ldab);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dpbtrf_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dpbtrf_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dpbtrs_n
 * Signature: (IBIII[DII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dpbtrs_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jint kd,
  jint nrhs,
  jdoubleArray ab,
  jint abOffset,
  jint ldab,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray aba = DoubleArray(env, ab, abOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dpbtrs(order, uplo, n, kd, nrhs, aba.ptr(), ldab,
            ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dpbtrs_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dpbtrs_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dpocon_n
 * Signature: (IBI[DIIDLorg/netlib/util/doubleW;[DI[IIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dpocon_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdouble anorm,
  jobject rcondDW,
  jdoubleArray work,
  jint workOffset,
  jintArray iwork,
  jint iworkOffset,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        double rcond = 0.0;
        int info = LAPACKE_dpocon_work(order, uplo, n, aa.ptr(), lda, anorm,
            &rcond, worka.ptr(), iworka.ptr());
        setDoubleWValue(env, rcondDW, rcond);
        return info;
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dpocon_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dpocon_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dposv_n
 * Signature: (IBII[DII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dposv_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jint nrhs,
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

        return LAPACKE_dposv(order, uplo, n, nrhs, aa.ptr(), lda, ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dposv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dposv_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dpotrf_n
 * Signature: (IBI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dpotrf_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);

        return LAPACKE_dpotrf(order, uplo, n, aa.ptr(), lda);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dpotrf_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dpotrf_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dpotrs_n
 * Signature: (IBII[DII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dpotrs_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jint nrhs,
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

        return LAPACKE_dpotrs(order, uplo, n, nrhs, aa.ptr(), lda,
            ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dpotrs_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dpotrs_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dppcon_n
 * Signature: (IBI[DIDLorg/netlib/util/doubleW;[DI[IIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dppcon_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jdoubleArray ap,
  jint apOffset,
  jdouble anorm,
  jobject rcondDW,
  jdoubleArray work,
  jint workOffset,
  jintArray iwork,
  jint iworkOffset,
  jboolean useCrit) {
    try {
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        double rcond = 0.0;
        int info = LAPACKE_dppcon_work(order, uplo, n, apa.ptr(), anorm,
            &rcond, worka.ptr(), iworka.ptr());
        setDoubleWValue(env, rcondDW, rcond);
        return info;
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dppcon_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dppcon_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dppsv_n
 * Signature: (IBII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dppsv_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jint nrhs,
  jdoubleArray ap,
  jint apOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dppsv(order, uplo, n, nrhs, apa.ptr(), ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dppsv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dppsv_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dpptrf_n
 * Signature: (IBI[DIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dpptrf_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jdoubleArray ap,
  jint apOffset,
  jboolean useCrit) {
    try {
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);

        return LAPACKE_dpptrf(order, uplo, n, apa.ptr());
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dpptrf_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dpptrf_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dpptrs_n
 * Signature: (IBII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dpptrs_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jint nrhs,
  jdoubleArray ap,
  jint apOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dpptrs(order, uplo, n, nrhs, apa.ptr(), ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dpptrs_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dpptrs_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dptsv_n
 * Signature: (III[DI[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dptsv_1n(JNIEnv* env, jclass,
  jint order,
  jint n,
  jint nrhs,
  jdoubleArray d,
  jint dOffset,
  jdoubleArray e,
  jint eOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray da = DoubleArray(env, d, dOffset, useCrit);
        DoubleArray ea = DoubleArray(env, e, eOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dptsv(order, n, nrhs, da.ptr(), ea.ptr(), ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dptsv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dptsv_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dsbevd_n
 * Signature: (IBBII[DII[DI[DII[DII[IIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dsbevd_1n(JNIEnv* env, jclass,
  jint order,
  jbyte jobz,
  jbyte uplo,
  jint n,
  jint kd,
  jdoubleArray ab,
  jint abOffset,
  jint ldab,
  jdoubleArray w,
  jint wOffset,
  jdoubleArray z,
  jint zOffset,
  jint ldz,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jintArray iwork,
  jint iworkOffset,
  jint liwork,
  jboolean useCrit) {
    try {
        DoubleArray aba = DoubleArray(env, ab, abOffset, useCrit);
        DoubleArray wa = DoubleArray(env, w, wOffset, useCrit);
        DoubleArray za = DoubleArray(env, z, zOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        return LAPACKE_dsbevd_work(order, jobz, uplo, n, kd, aba.ptr(),
            ldab, wa.ptr(), za.ptr(), ldz, worka.ptr(), lwork, iworka.ptr(), liwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsbevd_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsbevd_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dspevd_n
 * Signature: (IBBI[DI[DI[DII[DII[IIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dspevd_1n(JNIEnv* env, jclass,
  jint order,
  jbyte jobz,
  jbyte uplo,
  jint n,
  jdoubleArray ap,
  jint apOffset,
  jdoubleArray w,
  jint wOffset,
  jdoubleArray z,
  jint zOffset,
  jint ldz,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jintArray iwork,
  jint iworkOffset,
  jint liwork,
  jboolean useCrit) {
    try {
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);
        DoubleArray wa = DoubleArray(env, w, wOffset, useCrit);
        DoubleArray za = DoubleArray(env, z, zOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        return LAPACKE_dspevd_work(order, jobz, uplo, n, apa.ptr(), wa.ptr(),
            za.ptr(), ldz, worka.ptr(), lwork, iworka.ptr(), liwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dspevd_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dspevd_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dspsv_n
 * Signature: (IBII[DI[II[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dspsv_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jint nrhs,
  jdoubleArray ap,
  jint apOffset,
  jintArray ipiv,
  jint ipivOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);
        IntArray ipiva = IntArray(env, ipiv, ipivOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dspsv(order, uplo, n, nrhs, apa.ptr(), ipiva.ptr(),
            ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dspsv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dspsv_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dstevr_n
 * Signature: (IBBI[DI[DIDDIIDLorg/netlib/util/intW;[DI[DII[II[DII[IIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dstevr_1n(JNIEnv* env, jclass,
  jint order,
  jbyte jobz,
  jbyte range,
  jint n,
  jdoubleArray d,
  jint dOffset,
  jdoubleArray e,
  jint eOffset,
  jdouble vl,
  jdouble vu,
  jint il,
  jint iu,
  jdouble abstol,
  jobject mIW,
  jdoubleArray w,
  jint wOffset,
  jdoubleArray z,
  jint zOffset,
  jint ldz,
  jintArray isuppz,
  jint isuppzOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jintArray iwork,
  jint iworkOffset,
  jint liwork,
  jboolean useCrit) {
    try {
        DoubleArray da = DoubleArray(env, d, dOffset, useCrit);
        DoubleArray ea = DoubleArray(env, e, eOffset, useCrit);
        DoubleArray wa = DoubleArray(env, w, wOffset, useCrit);
        DoubleArray za = DoubleArray(env, z, zOffset, useCrit);
        IntArray isuppza = IntArray(env, isuppz, isuppzOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        int m = 0;
        int info = LAPACKE_dstevr_work(order, jobz, range, n, da.ptr(),
            ea.ptr(), vl, vu, il, iu, abstol, &m, wa.ptr(), za.ptr(),
            ldz, isuppza.ptr(), worka.ptr(), lwork, iworka.ptr(), liwork);
        setIntWValue(env, mIW, m);
        return info;
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dstevr_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dstevr_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dsyevr_n
 * Signature: (IBBBI[DIIDDIIDLorg/netlib/util/intW;[DI[DII[II[DII[IIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dsyevr_1n(JNIEnv* env, jclass,
  jint order,
  jbyte jobz,
  jbyte range,
  jbyte uplo,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdouble vl,
  jdouble vu,
  jint il,
  jint iu,
  jdouble abstol,
  jobject mIW,
  jdoubleArray w,
  jint wOffset,
  jdoubleArray z,
  jint zOffset,
  jint ldz,
  jintArray isuppz,
  jint isuppzOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jintArray iwork,
  jint iworkOffset,
  jint liwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray wa = DoubleArray(env, w, wOffset, useCrit);
        DoubleArray za = DoubleArray(env, z, zOffset, useCrit);
        IntArray isuppza = IntArray(env, isuppz, isuppzOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        int m = 0;
        int info = LAPACKE_dsyevr_work(order, jobz, range, uplo, n, aa.ptr(),
            lda, vl, vu, il, iu, abstol, &m, wa.ptr(), za.ptr(), ldz,
            isuppza.ptr(), worka.ptr(), lwork, iworka.ptr(), liwork);
        setIntWValue(env, mIW, m);
        return info;
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsyevr_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsyevr_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dsygvd_n
 * Signature: (IIBBI[DII[DII[DI[DII[IIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dsygvd_1n(JNIEnv* env, jclass,
  jint order,
  jint itype,
  jbyte jobz,
  jbyte uplo,
  jint n,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jdoubleArray w,
  jint wOffset,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jintArray iwork,
  jint iworkOffset,
  jint liwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);
        DoubleArray wa = DoubleArray(env, w, wOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);
        IntArray iworka = IntArray(env, iwork, iworkOffset, useCrit);

        return LAPACKE_dsygvd_work(order, itype, jobz, uplo, n, aa.ptr(),
            lda, ba.ptr(), ldb, wa.ptr(), worka.ptr(), lwork, iworka.ptr(), liwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsygvd_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsygvd_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dsysv_n
 * Signature: (IBII[DII[II[DII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dsysv_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jint n,
  jint nrhs,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jintArray ipiv,
  jint ipivOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jdoubleArray work,
  jint workOffset,
  jint lwork,
  jboolean useCrit) {
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        IntArray ipiva = IntArray(env, ipiv, ipivOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);
        DoubleArray worka = DoubleArray(env, work, workOffset, useCrit);

        return LAPACKE_dsysv_work(order, uplo, n, nrhs, aa.ptr(), lda,
            ipiva.ptr(), ba.ptr(), ldb, worka.ptr(), lwork);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dsysv_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dsysv_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dtbtrs_n
 * Signature: (IBBBIII[DII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dtbtrs_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jbyte trans,
  jbyte diag,
  jint n,
  jint kd,
  jint nrhs,
  jdoubleArray ab,
  jint abOffset,
  jint ldab,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray aba = DoubleArray(env, ab, abOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dtbtrs(order, uplo, trans, diag, n, kd, nrhs,
            aba.ptr(), ldab, ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dtbtrs_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dtbtrs_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dtptrs_n
 * Signature: (IBBBII[DI[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dtptrs_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jbyte trans,
  jbyte diag,
  jint n,
  jint nrhs,
  jdoubleArray ap,
  jint apOffset,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit) {
    try {
        DoubleArray apa = DoubleArray(env, ap, apOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dtptrs(order, uplo, trans, diag, n, nrhs,
            apa.ptr(), ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dtptrs_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dtptrs_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
 * Class:     net_dedekind_lapack_LapackN
 * Method:    dtrtrs_n
 * Signature: (IBBBII[DII[DIIZ)I
 */
JNIEXPORT jint JNICALL
Java_net_dedekind_lapack_LapackN_dtrtrs_1n(JNIEnv* env, jclass,
  jint order,
  jbyte uplo,
  jbyte trans,
  jbyte diag,
  jint n,
  jint nrhs,
  jdoubleArray a,
  jint aOffset,
  jint lda,
  jdoubleArray b,
  jint bOffset,
  jint ldb,
  jboolean useCrit){
    try {
        DoubleArray aa = DoubleArray(env, a, aOffset, useCrit);
        DoubleArray ba = DoubleArray(env, b, bOffset, useCrit);

        return LAPACKE_dtrtrs(order, uplo, trans, diag, n, nrhs, aa.ptr(),
            lda, ba.ptr(), ldb);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dtrtrs_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dtrtrs_n: caught unknown exception");
    }
    return NOT_REACHED;
}

/*
* Class:     net_dedekind_lapack_LapackN
* Method:    initialize_n
* Signature: ()V
*/
JNIEXPORT void JNICALL
Java_net_dedekind_lapack_LapackN_initialize_1n(JNIEnv* env, jclass) {
    try {
        initializeDescriptors(env);
    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "initialize_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "initialize_n: caught unknown exception");
    }
}


#ifdef __cplusplus
}
#endif // #ifdef __cplusplus