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

#ifndef _JAVASOFT_JNI_H_
#include <jni.h>
#endif /* _JAVASOFT_JNI_H_ */

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

#ifndef _MKL_H_
#include <mkl.h>
#endif /* _MKL_H_ */


#ifdef __cplusplus
extern "C" {
#endif

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    cgemm3m_n
 * Signature: (IIIIIIFF[FI[FIFF[FIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_cgemm3m_1n(JNIEnv* env, jclass,
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

            cblas_cgemm3m(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_TRANSPOSE>(transa),
                static_cast<CBLAS_TRANSPOSE>(transb), m, n, k, &alpha, pa, lda, pb, ldb, &beta, pc, ldc);

            long len = cac.complexLength();
            if (len > 0 && cac.hasCopy()) {
                floatCopy(len, ca.ptr(), pc);
            }
        }

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "cgemm3m_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "cgemm3m_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    zgemm3m_n
 * Signature: (IIIIIIDD[DI[DIDD[DIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_zgemm3m_1n(JNIEnv* env, jclass,
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

            cblas_zgemm3m(static_cast<CBLAS_LAYOUT>(order), static_cast<CBLAS_TRANSPOSE>(transa),
                static_cast<CBLAS_TRANSPOSE>(transb), m, n, k, &alpha, pa, lda, pb, ldb, &beta, pc, ldc);

            long len = cac.complexLength();
            if (len > 0 && cac.hasCopy()) {
                doubleCopy(len, ca.ptr(), pc);
            }
        }

    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "zgemm3m_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "zgemm3m_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    simatcopy_n
 * Signature: (BBIIF[FIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_simatcopy_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte trans,
  jint rows,
  jint cols,
  jfloat alpha,
  jfloatArray AB,
  jint lda,
  jint ldb,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "simatcopy_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "simatcopy_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    dimatcopy_n
 * Signature: (BBIID[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_dimatcopy_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte trans,
  jint rows,
  jint cols,
  jdouble alpha,
  jdoubleArray AB,
  jint lda,
  jint ldb,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "dimatcopy_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "dimatcopy_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    cimatcopy_n
 * Signature: (BBIIFF[FIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_cimatcopy_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte trans,
  jint rows,
  jint cols,
  jfloat alphar,
  jfloat alphai,
  jfloatArray AB,
  jint lda,
  jint ldb,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "cimatcopy_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "cimatcopy_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    zimatcopy_n
 * Signature: (BBIIDD[DIIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_zimatcopy_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte trans,
  jint rows,
  jint cols,
  jdouble alphar,
  jdouble alphai,
  jdoubleArray AB,
  jint lda,
  jint ldb,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "zimatcopy_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "zimatcopy_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    somatcopy_n
 * Signature: (BBIIF[FI[FIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_somatcopy_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte trans,
  jint rows,
  jint cols,
  jfloat alpha,
  jfloatArray A,
  jint lda,
  jfloatArray B,
  jint ldb,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "somatcopy_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "somatcopy_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    domatcopy_n
 * Signature: (BBIID[DI[DIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_domatcopy_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte trans,
  jint rows,
  jint cols,
  jdouble alpha,
  jdoubleArray A,
  jint lda,
  jdoubleArray B,
  jint ldb,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "domatcopy_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "domatcopy_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    comatcopy_n
 * Signature: (BBIIFF[FI[FIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_comatcopy_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte trans,
  jint rows,
  jint cols,
  jfloat alphar,
  jfloat alphai,
  jfloatArray A,
  jint lda,
  jfloatArray B,
  jint ldb,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "comatcopy_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "comatcopy_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    zomatcopy_n
 * Signature: (BBIIDD[DI[DIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_zomatcopy_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte trans,
  jint rows,
  jint cols,
  jdouble alphar,
  jdouble alphai,
  jdoubleArray A,
  jint lda,
  jdoubleArray B,
  jint ldb,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "zomatcopy_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "zomatcopy_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    somatadd_n
 * Signature: (BBBIIF[FIF[FI[FIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_somatadd_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte transa,
  jbyte transb,
  jint m,
  jint n,
  jfloat alpha,
  jfloatArray A,
  jint lda,
  jfloat beta,
  jfloatArray B,
  jint ldb,
  jfloatArray C,
  jint ldc,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "somatadd_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "somatadd_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    domatadd_n
 * Signature: (BBBIID[DID[DI[DIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_domatadd_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte transa,
  jbyte transb,
  jint m,
  jint n,
  jdouble alpha,
  jdoubleArray A,
  jint lda,
  jdouble beta,
  jdoubleArray B,
  jint ldb,
  jdoubleArray C,
  jint ldc,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "domatadd_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "domatadd_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    comatadd_n
 * Signature: (BBBIIFF[FIFF[FI[FIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_comatadd_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte transa,
  jbyte transb,
  jint m,
  jint n,
  jfloat alphar,
  jfloat alphai,
  jfloatArray A,
  jint lda,
  jfloat betar,
  jfloat betai,
  jfloatArray B,
  jint ldb,
  jfloatArray C,
  jint ldc,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "comatadd_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "comatadd_n: caught unknown exception");
    }
}

/*
 * Class:     net_dedekind_blas_BlasExt
 * Method:    zomatadd_n
 * Signature: (BBBIIDD[DIDD[DI[DIZ)V
 */
JNIEXPORT void JNICALL
Java_net_dedekind_blas_BlasExt_zomatadd_1n(JNIEnv* env, jclass,
  jbyte ordering,
  jbyte transa,
  jbyte transb,
  jint m,
  jint n,
  jdouble alphar,
  jdouble alphai,
  jdoubleArray A,
  jint lda,
  jdouble betar,
  jdouble betai,
  jdoubleArray B,
  jint ldb,
  jdoubleArray C,
  jint ldc,
  jboolean useCrit) {
    try {


    } catch (const JException& ex) {
        throwJavaRuntimeException(env, "%s %s", "zomatadd_n", ex.what());
    } catch (...) {
        throwJavaRuntimeException(env, "%s", "zomatadd_n: caught unknown exception");
    }
}

#ifdef __cplusplus
}
#endif // #ifdef __cplusplus
