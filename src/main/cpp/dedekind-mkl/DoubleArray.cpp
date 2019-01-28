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

#include "DoubleArray.h"

#ifndef _JAVASOFT_JNI_H_
#include <jni.h>
#endif /* _JAVASOFT_JNI_H_ */

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */



DoubleArray::DoubleArray(JNIEnv* env, jdoubleArray jarray, int offset, jboolean critical)
    : ctx(static_cast<Context*>(env)), jarray(jarray), offset(offset), critical(critical == JNI_TRUE)
{
    if (jarray) {
        jboolean isCopy = JNI_FALSE;
        if (critical) {
            carray = static_cast<double*>(ctx->GetPrimitiveArrayCritical(jarray, &isCopy));
        } else {
            carray = ctx->GetDoubleArrayElements(jarray, &isCopy);
        }
        if (carray == NULL) {
            throw JException("double* result: NULL");
        }
    } else {
        throw JException("jdoubleArray argument: null");
    }
}

double* DoubleArray::ptr() {
    return carray + offset;
}

DoubleArray::~DoubleArray() {
    if (critical) {
        ctx->ReleasePrimitiveArrayCritical(jarray, carray, 0);
    } else {
        ctx->ReleaseDoubleArrayElements(jarray, carray, 0);
    }
}
