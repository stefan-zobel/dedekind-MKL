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

#ifndef STDAFX_INCLUDED_
#include "stdafx.h"
#endif /* STDAFX_INCLUDED_ */

#ifndef DEDEKIND_INCLUDED_
#include "dedekind-mkl.h"
#endif /* DEDEKIND_INCLUDED_ */

#ifndef JVMPROVIDER_INCLUDED_
#include "JvmProvider.h"
#endif /* JVMPROVIDER_INCLUDED_ */

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */


static jclass intWclass = NULL;
static jclass doubleWclass = NULL;

static jfieldID intWval = NULL;
static jfieldID doubleWval = NULL;


/**
* JNI_OnLoad
*/
JNIEXPORT jint JNICALL JNI_OnLoad(JavaVM* vm, void* /*reserved*/) {
    JvmProvider::instance()->initializeJavaVM(vm);
    return JNI_VERSION_1_6;
}

/**
* JNI_OnUnload
*/
JNIEXPORT void JNICALL JNI_OnUnload(JavaVM* /*vm*/, void* /*reserved*/) {
    JvmProvider::instance()->clearOnJavaVMUnload();
}


void setDoubleWValue(JNIEnv* env, jobject doubleW, double value) {
    if (doubleW) {
        if (doubleWclass != NULL && doubleWval != NULL) {
            Context* ctx = static_cast<Context*>(env);
            ctx->SetDoubleField(doubleW, doubleWval, value);
        } else {
            throw JException("doubleW descriptors not initialized");
        }
    } else {
        throw JException("parameter doubleW == null");
    }
}


void setIntWValue(JNIEnv* env, jobject intW, int value) {
    if (intW) {
        if (intWclass != NULL && intWval != NULL) {
            Context* ctx = static_cast<Context*>(env);
            ctx->SetIntField(intW, intWval, value);
        } else {
            throw JException("intW descriptors not initialized");
        }
    } else {
        throw JException("parameter intW == null");
    }
}


jclass findClass(Context* pCtx, const char* className) {
    jclass localClass = pCtx->FindClass(className);
    if (localClass == NULL) {
        SlimString msg("Cache init findClass - ");
        msg.append(className).append(" class not found");
        throw JException(msg);
    }
    jclass globalClass = reinterpret_cast<jclass>(pCtx->NewGlobalRef(localClass));
    if (globalClass == NULL) {
        pCtx->DeleteLocalRef(localClass);
        SlimString msg("Cache init findClass - ");
        msg.append(className).append(" unable to create global reference (out of memory?)");
        throw JException(msg);
    }
    return globalClass;
}


jfieldID findIntFieldId(Context* pCtx, jclass clazz, const char* fieldName) {
    jfieldID fieldId = pCtx->GetFieldID(clazz, fieldName, "I");
    if (fieldId == NULL) {
        SlimString msg("Cache init findIntFieldId - ");
        msg.append(fieldName).append(" int member field not found");
        throw JException(msg);
    }
    return fieldId;
}


jfieldID findDoubleFieldId(Context* pCtx, jclass clazz, const char* fieldName) {
    jfieldID fieldId = pCtx->GetFieldID(clazz, fieldName, "D");
    if (fieldId == NULL) {
        SlimString msg("Cache init findDoubleFieldId - ");
        msg.append(fieldName).append(" double member field not found");
        throw JException(msg);
    }
    return fieldId;
}


void initializeDescriptors(JNIEnv* env) {
    Context* ctx = static_cast<Context*>(env);
    intWclass = findClass(ctx, "org/netlib/util/intW");
    doubleWclass = findClass(ctx, "org/netlib/util/doubleW");
    intWval = findIntFieldId(ctx, intWclass, "val");
    doubleWval = findDoubleFieldId(ctx, doubleWclass, "val");
}
