/*
 * Copyright (c) 2002 - 2011    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

//////////////////////////////////////////////////////////////////////
// Context.cpp: implementation of the Context class.
//
//////////////////////////////////////////////////////////////////////

#include "Context.h"

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */

#ifndef JEXCEPTIONUTILS_INCLUDED_
#include "JExceptionUtils.h"
#endif /* JEXCEPTIONUTILS_INCLUDED_ */


#if defined (__GNUG__)
// disable warning "deprecated conversion from string constant to 'char*'" on G++ compiler
#pragma GCC diagnostic ignored "-Wwrite-strings"
#endif /* __GNUG__ */



//////////////////////////////////////////////////////////////////////
// Private helper function for determining the Java stacktrace
//////////////////////////////////////////////////////////////////////
static void throwJException(JNIEnv* env, jthrowable error, char* contextMethod) {
    SlimString stackTrace;
    printStackTrace(env, error, stackTrace, contextMethod);
    if (env->ExceptionCheck() == JNI_TRUE) {
        env->ExceptionClear();
    }
    throw JException(stackTrace);
}

void Context::clearException() {
    if (JNIEnv_::ExceptionCheck() == JNI_TRUE) {
        JNIEnv_::ExceptionClear();
    }
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Context::Context() {
}


Context::~Context() {
}

//////////////////////////////////////////////////////////////////////
// Member methods
//////////////////////////////////////////////////////////////////////

jint Context::GetVersion()
{
    clearException();
    jint result = JNIEnv_::GetVersion();

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetVersion");
    }

    return result;
}


jclass Context::DefineClass(const char* name, jobject loader, const jbyte* buf, jsize len)
{
    clearException();
    jclass result = JNIEnv_::DefineClass(name, loader, buf, len);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::DefineClass");
    }

    return result;
}


jclass Context::FindClass(const char* name)
{
    clearException();
    jclass result = JNIEnv_::FindClass(name);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::FindClass");
    }

    return result;
}


jmethodID Context::FromReflectedMethod(jobject method)
{
    clearException();
    jmethodID result = JNIEnv_::FromReflectedMethod(method);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::FromReflectedMethod");
    }

    return result;
}


jfieldID Context::FromReflectedField(jobject field)
{
    clearException();
    jfieldID result = JNIEnv_::FromReflectedField(field);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::FromReflectedField");
    }

    return result;
}


jobject Context::ToReflectedMethod(jclass clazz, jmethodID methodID, jboolean isStatic)
{
    clearException();
    jobject result = JNIEnv_::ToReflectedMethod(clazz, methodID, isStatic);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::ToReflectedMethod");
    }

    return result;
}


jclass Context::GetSuperclass(jclass sub)
{
    clearException();
    jclass result = JNIEnv_::GetSuperclass(sub);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetSuperclass");
    }

    return result;
}


jboolean Context::IsAssignableFrom(jclass sub, jclass sup)
{
    clearException();
    jboolean result = JNIEnv_::IsAssignableFrom(sub, sup);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::IsAssignableFrom");
    }

    return result;
}


jobject Context::ToReflectedField(jclass clazz, jfieldID fieldID, jboolean isStatic)
{
    clearException();
    jobject result = JNIEnv_::ToReflectedField(clazz, fieldID, isStatic);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::ToReflectedField");
    }

    return result;
}


jint Context::Throw(jthrowable obj)
{
    // exception functions should never throw
    return JNIEnv_::Throw(obj);
}


jint Context::ThrowNew(jclass clazz, const char* msg)
{
    // exception functions should never throw
    return JNIEnv_::ThrowNew(clazz, msg);
}


jthrowable Context::ExceptionOccurred()
{
    // exception functions should never throw
    return JNIEnv_::ExceptionOccurred();
}


void Context::ExceptionDescribe()
{
    // exception functions should never throw
    JNIEnv_::ExceptionDescribe();
}


void Context::ExceptionClear()
{
    // exception functions should never throw
    JNIEnv_::ExceptionClear();
}


void Context::FatalError(const char* msg)
{
    // exception functions should never throw
    JNIEnv_::FatalError(msg);
}


jint Context::PushLocalFrame(jint capacity)
{
    clearException();
    jint result = JNIEnv_::PushLocalFrame(capacity);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::PushLocalFrame");
    }

    return result;
}


jobject Context::PopLocalFrame(jobject res)
{
    clearException();
    jobject result = JNIEnv_::PopLocalFrame(res);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::PopLocalFrame");
    }

    return result;
}


jobject Context::NewGlobalRef(jobject obj)
{
    clearException();
    jobject result = JNIEnv_::NewGlobalRef(obj);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewGlobalRef");
    }

    return result;
}


void Context::DeleteGlobalRef(jobject gref)
{
    // deallocation functions should never throw
    JNIEnv_::DeleteGlobalRef(gref);
}


void Context::DeleteLocalRef(jobject obj)
{
    // deallocation functions should never throw
    JNIEnv_::DeleteLocalRef(obj);
}


jboolean Context::IsSameObject(jobject obj1, jobject obj2)
{
    clearException();
    jboolean result = JNIEnv_::IsSameObject(obj1, obj2);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::IsSameObject");
    }

    return result;
}


jobject Context::NewLocalRef(jobject ref)
{
    clearException();
    jobject result = JNIEnv_::NewLocalRef(ref);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewLocalRef");
    }

    return result;
}


jint Context::EnsureLocalCapacity(jint capacity)
{
    clearException();
    jint result = JNIEnv_::EnsureLocalCapacity(capacity);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::EnsureLocalCapacity");
    }

    return result;
}


jobject Context::AllocObject(jclass clazz)
{
    clearException();
    jobject result = JNIEnv_::AllocObject(clazz);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::AllocObject");
    }

    return result;
}


jobject Context::NewObject(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jobject result = JNIEnv_::NewObjectV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewObject");
    }

    return result;
}


jobject Context::NewObjectV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jobject result = JNIEnv_::NewObjectV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewObjectV");
    }

    return result;
}


jobject Context::NewObjectA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jobject result = JNIEnv_::NewObjectA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewObjectA");
    }

    return result;
}


jclass Context::GetObjectClass(jobject obj)
{
    clearException();
    jclass result = JNIEnv_::GetObjectClass(obj);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetObjectClass");
    }

    return result;
}


jboolean Context::IsInstanceOf(jobject obj, jclass clazz)
{
    clearException();
    jboolean result = JNIEnv_::IsInstanceOf(obj, clazz);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::IsInstanceOf");
    }

    return result;
}


jmethodID Context::GetMethodID(jclass clazz, const char* name, const char* sig)
{
    clearException();
    jmethodID result = JNIEnv_::GetMethodID(clazz, name, sig);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetMethodID");
    }

    return result;
}


jobject Context::CallObjectMethod(jobject obj, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jobject result = JNIEnv_::CallObjectMethodV(obj, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallObjectMethod");
    }

    return result;
}


jobject Context::CallObjectMethodV(jobject obj, jmethodID methodID, va_list args)
{
    clearException();
    jobject result = JNIEnv_::CallObjectMethodV(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallObjectMethodV");
    }

    return result;
}


jobject Context::CallObjectMethodA(jobject obj, jmethodID methodID, jvalue* args)
{
    clearException();
    jobject result = JNIEnv_::CallObjectMethodA(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallObjectMethodA");
    }

    return result;
}


jboolean Context::CallBooleanMethod(jobject obj, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jboolean result = JNIEnv_::CallBooleanMethodV(obj, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallBooleanMethod");
    }

    return result;
}


jboolean Context::CallBooleanMethodV(jobject obj, jmethodID methodID, va_list args)
{
    clearException();
    jboolean result = JNIEnv_::CallBooleanMethodV(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallBooleanMethodV");
    }

    return result;
}


jboolean Context::CallBooleanMethodA(jobject obj, jmethodID methodID, jvalue* args)
{
    clearException();
    jboolean result = JNIEnv_::CallBooleanMethodA(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallBooleanMethodA");
    }

    return result;
}


jbyte Context::CallByteMethod(jobject obj, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jbyte result = JNIEnv_::CallByteMethodV(obj, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallByteMethod");
    }

    return result;
}


jbyte Context::CallByteMethodV(jobject obj, jmethodID methodID, va_list args)
{
    clearException();
    jbyte result = JNIEnv_::CallByteMethodV(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallByteMethodV");
    }

    return result;
}


jbyte Context::CallByteMethodA(jobject obj, jmethodID methodID, jvalue* args)
{
    clearException();
    jbyte result = JNIEnv_::CallByteMethodA(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallByteMethodA");
    }

    return result;
}


jchar Context::CallCharMethod(jobject obj, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jchar result = JNIEnv_::CallCharMethodV(obj, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallCharMethod");
    }

    return result;
}


jchar Context::CallCharMethodV(jobject obj, jmethodID methodID, va_list args)
{
    clearException();
    jchar result = JNIEnv_::CallCharMethodV(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallCharMethodV");
    }

    return result;
}


jchar Context::CallCharMethodA(jobject obj, jmethodID methodID, jvalue* args)
{
    clearException();
    jchar result = JNIEnv_::CallCharMethodA(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallCharMethodA");
    }

    return result;
}


jshort Context::CallShortMethod(jobject obj, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jshort result = JNIEnv_::CallShortMethodV(obj, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallShortMethod");
    }

    return result;
}


jshort Context::CallShortMethodV(jobject obj, jmethodID methodID, va_list args)
{
    clearException();
    jshort result = JNIEnv_::CallShortMethodV(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallShortMethodV");
    }

    return result;
}


jshort Context::CallShortMethodA(jobject obj, jmethodID methodID, jvalue* args)
{
    clearException();
    jshort result = JNIEnv_::CallShortMethodA(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallShortMethodA");
    }

    return result;
}


jint Context::CallIntMethod(jobject obj, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jint result = JNIEnv_::CallIntMethodV(obj, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallIntMethod");
    }

    return result;
}


jint Context::CallIntMethodV(jobject obj, jmethodID methodID, va_list args)
{
    clearException();
    jint result = JNIEnv_::CallIntMethodV(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallIntMethodV");
    }

    return result;
}


jint Context::CallIntMethodA(jobject obj, jmethodID methodID, jvalue* args)
{
    clearException();
    jint result = JNIEnv_::CallIntMethodA(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallIntMethodA");
    }

    return result;
}


jlong Context::CallLongMethod(jobject obj, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jlong result = JNIEnv_::CallLongMethodV(obj, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallLongMethod");
    }

    return result;
}


jlong Context::CallLongMethodV(jobject obj, jmethodID methodID, va_list args)
{
    clearException();
    jlong result = JNIEnv_::CallLongMethodV(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallLongMethodV");
    }

    return result;
}


jlong Context::CallLongMethodA(jobject obj, jmethodID methodID, jvalue* args)
{
    clearException();
    jlong result = JNIEnv_::CallLongMethodA(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallLongMethodA");
    }

    return result;
}


jfloat Context::CallFloatMethod(jobject obj, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jfloat result = JNIEnv_::CallFloatMethodV(obj, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallFloatMethod");
    }

    return result;
}


jfloat Context::CallFloatMethodV(jobject obj, jmethodID methodID, va_list args)
{
    clearException();
    jfloat result = JNIEnv_::CallFloatMethodV(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallFloatMethodV");
    }

    return result;
}


jfloat Context::CallFloatMethodA(jobject obj, jmethodID methodID, jvalue* args)
{
    clearException();
    jfloat result = JNIEnv_::CallFloatMethodA(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallFloatMethodA");
    }

    return result;
}


jdouble Context::CallDoubleMethod(jobject obj, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jdouble result = JNIEnv_::CallDoubleMethodV(obj, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallDoubleMethod");
    }

    return result;
}


jdouble Context::CallDoubleMethodV(jobject obj, jmethodID methodID, va_list args)
{
    clearException();
    jdouble result = JNIEnv_::CallDoubleMethodV(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallDoubleMethodV");
    }

    return result;
}


jdouble Context::CallDoubleMethodA(jobject obj, jmethodID methodID, jvalue* args)
{
    clearException();
    jdouble result = JNIEnv_::CallDoubleMethodA(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallDoubleMethodA");
    }

    return result;
}


void Context::CallVoidMethod(jobject obj, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    JNIEnv_::CallVoidMethodV(obj, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallVoidMethod");
    }
}


void Context::CallVoidMethodV(jobject obj, jmethodID methodID, va_list args)
{
    clearException();
    JNIEnv_::CallVoidMethodV(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallVoidMethodV");
    }
}


void Context::CallVoidMethodA(jobject obj, jmethodID methodID, jvalue* args)
{
    clearException();
    JNIEnv_::CallVoidMethodA(obj, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallVoidMethodA");
    }
}


jobject Context::CallNonvirtualObjectMethod(jobject obj, jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jobject result = JNIEnv_::CallNonvirtualObjectMethodV(obj, clazz,   methodID,args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualObjectMethod");
    }

    return result;
}


jobject Context::CallNonvirtualObjectMethodV(jobject obj, jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jobject result = JNIEnv_::CallNonvirtualObjectMethodV(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualObjectMethodV");
    }

    return result;
}


jobject Context::CallNonvirtualObjectMethodA(jobject obj, jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jobject result = JNIEnv_::CallNonvirtualObjectMethodA(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualObjectMethodA");
    }

    return result;
}


jboolean Context::CallNonvirtualBooleanMethod(jobject obj, jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jboolean result = JNIEnv_::CallNonvirtualBooleanMethodV(obj, clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualBooleanMethod");
    }

    return result;
}


jboolean Context::CallNonvirtualBooleanMethodV(jobject obj, jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jboolean result = JNIEnv_::CallNonvirtualBooleanMethodV(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualBooleanMethodV");
    }

    return result;
}


jboolean Context::CallNonvirtualBooleanMethodA(jobject obj, jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jboolean result = JNIEnv_::CallNonvirtualBooleanMethodA(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualBooleanMethodA");
    }

    return result;
}


jbyte Context::CallNonvirtualByteMethod(jobject obj, jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jbyte result = JNIEnv_::CallNonvirtualByteMethodV(obj, clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualByteMethod");
    }

    return result;
}


jbyte Context::CallNonvirtualByteMethodV(jobject obj, jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jbyte result = JNIEnv_::CallNonvirtualByteMethodV(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualByteMethodV");
    }

    return result;
}


jbyte Context::CallNonvirtualByteMethodA(jobject obj, jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jbyte result = JNIEnv_::CallNonvirtualByteMethodA(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualByteMethodA");
    }

    return result;
}


jchar Context::CallNonvirtualCharMethod(jobject obj, jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jchar result = JNIEnv_::CallNonvirtualCharMethodV(obj, clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualCharMethod");
    }

    return result;
}


jchar Context::CallNonvirtualCharMethodV(jobject obj, jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jchar result = JNIEnv_::CallNonvirtualCharMethodV(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualCharMethodV");
    }

    return result;
}


jchar Context::CallNonvirtualCharMethodA(jobject obj, jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jchar result = JNIEnv_::CallNonvirtualCharMethodA(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualCharMethodA");
    }

    return result;
}


jshort Context::CallNonvirtualShortMethod(jobject obj, jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jshort result = JNIEnv_::CallNonvirtualShortMethodV(obj, clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualShortMethod");
    }

    return result;
}


jshort Context::CallNonvirtualShortMethodV(jobject obj, jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jshort result = JNIEnv_::CallNonvirtualShortMethodV(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualShortMethodV");
    }

    return result;
}


jshort Context::CallNonvirtualShortMethodA(jobject obj, jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jshort result = JNIEnv_::CallNonvirtualShortMethodA(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualShortMethodA");
    }

    return result;
}


jint Context::CallNonvirtualIntMethod(jobject obj, jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jint result = JNIEnv_::CallNonvirtualIntMethodV(obj, clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualIntMethod");
    }

    return result;
}


jint Context::CallNonvirtualIntMethodV(jobject obj, jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jint result = JNIEnv_::CallNonvirtualIntMethodV(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualIntMethodV");
    }

    return result;
}


jint Context::CallNonvirtualIntMethodA(jobject obj, jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jint result = JNIEnv_::CallNonvirtualIntMethodA(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualIntMethodA");
    }

    return result;
}


jlong Context::CallNonvirtualLongMethod(jobject obj, jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jlong result = JNIEnv_::CallNonvirtualLongMethodV(obj, clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualLongMethod");
    }

    return result;
}


jlong Context::CallNonvirtualLongMethodV(jobject obj, jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jlong result = JNIEnv_::CallNonvirtualLongMethodV(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualLongMethodV");
    }

    return result;
}


jlong Context::CallNonvirtualLongMethodA(jobject obj, jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jlong result = JNIEnv_::CallNonvirtualLongMethodA(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualLongMethodA");
    }

    return result;
}


jfloat Context::CallNonvirtualFloatMethod(jobject obj, jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jfloat result = JNIEnv_::CallNonvirtualFloatMethodV(obj, clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualFloatMethod");
    }

    return result;
}


jfloat Context::CallNonvirtualFloatMethodV(jobject obj, jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jfloat result = JNIEnv_::CallNonvirtualFloatMethodV(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualFloatMethodV");
    }

    return result;
}


jfloat Context::CallNonvirtualFloatMethodA(jobject obj, jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jfloat result = JNIEnv_::CallNonvirtualFloatMethodA(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualFloatMethodA");
    }

    return result;
}


jdouble Context::CallNonvirtualDoubleMethod(jobject obj, jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jdouble result = JNIEnv_::CallNonvirtualDoubleMethodV(obj, clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualDoubleMethod");
    }

    return result;
}


jdouble Context::CallNonvirtualDoubleMethodV(jobject obj, jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jdouble result = JNIEnv_::CallNonvirtualDoubleMethodV(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualDoubleMethodV");
    }

    return result;
}


jdouble Context::CallNonvirtualDoubleMethodA(jobject obj, jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jdouble result = JNIEnv_::CallNonvirtualDoubleMethodA(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualDoubleMethodA");
    }

    return result;
}


void Context::CallNonvirtualVoidMethod(jobject obj, jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    JNIEnv_::CallNonvirtualVoidMethodV(obj, clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualVoidMethod");
    }
}


void Context::CallNonvirtualVoidMethodV(jobject obj, jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    JNIEnv_::CallNonvirtualVoidMethodV(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualVoidMethodV");
    }
}


void Context::CallNonvirtualVoidMethodA(jobject obj, jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    JNIEnv_::CallNonvirtualVoidMethodA(obj, clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallNonvirtualVoidMethodA");
    }
}


jfieldID Context::GetFieldID(jclass clazz, const char* name, const char* sig)
{
    clearException();
    jfieldID result = JNIEnv_::GetFieldID(clazz, name, sig);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetFieldID");
    }

    return result;
}


jobject Context::GetObjectField(jobject obj, jfieldID fieldID)
{
    clearException();
    jobject result = JNIEnv_::GetObjectField(obj, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetObjectField");
    }

    return result;
}


jboolean Context::GetBooleanField(jobject obj, jfieldID fieldID)
{
    clearException();
    jboolean result = JNIEnv_::GetBooleanField(obj, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetBooleanField");
    }

    return result;
}


jbyte Context::GetByteField(jobject obj, jfieldID fieldID)
{
    clearException();
    jbyte result = JNIEnv_::GetByteField(obj, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetByteField");
    }

    return result;
}


jchar Context::GetCharField(jobject obj, jfieldID fieldID)
{
    clearException();
    jchar result = JNIEnv_::GetCharField(obj, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetCharField");
    }

    return result;
}


jshort Context::GetShortField(jobject obj, jfieldID fieldID)
{
    clearException();
    jshort result = JNIEnv_::GetShortField(obj, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetShortField");
    }

    return result;
}


jint Context::GetIntField(jobject obj, jfieldID fieldID)
{
    clearException();
    jint result = JNIEnv_::GetIntField(obj, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetIntField");
    }

    return result;
}


jlong Context::GetLongField(jobject obj, jfieldID fieldID)
{
    clearException();
    jlong result = JNIEnv_::GetLongField(obj, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetLongField");
    }

    return result;
}


jfloat Context::GetFloatField(jobject obj, jfieldID fieldID)
{
    clearException();
    jfloat result = JNIEnv_::GetFloatField(obj, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetFloatField");
    }

    return result;
}


jdouble Context::GetDoubleField(jobject obj, jfieldID fieldID)
{
    clearException();
    jdouble result = JNIEnv_::GetDoubleField(obj, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetDoubleField");
    }

    return result;
}


void Context::SetObjectField(jobject obj, jfieldID fieldID, jobject val)
{
    clearException();
    JNIEnv_::SetObjectField(obj, fieldID, val);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetObjectField");
    }
}


void Context::SetBooleanField(jobject obj, jfieldID fieldID, jboolean val)
{
    clearException();
    JNIEnv_::SetBooleanField(obj, fieldID, val);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetBooleanField");
    }
}


void Context::SetByteField(jobject obj, jfieldID fieldID, jbyte val)
{
    clearException();
    JNIEnv_::SetByteField(obj, fieldID, val);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetByteField");
    }
}


void Context::SetCharField(jobject obj, jfieldID fieldID, jchar val)
{
    clearException();
    JNIEnv_::SetCharField(obj, fieldID, val);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetCharField");
    }
}


void Context::SetShortField(jobject obj, jfieldID fieldID, jshort val)
{
    clearException();
    JNIEnv_::SetShortField(obj, fieldID, val);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetShortField");
    }
}


void Context::SetIntField(jobject obj, jfieldID fieldID, jint val)
{
    clearException();
    JNIEnv_::SetIntField(obj, fieldID, val);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetIntField");
    }
}


void Context::SetLongField(jobject obj, jfieldID fieldID, jlong val)
{
    clearException();
    JNIEnv_::SetLongField(obj, fieldID, val);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetLongField");
    }
}


void Context::SetFloatField(jobject obj, jfieldID fieldID, jfloat val)
{
    clearException();
    JNIEnv_::SetFloatField(obj, fieldID, val);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetFloatField");
    }
}


void Context::SetDoubleField(jobject obj, jfieldID fieldID, jdouble val)
{
    clearException();
    JNIEnv_::SetDoubleField(obj, fieldID, val);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetDoubleField");
    }
}


jmethodID Context::GetStaticMethodID(jclass clazz, const char* name, const char* sig)
{
    clearException();
    jmethodID result = JNIEnv_::GetStaticMethodID(clazz, name, sig);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticMethodID");
    }

    return result;
}


jobject Context::CallStaticObjectMethod(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jobject result = JNIEnv_::CallStaticObjectMethodV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticObjectMethod");
    }

    return result;
}


jobject Context::CallStaticObjectMethodV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jobject result = JNIEnv_::CallStaticObjectMethodV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticObjectMethodV");
    }

    return result;
}


jobject Context::CallStaticObjectMethodA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jobject result = JNIEnv_::CallStaticObjectMethodA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticObjectMethodA");
    }

    return result;
}


jboolean Context::CallStaticBooleanMethod(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jboolean result = JNIEnv_::CallStaticBooleanMethodV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticBooleanMethod");
    }

    return result;
}


jboolean Context::CallStaticBooleanMethodV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jboolean result = JNIEnv_::CallStaticBooleanMethodV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticBooleanMethodV");
    }

    return result;
}


jboolean Context::CallStaticBooleanMethodA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jboolean result = JNIEnv_::CallStaticBooleanMethodA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticBooleanMethodA");
    }

    return result;
}


jbyte Context::CallStaticByteMethod(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jbyte result = JNIEnv_::CallStaticByteMethodV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticByteMethod");
    }

    return result;
}


jbyte Context::CallStaticByteMethodV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jbyte result = JNIEnv_::CallStaticByteMethodV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticByteMethodV");
    }

    return result;
}


jbyte Context::CallStaticByteMethodA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jbyte result = JNIEnv_::CallStaticByteMethodA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticByteMethodA");
    }

    return result;
}


jchar Context::CallStaticCharMethod(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jchar result = JNIEnv_::CallStaticCharMethodV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticCharMethod");
    }

    return result;
}


jchar Context::CallStaticCharMethodV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jchar result = JNIEnv_::CallStaticCharMethodV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticCharMethodV");
    }

    return result;
}


jchar Context::CallStaticCharMethodA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jchar result = JNIEnv_::CallStaticCharMethodA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticCharMethodA");
    }

    return result;
}


jshort Context::CallStaticShortMethod(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jshort result = JNIEnv_::CallStaticShortMethodV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticShortMethod");
    }

    return result;
}


jshort Context::CallStaticShortMethodV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jshort result = JNIEnv_::CallStaticShortMethodV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticShortMethodV");
    }

    return result;
}


jshort Context::CallStaticShortMethodA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jshort result = JNIEnv_::CallStaticShortMethodA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticShortMethodA");
    }

    return result;
}


jint Context::CallStaticIntMethod(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jint result = JNIEnv_::CallStaticIntMethodV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticIntMethod");
    }

    return result;
}


jint Context::CallStaticIntMethodV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jint result = JNIEnv_::CallStaticIntMethodV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticIntMethodV");
    }

    return result;
}


jint Context::CallStaticIntMethodA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jint result = JNIEnv_::CallStaticIntMethodA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticIntMethodA");
    }

    return result;
}


jlong Context::CallStaticLongMethod(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jlong result = JNIEnv_::CallStaticLongMethodV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticLongMethod");
    }

    return result;
}


jlong Context::CallStaticLongMethodV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jlong result = JNIEnv_::CallStaticLongMethodV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticLongMethodV");
    }

    return result;
}


jlong Context::CallStaticLongMethodA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jlong result = JNIEnv_::CallStaticLongMethodA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticLongMethodA");
    }

    return result;
}


jfloat Context::CallStaticFloatMethod(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jfloat result = JNIEnv_::CallStaticFloatMethodV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticFloatMethod");
    }

    return result;
}


jfloat Context::CallStaticFloatMethodV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jfloat result = JNIEnv_::CallStaticFloatMethodV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticFloatMethodV");
    }

    return result;
}


jfloat Context::CallStaticFloatMethodA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jfloat result = JNIEnv_::CallStaticFloatMethodA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticFloatMethodA");
    }

    return result;
}


jdouble Context::CallStaticDoubleMethod(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    jdouble result = JNIEnv_::CallStaticDoubleMethodV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticDoubleMethod");
    }

    return result;
}


jdouble Context::CallStaticDoubleMethodV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    jdouble result = JNIEnv_::CallStaticDoubleMethodV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticDoubleMethodV");
    }

    return result;
}


jdouble Context::CallStaticDoubleMethodA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    jdouble result = JNIEnv_::CallStaticDoubleMethodA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticDoubleMethodA");
    }

    return result;
}


void Context::CallStaticVoidMethod(jclass clazz, jmethodID methodID, ...)
{
    clearException();
    va_list args;
    va_start(args, methodID);

    JNIEnv_::CallStaticVoidMethodV(clazz, methodID, args);

    va_end(args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticVoidMethod");
    }
}


void Context::CallStaticVoidMethodV(jclass clazz, jmethodID methodID, va_list args)
{
    clearException();
    JNIEnv_::CallStaticVoidMethodV(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticVoidMethodV");
    }
}


void Context::CallStaticVoidMethodA(jclass clazz, jmethodID methodID, jvalue* args)
{
    clearException();
    JNIEnv_::CallStaticVoidMethodA(clazz, methodID, args);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::CallStaticVoidMethodA");
    }
}


jfieldID Context::GetStaticFieldID(jclass clazz, const char* name, const char* sig)
{
    clearException();
    jfieldID result = JNIEnv_::GetStaticFieldID(clazz, name, sig);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticFieldID");
    }

    return result;
}


jobject Context::GetStaticObjectField(jclass clazz, jfieldID fieldID)
{
    clearException();
    jobject result = JNIEnv_::GetStaticObjectField(clazz, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticObjectField");
    }

    return result;
}


jboolean Context::GetStaticBooleanField(jclass clazz, jfieldID fieldID)
{
    clearException();
    jboolean result = JNIEnv_::GetStaticBooleanField(clazz, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticBooleanField");
    }

    return result;
}


jbyte Context::GetStaticByteField(jclass clazz, jfieldID fieldID)
{
    clearException();
    jbyte result = JNIEnv_::GetStaticByteField(clazz, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticByteField");
    }

    return result;
}


jchar Context::GetStaticCharField(jclass clazz, jfieldID fieldID)
{
    clearException();
    jchar result = JNIEnv_::GetStaticCharField(clazz, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticCharField");
    }

    return result;
}


jshort Context::GetStaticShortField(jclass clazz, jfieldID fieldID)
{
    clearException();
    jshort result = JNIEnv_::GetStaticShortField(clazz, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticShortField");
    }

    return result;
}


jint Context::GetStaticIntField(jclass clazz, jfieldID fieldID)
{
    clearException();
    jint result = JNIEnv_::GetStaticIntField(clazz, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticIntField");
    }

    return result;
}


jlong Context::GetStaticLongField(jclass clazz, jfieldID fieldID)
{
    clearException();
    jlong result = JNIEnv_::GetStaticLongField(clazz, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticLongField");
    }

    return result;
}


jfloat Context::GetStaticFloatField(jclass clazz, jfieldID fieldID)
{
    clearException();
    jfloat result = JNIEnv_::GetStaticFloatField(clazz, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticFloatField");
    }

    return result;
}


jdouble Context::GetStaticDoubleField(jclass clazz, jfieldID fieldID)
{
    clearException();
    jdouble result = JNIEnv_::GetStaticDoubleField(clazz, fieldID);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStaticDoubleField");
    }

    return result;
}


void Context::SetStaticObjectField(jclass clazz, jfieldID fieldID, jobject value)
{
    clearException();
    JNIEnv_::SetStaticObjectField(clazz, fieldID, value);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetStaticObjectField");
    }
}


void Context::SetStaticBooleanField(jclass clazz, jfieldID fieldID, jboolean value)
{
    clearException();
    JNIEnv_::SetStaticBooleanField(clazz, fieldID, value);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetStaticBooleanField");
    }
}


void Context::SetStaticByteField(jclass clazz, jfieldID fieldID, jbyte value)
{
    clearException();
    JNIEnv_::SetStaticByteField(clazz, fieldID, value);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetStaticByteField");
    }
}


void Context::SetStaticCharField(jclass clazz, jfieldID fieldID, jchar value)
{
    clearException();
    JNIEnv_::SetStaticCharField(clazz, fieldID, value);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetStaticCharField");
    }
}


void Context::SetStaticShortField(jclass clazz, jfieldID fieldID, jshort value)
{
    clearException();
    JNIEnv_::SetStaticShortField(clazz, fieldID, value);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetStaticShortField");
    }
}


void Context::SetStaticIntField(jclass clazz, jfieldID fieldID, jint value)
{
    clearException();
    JNIEnv_::SetStaticIntField(clazz, fieldID, value);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetStaticIntField");
    }
}


void Context::SetStaticLongField(jclass clazz, jfieldID fieldID, jlong value)
{
    clearException();
    JNIEnv_::SetStaticLongField(clazz, fieldID, value);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetStaticLongField");
    }
}


void Context::SetStaticFloatField(jclass clazz, jfieldID fieldID, jfloat value)
{
    clearException();
    JNIEnv_::SetStaticFloatField(clazz, fieldID, value);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetStaticFloatField");
    }
}


void Context::SetStaticDoubleField(jclass clazz, jfieldID fieldID, jdouble value)
{
    clearException();
    JNIEnv_::SetStaticDoubleField(clazz, fieldID, value);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetStaticDoubleField");
    }
}


jstring Context::NewString(const jchar* unicode, jsize len)
{
    clearException();
    jstring result = JNIEnv_::NewString(unicode, len);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewString");
    }

    return result;
}


jsize Context::GetStringLength(jstring str)
{
    clearException();
    jsize result = JNIEnv_::GetStringLength(str);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStringLength");
    }

    return result;
}


const jchar* Context::GetStringChars(jstring str, jboolean* isCopy)
{
    clearException();
    const jchar* result = JNIEnv_::GetStringChars(str, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStringChars");
    }

    return result;
}


void Context::ReleaseStringChars(jstring str, const jchar* chars)
{
    // release funtions should never throw
    JNIEnv_::ReleaseStringChars(str, chars);
}


jstring Context::NewStringUTF(const char* utf)
{
    clearException();
    jstring result = JNIEnv_::NewStringUTF(utf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewStringUTF");
    }

    return result;
}


jsize Context::GetStringUTFLength(jstring str)
{
    clearException();
    jsize result = JNIEnv_::GetStringUTFLength(str);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStringUTFLength");
    }

    return result;
}


const char* Context::GetStringUTFChars(jstring str, jboolean* isCopy)
{
    clearException();
    const char* result = JNIEnv_::GetStringUTFChars(str, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStringUTFChars");
    }

    return result;
}


void Context::ReleaseStringUTFChars(jstring str, const char* chars)
{
    // release funtions should never throw
    JNIEnv_::ReleaseStringUTFChars(str, chars);
}


jsize Context::GetArrayLength(jarray array)
{
    clearException();
    jsize result = JNIEnv_::GetArrayLength(array);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetArrayLength");
    }

    return result;
}


jobjectArray Context::NewObjectArray(jsize len, jclass clazz, jobject init)
{
    clearException();
    jobjectArray result = JNIEnv_::NewObjectArray(len, clazz, init);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewObjectArray");
    }

    return result;
}


jobject Context::GetObjectArrayElement(jobjectArray array, jsize index)
{
    clearException();
    jobject result = JNIEnv_::GetObjectArrayElement(array, index);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetObjectArrayElement");
    }

    return result;
}


void Context::SetObjectArrayElement(jobjectArray array, jsize index, jobject val)
{
    clearException();
    JNIEnv_::SetObjectArrayElement(array, index, val);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetObjectArrayElement");
    }
}


jbooleanArray Context::NewBooleanArray(jsize len)
{
    clearException();
    jbooleanArray result = JNIEnv_::NewBooleanArray(len);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewBooleanArray");
    }

    return result;
}


jbyteArray Context::NewByteArray(jsize len)
{
    clearException();
    jbyteArray result = JNIEnv_::NewByteArray(len);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewByteArray");
    }

    return result;
}


jcharArray Context::NewCharArray(jsize len)
{
    clearException();
    jcharArray result = JNIEnv_::NewCharArray(len);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewCharArray");
    }

    return result;
}


jshortArray Context::NewShortArray(jsize len)
{
    clearException();
    jshortArray result = JNIEnv_::NewShortArray(len);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewShortArray");
    }

    return result;
}


jintArray Context::NewIntArray(jsize len)
{
    clearException();
    jintArray result = JNIEnv_::NewIntArray(len);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewIntArray");
    }

    return result;
}


jlongArray Context::NewLongArray(jsize len)
{
    clearException();
    jlongArray result = JNIEnv_::NewLongArray(len);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewLongArray");
    }

    return result;
}


jfloatArray Context::NewFloatArray(jsize len)
{
    clearException();
    jfloatArray result = JNIEnv_::NewFloatArray(len);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewFloatArray");
    }

    return result;
}


jdoubleArray Context::NewDoubleArray(jsize len)
{
    clearException();
    jdoubleArray result = JNIEnv_::NewDoubleArray(len);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewDoubleArray");
    }

    return result;
}


jboolean* Context::GetBooleanArrayElements(jbooleanArray array, jboolean* isCopy)
{
    clearException();
    jboolean* result = JNIEnv_::GetBooleanArrayElements(array, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetBooleanArrayElements");
    }

    return result;
}


jbyte* Context::GetByteArrayElements(jbyteArray array, jboolean* isCopy)
{
    clearException();
    jbyte* result = JNIEnv_::GetByteArrayElements(array, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetByteArrayElements");
    }

    return result;
}


jchar* Context::GetCharArrayElements(jcharArray array, jboolean* isCopy)
{
    clearException();
    jchar* result = JNIEnv_::GetCharArrayElements(array, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetCharArrayElements");
    }

    return result;
}


jshort* Context::GetShortArrayElements(jshortArray array, jboolean* isCopy)
{
    clearException();
    jshort* result = JNIEnv_::GetShortArrayElements(array, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetShortArrayElements");
    }

    return result;
}


jint* Context::GetIntArrayElements(jintArray array, jboolean* isCopy)
{
    clearException();
    jint* result = JNIEnv_::GetIntArrayElements(array, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetIntArrayElements");
    }

    return result;
}


jlong* Context::GetLongArrayElements(jlongArray array, jboolean* isCopy)
{
    clearException();
    jlong* result = JNIEnv_::GetLongArrayElements(array, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetLongArrayElements");
    }

    return result;
}


jfloat* Context::GetFloatArrayElements(jfloatArray array, jboolean* isCopy)
{
    clearException();
    jfloat* result = JNIEnv_::GetFloatArrayElements(array, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetFloatArrayElements");
    }

    return result;
}


jdouble* Context::GetDoubleArrayElements(jdoubleArray array, jboolean* isCopy)
{
    clearException();
    jdouble* result = JNIEnv_::GetDoubleArrayElements(array, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetDoubleArrayElements");
    }

    return result;
}


void Context::ReleaseBooleanArrayElements(jbooleanArray array, jboolean* elems, jint mode)
{
    // release funtions should never throw
    JNIEnv_::ReleaseBooleanArrayElements(array, elems, mode);
}


void Context::ReleaseByteArrayElements(jbyteArray array, jbyte* elems, jint mode)
{
    // release funtions should never throw
    JNIEnv_::ReleaseByteArrayElements(array, elems, mode);
}


void Context::ReleaseCharArrayElements(jcharArray array, jchar* elems, jint mode)
{
    // release funtions should never throw
    JNIEnv_::ReleaseCharArrayElements(array, elems, mode);
}


void Context::ReleaseShortArrayElements(jshortArray array, jshort* elems, jint mode)
{
    // release funtions should never throw
    JNIEnv_::ReleaseShortArrayElements(array, elems, mode);
}


void Context::ReleaseIntArrayElements(jintArray array, jint* elems, jint mode)
{
    // release funtions should never throw
    JNIEnv_::ReleaseIntArrayElements(array, elems, mode);
}


void Context::ReleaseLongArrayElements(jlongArray array, jlong* elems, jint mode)
{
    // release funtions should never throw
    JNIEnv_::ReleaseLongArrayElements(array, elems, mode);
}


void Context::ReleaseFloatArrayElements(jfloatArray array, jfloat* elems, jint mode)
{
    // release funtions should never throw
    JNIEnv_::ReleaseFloatArrayElements(array, elems, mode);
}


void Context::ReleaseDoubleArrayElements(jdoubleArray array, jdouble* elems, jint mode)
{
    // release funtions should never throw
    JNIEnv_::ReleaseDoubleArrayElements(array, elems, mode);
}


void Context::GetBooleanArrayRegion(jbooleanArray array, jsize start, jsize len, jboolean* buf)
{
    clearException();
    JNIEnv_::GetBooleanArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetBooleanArrayRegion");
    }
}


void Context::GetByteArrayRegion(jbyteArray array, jsize start, jsize len, jbyte* buf)
{
    clearException();
    JNIEnv_::GetByteArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetByteArrayRegion");
    }
}


void Context::GetCharArrayRegion(jcharArray array, jsize start, jsize len, jchar* buf)
{
    clearException();
    JNIEnv_::GetCharArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetCharArrayRegion");
    }
}


void Context::GetShortArrayRegion(jshortArray array, jsize start, jsize len, jshort* buf)
{
    clearException();
    JNIEnv_::GetShortArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetShortArrayRegion");
    }
}


void Context::GetIntArrayRegion(jintArray array, jsize start, jsize len, jint* buf)
{
    clearException();
    JNIEnv_::GetIntArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetIntArrayRegion");
    }
}


void Context::GetLongArrayRegion(jlongArray array, jsize start, jsize len, jlong* buf)
{
    clearException();
    JNIEnv_::GetLongArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetLongArrayRegion");
    }
}


void Context::GetFloatArrayRegion(jfloatArray array, jsize start, jsize len, jfloat* buf)
{
    clearException();
    JNIEnv_::GetFloatArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetFloatArrayRegion");
    }
}


void Context::GetDoubleArrayRegion(jdoubleArray array, jsize start, jsize len, jdouble* buf)
{
    clearException();
    JNIEnv_::GetDoubleArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetDoubleArrayRegion");
    }
}


void Context::SetBooleanArrayRegion(jbooleanArray array, jsize start, jsize len, jboolean* buf)
{
    clearException();
    JNIEnv_::SetBooleanArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetBooleanArrayRegion");
    }
}


void Context::SetByteArrayRegion(jbyteArray array, jsize start, jsize len, jbyte* buf)
{
    clearException();
    JNIEnv_::SetByteArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetByteArrayRegion");
    }
}


void Context::SetCharArrayRegion(jcharArray array, jsize start, jsize len, jchar* buf)
{
    clearException();
    JNIEnv_::SetCharArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetCharArrayRegion");
    }
}


void Context::SetShortArrayRegion(jshortArray array, jsize start, jsize len, jshort* buf)
{
    clearException();
    JNIEnv_::SetShortArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetShortArrayRegion");
    }
}


void Context::SetIntArrayRegion(jintArray array, jsize start, jsize len, jint* buf)
{
    clearException();
    JNIEnv_::SetIntArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetIntArrayRegion");
    }
}


void Context::SetLongArrayRegion(jlongArray array, jsize start, jsize len, jlong* buf)
{
    clearException();
    JNIEnv_::SetLongArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetLongArrayRegion");
    }
}


void Context::SetFloatArrayRegion(jfloatArray array, jsize start, jsize len, jfloat* buf)
{
    clearException();
    JNIEnv_::SetFloatArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetFloatArrayRegion");
    }
}


void Context::SetDoubleArrayRegion(jdoubleArray array, jsize start, jsize len, jdouble* buf)
{
    clearException();
    JNIEnv_::SetDoubleArrayRegion(array, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::SetDoubleArrayRegion");
    }
}


jint Context::RegisterNatives(jclass clazz, const JNINativeMethod* methods, jint nMethods)
{
    clearException();
    jint result = JNIEnv_::RegisterNatives(clazz, methods, nMethods);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::RegisterNatives");
    }

    return result;
}


jint Context::UnregisterNatives(jclass clazz)
{
    clearException();
    jint result = JNIEnv_::UnregisterNatives(clazz);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::UnregisterNatives");
    }

    return result;
}


jint Context::MonitorEnter(jobject obj)
{
    clearException();
    // should return 0 on success, not sure if it ever throws
    jint result = JNIEnv_::MonitorEnter(obj);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::MonitorEnter");
    }

    return result;
}


jint Context::MonitorExit(jobject obj)
{
    clearException();
    // should return 0 on success, not sure if it ever throws
    jint result = JNIEnv_::MonitorExit(obj);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::MonitorExit");
    }

    return result;
}


jint Context::GetJavaVM(JavaVM** vm)
{
    clearException();
    jint result = JNIEnv_::GetJavaVM(vm);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetJavaVM");
    }

    return result;
}


void Context::GetStringRegion(jstring str, jsize start, jsize len, jchar* buf)
{
    clearException();
    JNIEnv_::GetStringRegion(str, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStringRegion");
    }
}


void Context::GetStringUTFRegion(jstring str, jsize start, jsize len, char* buf)
{
    clearException();
    JNIEnv_::GetStringUTFRegion(str, start, len, buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStringUTFRegion");
    }
}


void* Context::GetPrimitiveArrayCritical(jarray array, jboolean* isCopy)
{
    clearException();
    void* result = JNIEnv_::GetPrimitiveArrayCritical(array, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetPrimitiveArrayCritical");
    }

    return result;
}


void Context::ReleasePrimitiveArrayCritical(jarray array, void* carray, jint mode)
{
    // release funtions should never throw
    JNIEnv_::ReleasePrimitiveArrayCritical(array, carray, mode);
}


const jchar* Context::GetStringCritical(jstring string, jboolean* isCopy)
{
    clearException();
    const jchar* result = JNIEnv_::GetStringCritical(string, isCopy);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetStringCritical");
    }

    return result;
}


void Context::ReleaseStringCritical(jstring string, const jchar* cstring)
{
    // release funtions should never throw
    JNIEnv_::ReleaseStringCritical(string, cstring);
}


jweak Context::NewWeakGlobalRef(jobject obj)
{
    clearException();
    jweak result = JNIEnv_::NewWeakGlobalRef(obj);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewWeakGlobalRef");
    }

    return result;
}


void Context::DeleteWeakGlobalRef(jweak ref)
{
    // deallocation functions should never throw
    JNIEnv_::DeleteWeakGlobalRef(ref);
}


jboolean Context::ExceptionCheck()
{
    // exception functions should never throw
    return JNIEnv_::ExceptionCheck();
}


jobject Context::NewDirectByteBuffer(void* address, jlong capacity)
{
    clearException();
    jobject result = JNIEnv_::NewDirectByteBuffer(address, capacity);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::NewDirectByteBuffer");
    }

    return result;
}


void* Context::GetDirectBufferAddress(jobject buf)
{
    clearException();
    void* result = JNIEnv_::GetDirectBufferAddress(buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetDirectBufferAddress");
    }

    return result;
}


jlong Context::GetDirectBufferCapacity(jobject buf)
{
    clearException();
    jlong result = JNIEnv_::GetDirectBufferCapacity(buf);

    jthrowable error = JNIEnv_::ExceptionOccurred();
    if (error)
    {
        throwJException(this, error, "Context::GetDirectBufferCapacity");
    }

    return result;
}
