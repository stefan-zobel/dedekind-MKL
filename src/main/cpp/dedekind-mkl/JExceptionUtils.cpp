/*
 * Copyright (c) 2010    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

#include "JExceptionUtils.h"


#ifndef _JAVASOFT_JNI_H_
#include "jni.h"
#endif /* _JAVASOFT_JNI_H_ */

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */

#ifndef __CONTEXT_H_INCLUDED_
#include "Context.h"
#endif /* __CONTEXT_H_INCLUDED_ */

#ifndef POTABILITY_INCLUDED_
#include "Portability.h" // for snprintf_portable
#endif /* POTABILITY_INCLUDED_ */


#include <string.h> // for strncpy



#if defined (_WIN64) || defined (_WIN32)
// disable "This function may be unsafe" warnings for older C API functions
#pragma warning( disable: 4996 ) // instead of _CRT_SECURE_NO_WARNINGS, which doesn't work
#endif /* (_WIN64) || (_WIN32) */




/*
 * Get a human-readable summary of an exception object. The buffer will
 * be populated with the class name and, if present, the exception message.
 */
static bool getExceptionSummary(JNIEnv* env, jthrowable exception, char* buf, size_t bufLen);




/*
 * Formats an exception as a string with its stack trace.
 */
bool printStackTrace(JNIEnv* env, jthrowable exception, char* buf, size_t bufLen)
{
    bool success = false;

    if (env->ExceptionCheck()) {
        env->ExceptionClear();
    }

    jclass stringWriterClazz = env->FindClass("java/io/StringWriter");
    if (stringWriterClazz != NULL) {
        jmethodID stringWriterCtor = env->GetMethodID(stringWriterClazz,
                "<init>", "()V");
        jmethodID stringWriterToStringMethod = env->GetMethodID(stringWriterClazz,
                "toString", "()Ljava/lang/String;");

        jclass printWriterClazz = env->FindClass("java/io/PrintWriter");
        if (printWriterClazz != NULL) {
            jmethodID printWriterCtor = env->GetMethodID(printWriterClazz, "<init>",
                "(Ljava/io/Writer;)V");

            jobject stringWriterObj = env->NewObject(stringWriterClazz, stringWriterCtor);
            if (stringWriterObj != NULL) {
                jobject printWriterObj = env->NewObject(printWriterClazz, printWriterCtor,
                        stringWriterObj);
                if (printWriterObj != NULL) {
                    jclass exceptionClazz = env->GetObjectClass(exception); // can't fail
                    jmethodID printStackTraceMethod = env->GetMethodID(exceptionClazz,
                        "printStackTrace", "(Ljava/io/PrintWriter;)V");

                    env->CallVoidMethod(exception, printStackTraceMethod, printWriterObj);

                    if (! env->ExceptionCheck() ) {
                        jstring messageStr = static_cast<jstring>(env->CallObjectMethod(
                                stringWriterObj, stringWriterToStringMethod));
                        if (messageStr != NULL) {
                            jsize messageStrLength = env->GetStringLength(messageStr);
                            if (messageStrLength >= static_cast<jsize>(bufLen)) {
                                messageStrLength = static_cast<jsize>(bufLen) - 1;
                            }
                            env->GetStringUTFRegion(messageStr, 0, messageStrLength, buf);
                            env->DeleteLocalRef(messageStr);
                            buf[messageStrLength] = '\0';
                            success = true;
                        }
                    }
                    env->DeleteLocalRef(exceptionClazz);
                    env->DeleteLocalRef(printWriterObj);
                }
                env->DeleteLocalRef(stringWriterObj);
            }
            env->DeleteLocalRef(printWriterClazz);
        }
        env->DeleteLocalRef(stringWriterClazz);
    }

    if (!success) {
        if (env->ExceptionCheck()) {
            env->ExceptionClear();
        }
        return getExceptionSummary(env, exception, buf, bufLen);
    }
    return success;
}




/*
 * Formats an exception as a string with its stack trace and additional context information.
 */
bool printStackTrace(JNIEnv* env, jthrowable exception, SlimString& msg, char* context) {
    const int BUFFER_LEN = 4096; // max size of gathered stacktrace is 4K
    char buffer[BUFFER_LEN] = {0};
    msg.append(context).append("\n");
    bool success = printStackTrace(env, exception, buffer, BUFFER_LEN);
    msg.append(buffer);
    return success;
}




/*
 * Get a human-readable summary of an exception object. The buffer will
 * be populated with the class name and, if present, the exception message.
 */
bool getExceptionSummary(JNIEnv* env, jthrowable exception, char* buf, size_t bufLen)
{
    bool success = false;

    /* get the name of the exception's class */
    jclass exceptionClazz = env->GetObjectClass(exception); // can't fail
    jclass classClazz = env->GetObjectClass(exceptionClazz); // java.lang.Class, can't fail
    jmethodID classGetNameMethod = env->GetMethodID(classClazz, "getName", "()Ljava/lang/String;");
    jstring classNameStr = static_cast<jstring>(env->CallObjectMethod(exceptionClazz,
        classGetNameMethod));

    if (classNameStr != NULL) {
        /* get printable string */
        const char* classNameChars = env->GetStringUTFChars(classNameStr, NULL);
        if (classNameChars != NULL) {
            /* if the exception has a message string, get that */
            jmethodID throwableGetMessageMethod = env->GetMethodID(
                    exceptionClazz, "getMessage", "()Ljava/lang/String;");
            jstring messageStr = static_cast<jstring>(env->CallObjectMethod(exception,
                throwableGetMessageMethod));

            if (messageStr != NULL) {
                const char* messageChars = env->GetStringUTFChars(messageStr, NULL);
                if (messageChars != NULL) {
                    snprintf_portable(buf, bufLen, "%s: %s", classNameChars, messageChars);
                    env->ReleaseStringUTFChars(messageStr, messageChars);
                } else {
                    if (env->ExceptionCheck()) {
                        env->ExceptionClear(); // clear potential OOM
                    }
                    snprintf_portable(buf, bufLen, "%s: <error getting message>", classNameChars);
                }
                env->DeleteLocalRef(messageStr);
            } else {
                strncpy(buf, classNameChars, bufLen);
                buf[bufLen - 1] = '\0';
            }

            env->ReleaseStringUTFChars(classNameStr, classNameChars);
            success = true;
        }
        env->DeleteLocalRef(classNameStr);
    }
    env->DeleteLocalRef(classClazz);
    env->DeleteLocalRef(exceptionClazz);

    if (!success) {
        if (env->ExceptionCheck()) {
            env->ExceptionClear();
        }
        snprintf_portable(buf, bufLen, "%s", "<error getting class name>");
    }
    return success;
}




bool throwJavaRuntimeException(JNIEnv* env, const char* format, ...) {
    try {
        Context* pCtx = static_cast<Context*>(env);
        jclass exceptClass = pCtx->FindClass("java/lang/RuntimeException");
        if (exceptClass == NULL) {
            return false;
        }
        const int MAX_MSG_SIZE = 4096;
        char message[MAX_MSG_SIZE] = {0};
        va_list args;
        va_start(args, format);
        int rc = vsprintf(message, format, args);
        va_end(args);
        if (rc < 0) {
            return false;
        }
        rc = pCtx->ThrowNew(exceptClass, message);
        pCtx->DeleteLocalRef(exceptClass);
        if (rc < 0) {
            return false;
        }
    } catch (const JException& /*ignore*/) {
//      __LOG_WARN __LARG("JExceptionUtils::throwJavaRuntimeException : ") __LARG(ignore.what());
        return false;
    } catch (...) {
//      __LOG_WARN __LARG("JExceptionUtils::throwJavaRuntimeException") __LARG(UNEXPECTED_ERR);
        return false;
    }
    return true;
}
