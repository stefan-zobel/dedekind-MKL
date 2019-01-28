/*
 * Copyright (c) 2010 - 2011    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

#include "ThreadAutoAttach.h"

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */

#ifndef JVMPROVIDER_INCLUDED_
#include "JvmProvider.h"
#endif /* JVMPROVIDER_INCLUDED_ */

#ifndef LOGGER_INCLUDED_
//#include "Logger.h"
#endif /* LOGGER_INCLUDED_ */


#if !defined (_WIN64) && !defined (_WIN32) // POSIX threads
#include <pthread.h>
#endif /* !(_WIN64) && !(_WIN32) */


#if defined (__GNUG__)
// disable warning "deprecated conversion from string constant to 'char*'" on G++ compiler
#pragma GCC diagnostic ignored "-Wwrite-strings"
#endif /* __GNUG__ */




// static global mutex for serialization of thread's access to the Attach/Detach API functions
// (which are _not_ thread-safe on Sun VMs!) 
ThreadMutex ThreadAutoAttach::globalMutex;




// constructor
ThreadAutoAttach::ThreadAutoAttach()
 : m_vm(NULL), m_env(NULL)
{
    m_vm = JvmProvider::instance()->getJavaVM();
    if (m_vm) {
        JavaVMAttachArgs attachArgs = {JNI_VERSION_1_6, "MKL native thread", NULL};
        jint rc = 0;
        try {
            // disable POSIX thread cancellation if on Unix
            disableCancellation();
            // the Attach/Detach calls themselves (not the JNI calls in between) need to
            // be serialized on Sun VMs, otherwise Windows will leak event handles like crazy!
            ThreadAutoAttach::globalMutex.acquire();
            rc = m_vm->AttachCurrentThreadAsDaemon(reinterpret_cast<void**>(&m_env), &attachArgs);
            ThreadAutoAttach::globalMutex.release();
        } catch (const JException& e) {
            m_vm = NULL;
            throw e;
        } catch (...) {
            m_vm = NULL;
            throw JException("ThreadAutoAttach::ThreadAutoAttach() - JNI AttachCurrentThreadAsDaemon failed at catch(...)");
        }
        if (rc != 0) {
            m_vm = NULL;
            throw JException("ThreadAutoAttach::ThreadAutoAttach() - JNI AttachCurrentThreadAsDaemon failed at rc != 0");
        }
        if (!m_env) {
            detach();
            throw JException("ThreadAutoAttach::ThreadAutoAttach() - no JNI Environment (JNIEnv) acquired");
        }
    }
}


// destructor
ThreadAutoAttach::~ThreadAutoAttach()
{
    if (m_vm) {
        detach();
    } else {
//        __LOG_ERROR __LARG("ThreadAutoAttach::~ThreadAutoAttach() - m_vm member is NULL");
    }
}


// public instance methods
Context* ThreadAutoAttach::getContext() {
    if (m_env) {
        return static_cast<Context*>(m_env);
    } else {
        throw JException("ThreadAutoAttach::getContext - unexpected: JNIEnv* member is NULL");
    }
}

void ThreadAutoAttach::prematureDetach() {
    if (m_vm) {
        // delegate to private static detach method
        detach(m_vm, "(premature detach)");
    }
}


// public static method(s)
void ThreadAutoAttach::prematureDetach(JavaVM* javaVM, const char* logMsg) {
    if (javaVM) {
        // delegate to private static detach method
        detach(javaVM, logMsg);
    }
}


void ThreadAutoAttach::disableCancellation() {

#if !defined (_WIN64) && !defined (_WIN32)
    int oldState = 0; // unused
    int rc = pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldState);
    if (rc != 0) {
        const char* errMsg = "ThreadAutoAttach::disableCancellation() - failed to disable thread cancellation!";
        __LOG_FATAL __LARG(errMsg);
        throw JException(errMsg);
    }
#endif /* !(_WIN64) && !(_WIN32) */

}


void ThreadAutoAttach::enableCancellation() {

#if !defined (_WIN64) && !defined (_WIN32)
    int oldState = 0; // unused
    int rc = pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &oldState);
    if (rc != 0) {
        const char* warn = "ThreadAutoAttach::enableCancellation() - failed to re-enable thread cancellation!";
        __LOG_WARN __LARG(warn);
        // never throw from a destructor!
    }
#endif /* !(_WIN64) && !(_WIN32) */

}


// private static method(s)
void ThreadAutoAttach::detach(JavaVM* javaVM, const char* /*logMsg*/) {
    if (javaVM) {
        try {
            // the Attach/Detach calls themselves (not the JNI calls in between) need to
            // be serialized on Sun VMs, otherwise Windows will leak event handles like crazy!
            ThreadAutoAttach::globalMutex.acquire();
            jint rc = javaVM->DetachCurrentThread();
            ThreadAutoAttach::globalMutex.release();
            if (rc != 0) {
//                __LOG_WARN __LARG("ThreadAutoAttach::detach() - failed to DetachCurrentThread (rc != JNI_OK) ") __LARG(logMsg);
            }/* else if (Logger::isLogLevelEnabled(LOG_FINER)) {
//                __LOG_FINER __LARG("ThreadAutoAttach::detach() successful ") __LARG(logMsg);
            }*/
            // re-enable POSIX thread cancellation if on Unix
            enableCancellation();
        } catch (const JException& /*e*/) {
//            __LOG_FATAL __LARG("ThreadAutoAttach::detach() - unable to DetachCurrentThread() : ") __LARG(e.what());
        } catch (...) {
//            __LOG_FATAL __LARG("ThreadAutoAttach::detach() - failed to DetachCurrentThread (catch ...)");
        }
    }
}


// private methods
void ThreadAutoAttach::detach() {
    if (m_vm) {
        // delegate to private static detach method
        detach(m_vm, "(destructor)");
        // and clear members
        m_vm = NULL;
        m_env = NULL;
    }
}
