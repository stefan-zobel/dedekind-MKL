/*
 * Copyright (c) 2010 - 2011    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

#ifndef THREADAUTOATTACH_INCLUDED_
#define THREADAUTOATTACH_INCLUDED_


#ifndef STDAFX_INCLUDED_
#include "stdafx.h"
#endif /* STDAFX_INCLUDED_ */

#ifndef _JAVASOFT_JNI_H_
#include "jni.h"
#endif /* _JAVASOFT_JNI_H_ */

#ifndef __CONTEXT_H_INCLUDED_
#include "Context.h"
#endif /* __CONTEXT_H_INCLUDED_ */

#ifndef THREADMUTEX_INCLUDED__
#include "ThreadMutex.h"
#endif /* THREADMUTEX_INCLUDED__ */



//////////////////////////////////////////////////////////////////////
// Use the RAII idiom to guarantee deterministic native thread attach
// and detach.
//
// The idea is that the object's destructor is responsible for detaching
// the native thread. When the object is stack allocated, the object's
// destructor will be called whenever the object is popped off the stack
// (when the object cleanly goes out of scope or during stack unwinding
// in the presence of an exception). Thus, the detach of the native thread
// should be guaranteed to occur.
//////////////////////////////////////////////////////////////////////

class __GCC_DONT_EXPORT ThreadAutoAttach {
public:
    ThreadAutoAttach();
    ~ThreadAutoAttach();

    Context* getContext();

    // for prematurely detaching the native thread before the destructor is run
    void prematureDetach();

public: // static methods

    // for prematurely detaching the native thread without a reference
    // to a ThreadAutoAttach instance
    static void prematureDetach(JavaVM* javaVM, const char* logMsg);
    // static methods for disabling and re-enabling of
    // POSIX thread cancellation (only used on Unix/Linux)
    static void disableCancellation();
    static void enableCancellation();

private:
    // prohibit heap allocation (must be stack allocated)
    // Note: these operators shouldn't have any implementation,
    // but we couldn't manage to link a loadable sl with aCC
    // on HP/UX without an implementation!
    void* operator new(size_t) { throw 0; }
    void operator delete(void*, size_t) { throw 0; }
    void* operator new[](size_t) { throw 0; }
    void operator delete[](void*, size_t) { throw 0; }

    // prohibit change of ownership (each object is the exclusive owner of its resources)
    ThreadAutoAttach(const ThreadAutoAttach&); // copy constructor
    ThreadAutoAttach& operator= (const ThreadAutoAttach&); // assignment operator

    // internal thread detach
    void detach();

private: // static method(s)
    static void detach(JavaVM* javaVM, const char* logMsg);

private: // data
    JavaVM* m_vm;
    JNIEnv* m_env;

    // Global mutex to serialize thread access to the Attach/Detach API functions.
    // This is unfortunately necessary on Sun VMs, otherwise they will leak Windows
    // event handles like hell!
    static ThreadMutex globalMutex;
};



#endif /* THREADAUTOATTACH_INCLUDED_ */
