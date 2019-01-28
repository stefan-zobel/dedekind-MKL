/*
 * Copyright (c) 2002           Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

//////////////////////////////////////////////////////////////////////
// ThreadMutex.h: interface for the ThreadMutex class.
//
//////////////////////////////////////////////////////////////////////

#ifndef THREADMUTEX_INCLUDED__
#define THREADMUTEX_INCLUDED__


#ifndef STDAFX_INCLUDED_
#include "stdafx.h"
#endif /* STDAFX_INCLUDED_ */


// POSIX threads
#if !defined (_WIN64) && !defined (_WIN32)

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */

#include <stdio.h> // sprintf

#include <pthread.h>
#endif




class __GCC_DONT_EXPORT ThreadMutex {
public:

    // constructor
    ThreadMutex();

    // destructor
    inline ~ThreadMutex() {

#if defined (_WIN64) || defined (_WIN32)
        DeleteCriticalSection(&m_cs);
#else /* POSIX threads */
        pthread_mutex_destroy(&m_cs);
#endif /* (_WIN64) || (_WIN32) */

    }

    // acquire() method
    void acquire() const;

    // release() method
    inline void release() const {

#if defined (_WIN64) || defined (_WIN32)
        LeaveCriticalSection(&m_cs);
#else /* POSIX threads */
        int rc = pthread_mutex_unlock(&m_cs);
        if (rc != 0) {
            const int MAX_LEN = 64;
            char code[MAX_LEN] = {0};
            sprintf(code, "%d", rc);
            const char* errMsg = "ThreadMutex::release() failed unexpectedly"
                " with error code: ";
            SlimString msg(errMsg);
            msg.append(code);
            throw JException(msg);
        }
#endif /* (_WIN64) || (_WIN32) */

    }


#if defined (_WIN64) || defined (_WIN32)
    // Additional member functions on Windows

    /**
     * @return true if the mutex could be acquired, false if not.
     */
    inline bool tryAcquire() const {return !!TryEnterCriticalSection(&m_cs);}

    /**
     * @return true if the current thread is the owner of this
     *         mutex, false if the mutex is owned by another thread
     *         or is not owned by any thread.
     */
    inline bool threadHoldsLock() const {
        return reinterpret_cast<DWORD>(m_cs.OwningThread) == GetCurrentThreadId();
    }

#endif /* (_WIN64) || (_WIN32) */

private:
    // prevent assignment and copying of a ThreadMutex
    ThreadMutex(const ThreadMutex&);
    ThreadMutex& operator=(const ThreadMutex&);

private:

#if defined (_WIN64) || defined (_WIN32)

    // the spincount (only used on SMP systems)
    enum {SPINCOUNT = 4000uL};

    // Windows-specific locking primitive
    mutable CRITICAL_SECTION m_cs;

#else /* POSIX threads */

    mutable pthread_mutex_t m_cs;

#endif /* (_WIN64) || (_WIN32) */

};

#endif // THREADMUTEX_INCLUDED__
