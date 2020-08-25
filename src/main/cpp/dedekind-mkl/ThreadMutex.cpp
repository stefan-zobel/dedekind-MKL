/*
 * Copyright (c) 2002           Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

//////////////////////////////////////////////////////////////////////
// ThreadMutex.cpp: implementation of the ThreadMutex class.
//
//////////////////////////////////////////////////////////////////////

#include "ThreadMutex.h"

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */

#ifndef SLIMSTRING_INCLUDED_
#include "SlimString.h"
#endif /* SLIMSTRING_INCLUDED_ */




ThreadMutex::ThreadMutex() {
#if defined (_WIN64) || defined (_WIN32)

    // Note: The spincount will be ignored on single-processor machines
    if (!InitializeCriticalSectionAndSpinCount(&m_cs, SPINCOUNT)) {
        // this is unexpected
        unsigned long err = GetLastError();
        SlimString msg("ThreadMutex::ThreadMutex() failed unexpectedly (low memory situation?) - Error Code: ");
        msg.append(err);
        throw JException(msg);
    }

#else /* POSIX threads */

    int rc = pthread_mutex_init(&m_cs, NULL);

    if (rc != 0) {
        SlimString msg("ThreadMutex::ThreadMutex() failed unexpectedly with error code: ");
        msg.append(rc);
        throw JException(msg);
    }

#endif /* (_WIN64) || (_WIN32) */
}


void ThreadMutex::acquire() const {
#if defined (_WIN64) || defined (_WIN32)

    __try {
        EnterCriticalSection(&m_cs);
    }
    __except (GetExceptionCode() == EXCEPTION_INVALID_HANDLE
    ? EXCEPTION_EXECUTE_HANDLER : EXCEPTION_CONTINUE_SEARCH) {
        // On Windows 2000 EXCEPTION_INVALID_HANDLE could happen
        // in very low memory situations. I translate this into a C++
        // exception here (Note that due to EXCEPTION_CONTINUE_SEARCH
        // any other SE will be 'rethrown' untranslated!).
        // On XP or later EnterCriticalSection is guaranteed not to
        // fail due to lack of resources.
        const char* pszMsg = "ThreadMutex::acquire() failed with EXCEPTION_INVALID_HANDLE"
            " - this may be an indicator of a very low memory situation";
        throw JException(pszMsg);
    }

#else /* POSIX threads*/

    int rc = pthread_mutex_lock(&m_cs);

    if (rc != 0) {
        SlimString msg("ThreadMutex::acquire() failed with return code ");
        msg.append(rc);
        throw JException(msg);
    }

#endif /* (_WIN64) || (_WIN32) */
}
