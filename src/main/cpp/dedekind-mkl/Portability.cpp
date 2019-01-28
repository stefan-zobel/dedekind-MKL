/*
 * Copyright (c) 2010    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

#include "Portability.h"

#include <string.h>
#include <time.h>       // _tzset or gettimeofday
#include <stdio.h>      // vsnprintf



#if !defined (_WIN64) && !defined (_WIN32)
#include <sys/time.h>       // gettimeofday
#include <stdarg.h>         // va_list for vsnprintf (on Unix)

#if defined(__linux)
#include <sys/syscall.h>    // SYS_gettid
#include <unistd.h>         // syscall
#endif /* __linux */


#if !defined (__GNUG__)
// GCC no longer implements <varargs.h>
#include <varargs.h>        // Unix V compatibility
#endif /* !__GNUG__ */


#if defined (__hpux)
#include <pthread.h>        // pthread_self
#endif /* __hpux */

#endif



#if defined (_WIN64) || defined (_WIN32)
// disable "This function may be unsafe" warnings for older C API functions
#pragma warning( disable: 4996 ) // instead of _CRT_SECURE_NO_WARNINGS, which doesn't work
#include "windows.h"

#if defined (_MSC_VER) || defined (_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif /* (_MSC_VER) || (_MSC_EXTENSIONS) */

#else /* Unix or Linux */

#endif /* (_WIN64) || (_WIN32) */




const SlimString getCurrentThreadId_portable() {
    const int MAX_LEN = 64;
    char id_[MAX_LEN] = {0};

#if defined (_WIN64) || defined (_WIN32)
    sprintf(id_, "%d", GetCurrentThreadId());
    return SlimString(id_);
#elif defined (__linux) && defined (__GNUG__)
    sprintf(id_, "%ld", syscall( SYS_gettid ));
    return SlimString(id_);
#elif defined (__hpux)
    sprintf(id_, "%d", pthread_self());
    return SlimString(id_);
#else
    return SlimString("0");
#endif

}




#ifdef __cplusplus
extern "C" {
#endif



size_t strnlen_portable(const char* str, size_t maxsize) {

#if defined (_WIN64) || defined (_WIN32)

    return strnlen(str, maxsize);

#else
    size_t n;

    /* Note that we do not check if s == NULL ! */

    for (n = 0; n < maxsize && *str; n++, str++)
        ;

    return n;
#endif

}


// forward to vsnprintf which should also be available on Unix/Linux.
int snprintf_portable(char* buffer, size_t count, const char* format, ...) {
    va_list argslist;
    va_start(argslist, format);
    int num = vsnprintf(buffer, count, format, argslist);
    va_end(argslist);
    return num;
}



int gettimeofday_portable(struct timeval* tv, struct timezone* tz) {

#if !defined (_WIN64) && !defined (_WIN32)
    return gettimeofday(tv, tz);
#else

    FILETIME ft;
    unsigned __int64 tmpres = 0;
    static int tzflag = 0;

    if (tv) {
        GetSystemTimeAsFileTime(&ft);

        tmpres |= ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |= ft.dwLowDateTime;

        // convert file time to unix epoch
        tmpres /= 10;  // convert into microseconds
        tmpres -= DELTA_EPOCH_IN_MICROSECS;

        tv->tv_sec = (long) (tmpres / 1000000UL);
        tv->tv_usec = (long) (tmpres % 1000000UL);
    }

    if (tz) {
        if (!tzflag) {
            _tzset();
            ++tzflag;
        }
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
    }

    return 0;

#endif

}



#ifdef __cplusplus
}
#endif /* __cplusplus */
