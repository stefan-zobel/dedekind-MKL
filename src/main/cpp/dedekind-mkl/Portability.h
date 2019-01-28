/*
 * Copyright (c) 2010    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

#ifndef POTABILITY_INCLUDED_
#define POTABILITY_INCLUDED_


#include <stddef.h> // size_t



#ifndef STDAFX_INCLUDED_
#include "stdafx.h"
#endif /* STDAFX_INCLUDED_ */

#ifndef SLIMSTRING_INCLUDED_
#include "SlimString.h"
#endif /* SLIMSTRING_INCLUDED_ */

#if defined (_WIN64) || defined (_WIN32)
#include <WinSock2.h>
#endif /* (_WIN64) || (_WIN32) */




#if defined (__linux)

#define CATCHALL catch (const abi::__forced_unwind&) {\
    throw;\
} catch (...)

#else

#define CATCHALL catch (...)

#endif




__GCC_DONT_EXPORT const SlimString getCurrentThreadId_portable();




#ifdef __cplusplus
extern "C" {
#endif



#if defined (_WIN64) || defined (_WIN32)


struct timezone 
{
  int  tz_minuteswest; /* minutes west of Greenwich */
  int  tz_dsttime;     /* type of DST correction */
};


#endif /* (_WIN64) || (_WIN32) */



/**
 * strnlen_portable - return the length of a null-terminated string
 *
 * Purpose:
 *   Finds the length in bytes of the given string, not including
 *   the final null character. Only the first maxsize characters
 *   are inspected: if the null character is not found, maxsize is
 *   returned.
 *
 * Entry:
 *   const char* str - string whose length is to be computed
 *   size_t maxsize
 *
 * Exit:
 *   Length of the string "str", exclusive of the final null byte, or
 *   maxsize if the null character is not found.
 *
 * Exceptions:
 *   Access violation / segfault if "str" is NULL
 *
*******************************************************************************/
__GCC_DONT_EXPORT size_t strnlen_portable(const char* str, size_t maxsize);



__GCC_DONT_EXPORT int snprintf_portable(char* buffer, size_t count, const char* format, ...);



__GCC_DONT_EXPORT int gettimeofday_portable(struct timeval* tv, struct timezone* tz);



#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* POTABILITY_INCLUDED_ */

