/*
 * Copyright (c) 2002 - 2004    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

//////////////////////////////////////////////////////////////////////
// SlimString.cpp: implementation of the SlimString class.
//
//////////////////////////////////////////////////////////////////////


#include "SlimString.h"

#include <stdio.h>  // for sprintf
#include <string.h> // for memcpy, strlen


#if defined (_WIN64) || defined (_WIN32)
// disable "This function may be unsafe" warning for sprintf()
#pragma warning( disable: 4996 )
#endif /* (_WIN64) || (_WIN32) */




SlimString& SlimString::append(const char* s)
{
    return append(const_cast<char*>(s));
}

SlimString& SlimString::append(char* s)
{
    if (s) {
        size_t len1 = 1; // initialize in case that m_pStr == NULL
        const size_t len2 = strlen(s) + 1;
        char* tmp;
        if (m_pStr) {
            len1 = strlen(m_pStr) + 1;
        }
        tmp = new char[len1 + len2 - 1];
        if (m_pStr) {
            memcpy(tmp, m_pStr, len1 * sizeof(char));
            delete [] m_pStr;
        }
        memcpy(tmp + (len1 - 1), s, len2 * sizeof(char));
        m_pStr = tmp;
    }
    return *this;
}

SlimString& SlimString::append(long value)
{
    char buffer[64] = {0};
    sprintf(buffer, "%ld", value);
    return this->append(buffer);
}

void SlimString::copy(const char* s)
{
    size_t len = strlen(s) + 1;
    m_pStr = new char[len];
    memcpy(m_pStr, s, len * sizeof(char));
}
