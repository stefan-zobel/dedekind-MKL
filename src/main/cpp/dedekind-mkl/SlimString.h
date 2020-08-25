/*
 * Copyright (c) 2002 - 2004    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

//////////////////////////////////////////////////////////////////////
// SlimString.h: interface for the SlimString class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SLIMSTRING_INCLUDED_)
#define SLIMSTRING_INCLUDED_

#ifndef STDAFX_INCLUDED_
#include "stdafx.h"
#endif /* STDAFX_INCLUDED_ */



class __GCC_DONT_EXPORT SlimString  
{
public:
    explicit SlimString(const char* s = 0);
    SlimString(const SlimString&);
    void operator=(const char*);
    SlimString& append(const char*);
    SlimString& append(char*);
    SlimString& append(long);
    inline const char* c_str() const {return m_pStr;}
    inline ~SlimString() {if (m_pStr) delete [] m_pStr;}

private:
    void copy(const char*);

    SlimString& operator=(const SlimString&);

private:
    char* m_pStr;
};



inline
SlimString::SlimString(const char* s)
 : m_pStr(0)
{
    if (s) {
        copy(s);
    }
}

inline
SlimString::SlimString(const SlimString& s)
 : m_pStr(0)
{
    if (s.m_pStr) {
        copy(s.m_pStr);
    }
}

inline
void SlimString::operator=(const char* s)
{
    if (m_pStr != s) {
        if (m_pStr) {
            delete [] m_pStr;
        }
        if (s) {
            copy(s);
        } else {
            m_pStr = 0;
        }
    }
}

#endif // !defined(SLIMSTRING_INCLUDED_)
