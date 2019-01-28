/*
 * Copyright (c) 2002    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

//////////////////////////////////////////////////////////////////////
// JException.h: interface for the JException class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(JEXCEPTION_INCLUDED_)
#define JEXCEPTION_INCLUDED_

#ifndef SLIMSTRING_INCLUDED_
#include "SlimString.h"
#endif /* SLIMSTRING_INCLUDED_ */

#ifndef STDAFX_INCLUDED_
#include "stdafx.h"
#endif /* STDAFX_INCLUDED_ */




class __GCC_DONT_EXPORT JException {

public:
    JException();
    JException(char* message);
    JException(const char* message);
    JException(SlimString& message);
    virtual ~JException();

    const char* what() const;
    inline const char* name() const {return "JException";}

private:
    SlimString msg;
};

#endif // !defined(JEXCEPTION_INCLUDED_)
