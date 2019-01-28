/*
 * Copyright (c) 2002    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

//////////////////////////////////////////////////////////////////////
// JException.cpp: implementation of the JException class.
//
//////////////////////////////////////////////////////////////////////

#include "JException.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

JException::JException() {
}

JException::JException(char* message) {
    msg.append(message);
}

JException::JException(const char* message) {
    msg.append(message);
}

JException::JException(SlimString& message) {
    msg.append(message.c_str());
}

JException::~JException() {
}

const char* JException::what() const {
    return msg.c_str();
}
