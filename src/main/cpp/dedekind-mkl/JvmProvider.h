/*
 * Copyright (c) 2010    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

#ifndef JVMPROVIDER_INCLUDED_
#define JVMPROVIDER_INCLUDED_

#ifndef STDAFX_INCLUDED_
#include "stdafx.h"
#endif /* STDAFX_INCLUDED_ */

#ifndef _JAVASOFT_JNI_H_
#include "jni.h"
#endif /* _JAVASOFT_JNI_H_ */



class __GCC_DONT_EXPORT JvmProvider {
public:
    ~JvmProvider();

    static JvmProvider* instance();
    static void clear();

    JavaVM* getJavaVM();
    void initializeJavaVM(JavaVM* const jvm);
    void clearOnJavaVMUnload();

private: // constructor
    JvmProvider();

private: // data
    static JvmProvider* m_instance;
    JavaVM* m_vm;
};


#endif /* JVMPROVIDER_INCLUDED_ */
