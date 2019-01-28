/*
 * Copyright (c) 2010    Stefan Zobel
 *
 * http://www.opensource.org/licenses/mit-license.php
 */

#include "JvmProvider.h"
#include "JException.h"


// static data
JvmProvider* JvmProvider::m_instance = NULL;


// constructor
JvmProvider::JvmProvider()
 : m_vm(NULL)
{
}

// destructor
JvmProvider::~JvmProvider()
{
    m_vm = NULL;
}

// member methods
JavaVM* JvmProvider::getJavaVM() {
    if (!m_vm) {
        throw JException("JvmProvider::getJavaVM - JavaVM* member is not initialized");
    }
    return m_vm;
}

void JvmProvider::initializeJavaVM(JavaVM* const jvm) {
    if (!jvm) {
        throw JException("JvmProvider::initializeJavaVM - NULL JavaVM* argument provided");
    }
    if (m_vm) {
        throw JException("JvmProvider::initializeJavaVM - JavaVM* member is already initialized");
    }
    m_vm = jvm;
}

void JvmProvider::clearOnJavaVMUnload() {
    m_vm = NULL;
}

// static member methods
JvmProvider* JvmProvider::instance() {
    if (m_instance == NULL) {
        m_instance = new JvmProvider();
    }
    return m_instance;
}

void JvmProvider::clear() {
    if (m_instance) {
        delete m_instance;
        m_instance = NULL;
    }
}

