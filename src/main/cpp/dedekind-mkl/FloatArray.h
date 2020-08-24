/*
 * Copyright 2019, 2020 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef FLOATARRAY_INCLUDED_
#define FLOATARRAY_INCLUDED_

#ifndef STDAFX_INCLUDED_
#include "stdafx.h"
#endif /* STDAFX_INCLUDED_ */

#ifndef __CONTEXT_H_INCLUDED_
#include "Context.h"
#endif /* __CONTEXT_H_INCLUDED_ */

class __GCC_DONT_EXPORT FloatArray
{
public:
    FloatArray(JNIEnv* env, jfloatArray jarray, int offset, jboolean critical);
    ~FloatArray();
    float* ptr();
    long length();
private:
    Context* ctx;
    jfloatArray jarray;
    float* carray;
    int offset;
    bool critical;
};

#endif /* FLOATARRAY_INCLUDED_ */
