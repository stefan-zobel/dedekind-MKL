/*
 * Copyright 2019 Stefan Zobel
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

#ifndef DEDEKIND_INCLUDED_
#define DEDEKIND_INCLUDED_

#ifndef __CONTEXT_H_INCLUDED_
#include "Context.h"
#endif /* __CONTEXT_H_INCLUDED_ */


#ifdef __cplusplus
extern "C" {
#endif

    __GCC_DONT_EXPORT void initializeDescriptors(JNIEnv* env);

    __GCC_DONT_EXPORT void setIntWValue(JNIEnv* env, jobject intW, int value);

    __GCC_DONT_EXPORT void setDoubleWValue(JNIEnv* env, jobject doubleW, double value);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* DEDEKIND_INCLUDED_ */
