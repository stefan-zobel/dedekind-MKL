/*
 * Copyright 2020 Stefan Zobel
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

#ifndef COMPLEX_DOUBLEARRAY_INCLUDED_
#define COMPLEX_DOUBLEARRAY_INCLUDED_

#ifndef STDAFX_INCLUDED_
#include "stdafx.h"
#endif /* STDAFX_INCLUDED_ */

#ifndef DOUBLEARRAY_INCLUDED_
#include "DoubleArray.h"
#endif /* DOUBLEARRAY_INCLUDED_ */

#ifndef _MKL_H_
#include <mkl.h>
#endif /* _MKL_H_ */


__GCC_DONT_EXPORT void doubleCopy(long len, double* mixed, MKL_Complex16* complex);

class __GCC_DONT_EXPORT ComplexDoubleArray
{
public:
    ComplexDoubleArray(DoubleArray& array_, bool copy = false);
    ~ComplexDoubleArray();
    MKL_Complex16* ptr();
    long complexLength();
    bool hasCopy();
private:
    ComplexDoubleArray& operator=(const ComplexDoubleArray&);
    DoubleArray& array;
    long complex_array_len;
    MKL_Complex16* complex_array;
    bool isCopy;
};

#endif /* COMPLEX_DOUBLEARRAY_INCLUDED_ */
