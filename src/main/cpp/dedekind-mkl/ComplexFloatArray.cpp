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

#include "ComplexFloatArray.h"

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */

#ifndef SLIMSTRING_INCLUDED_
#include "SlimString.h"
#endif /* SLIMSTRING_INCLUDED_ */



ComplexFloatArray::ComplexFloatArray(FloatArray& array_)
    : array(array_), complex_array(NULL), complex_array_len(0)
{
    long length = array.length();
    if (length > 0) {
        if (length % 2 != 0) {
            SlimString msg("complex arrays must have even length: ");
            msg.append(length);
            throw JException(msg);
        }
        length /= 2;
        complex_array_len = length;
        complex_array = (MKL_Complex8*) mkl_malloc(length * sizeof(MKL_Complex8), 64);
        if (!complex_array) {
            SlimString msg("couldn't allocate MKL_Complex8 array of length ");
            msg.append(length);
            throw JException(msg);
        }
        float* mixed = array.ptr();
        cblas_scopy(length, &(mixed[0]), 2, &(complex_array[0].real), 2);
        cblas_scopy(length, &(mixed[1]), 2, &(complex_array[0].imag), 2);
    }
}

MKL_Complex8* ComplexFloatArray::ptr() {
    return complex_array;
}

long ComplexFloatArray::complexLength() {
    return complex_array_len;
}

ComplexFloatArray::~ComplexFloatArray() {
    if (complex_array) {
        mkl_free(complex_array);
    }
}
