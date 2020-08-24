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

#include "ComplexDoubleArray.h"

#ifndef JEXCEPTION_INCLUDED_
#include "JException.h"
#endif /* JEXCEPTION_INCLUDED_ */



ComplexDoubleArray::ComplexDoubleArray(DoubleArray& array_)
    : array(array_), complex_array(NULL), complex_array_len(0)
{
    long length = array.length();
    if (length > 0) {
        if (length % 2 != 0) {
            throw JException("complex arrays must have even length"); // TODO: add actual length
        }
        length /= 2;
        complex_array_len = length;
        complex_array = (MKL_Complex16*) mkl_malloc(length * sizeof(MKL_Complex16), 64);
        if (!complex_array) {
            throw JException("couldn't allocate MKL_Complex16 array"); // TODO: add requested length
        }
        double* mixed = array.ptr();
        cblas_dcopy(length, &(mixed[0]), 2, &(complex_array[0].real), 2);
        cblas_dcopy(length, &(mixed[1]), 2, &(complex_array[0].imag), 2);
    }
}

MKL_Complex16* ComplexDoubleArray::ptr() {
    return complex_array;
}

long ComplexDoubleArray::complexLength() {
    return complex_array_len;
}

ComplexDoubleArray::~ComplexDoubleArray() {
    if (complex_array) {
        mkl_free(complex_array);
    }
}
