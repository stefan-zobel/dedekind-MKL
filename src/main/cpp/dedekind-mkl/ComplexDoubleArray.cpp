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

#ifndef SLIMSTRING_INCLUDED_
#include "SlimString.h"
#endif /* SLIMSTRING_INCLUDED_ */

#include <string.h> // for memcpy



__GCC_DONT_EXPORT void doubleCopy(long len, double* mixed, MKL_Complex16* complex) {
    if (len > 0 && mixed && complex) {
        memcpy(mixed, complex, len * sizeof(MKL_Complex16));
//        cblas_dcopy(len, &(complex[0].real), 2, &(mixed[0]), 2);
//        cblas_dcopy(len, &(complex[0].imag), 2, &(mixed[1]), 2);
    }
}

ComplexDoubleArray::ComplexDoubleArray(DoubleArray& array_, bool copy)
    : array(array_), complex_array_len(0), complex_array(NULL), isCopy(copy)
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
        if (isCopy) {
            complex_array = (MKL_Complex16*) mkl_malloc(length * sizeof(MKL_Complex16), 64);
            if (!complex_array) {
                SlimString msg("couldn't allocate MKL_Complex16 array of length ");
                msg.append(length);
                throw JException(msg);
            }
            memcpy(complex_array, array.ptr(), length * sizeof(MKL_Complex16));
//            double* mixed = array.ptr();
//            cblas_dcopy(length, &(mixed[0]), 2, &(complex_array[0].real), 2);
//            cblas_dcopy(length, &(mixed[1]), 2, &(complex_array[0].imag), 2);
        }
        else {
            complex_array = (MKL_Complex16*) array.ptr();
        }
    }
}

MKL_Complex16* ComplexDoubleArray::ptr() {
    return complex_array;
}

long ComplexDoubleArray::complexLength() {
    return complex_array_len;
}

bool ComplexDoubleArray::hasCopy() {
    return isCopy;
}

ComplexDoubleArray::~ComplexDoubleArray() {
    if (isCopy && complex_array) {
        mkl_free(complex_array);
    }
}
