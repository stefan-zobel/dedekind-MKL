# dedekind-MKL

Selected BLAS and LAPACK Java bindings for Intel's MKL (Math Kernel Library)


### What is included?

The current implementation covers 76 prominent functions for column-major double[] arrays / dense matrices (about two-thirds from LAPACK, the rest focusing on level 2 / level 3 operations from BLAS and BLAS extensions). Roughly a dozen of those functions are also available in float[] and complex float / complex double versions adding up to about 110 functions overall.

In particular, the library has implementations for all 4 data types of

* `?gemm, ?imatcopy, ?omatadd, ?omatcopy`
* `?geev`
* `?gels`
* `?geqrf`
* `?gesdd`
* `?gesv`
* `?getrf`

In addition, the special complex functions

* `cgemm3m, zgemm3m`
* `cungqr, zungqr` (complex versions of `sorgqr / dorgqr`)

are also provided. Hence, the basic operations for general dense matrices, as well as the `LU`, `QR`, `EVD` and `SVD` decompositions are covered for all 4 data types (and there is more for double[] only).


### How to use?

The intended purpose of the library is to be a Java JNI wrapper of Intel's native *Math Kernel Library* (`MKL`). However, the float and double functions can also be used without an MKL installation as the pure Java *org.netlib* implementation of those functions is included as a dependency.

Given that the original Fortran API is very low-level, there is [PlainLapack](https://github.com/stefan-zobel/dedekind-MKL/blob/master/src/main/java/net/frobenius/lapack/PlainLapack.java) which is a very thin wrapper that provides automatic work array allocation and extensive parameter checking (only partially implemented by now).
Since parameter checking is crucial for the JNI implementation, the safest approach is to go through the `PlainLapack` API.

Beyond that, no attempt has been made to hide or simplify the difficulties of the original BLAS / LAPACK API.
