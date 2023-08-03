[![CodeQL](https://github.com/stefan-zobel/dedekind-MKL/actions/workflows/codeql.yml/badge.svg)](https://github.com/stefan-zobel/dedekind-MKL/actions/workflows/codeql.yml)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/net.sourceforge.streamsupport/dedekind-mkl/badge.svg)](https://maven-badges.herokuapp.com/maven-central/net.sourceforge.streamsupport/dedekind-mkl)
[![javadoc.io](https://javadoc.io/badge2/net.sourceforge.streamsupport/dedekind-mkl/javadoc.svg)](https://javadoc.io/doc/net.sourceforge.streamsupport/dedekind-mkl)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# dedekind-MKL

Selected BLAS and LAPACK Java bindings for Intel's oneAPI Math Kernel Library on Windows and Linux


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

The intended purpose of the library is to be a Java JNI wrapper of Intel's native *oneAPI Math Kernel Library* (`oneMKL`). However, the float and double functions can also be used without an MKL installation as the pure Java *org.netlib* implementation of those functions is included as a dependency.

Given that the original Fortran API is very low-level, there is [PlainLapack](https://github.com/stefan-zobel/dedekind-MKL/blob/master/src/main/java/net/frobenius/lapack/PlainLapack.java) which is a very thin wrapper that provides automatic work array allocation and extensive parameter checking (the latter being only partially implemented by now).
Since parameter checking is crucial for the JNI implementation, the safest approach is to go through the `PlainLapack` API.

Beyond that, no attempt has been made to hide or simplify the difficulties of the original BLAS / LAPACK API.


### Maven:

```xml
<dependency>
    <groupId>net.sourceforge.streamsupport</groupId>
    <artifactId>dedekind-mkl</artifactId>
    <version>1.1.0</version>
</dependency>
```

### Enabling the Intel MKL library

*dedekind-MKL* doesn't ship with the Intel [oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) library. If MKL isn't already installed on your system an easy path is to download one of the redistributables [here](https://repo1.maven.org/maven2/org/bytedeco/mkl/) or [here](https://github.com/Anlon-Burke/intel-mkl-x64-redist/releases), extract the archive to a directory and add that to your `PATH` variable and you should be ready to go. There is also a [nuget](https://www.nuget.org/packages/intelmkl.redist.win-x64/) redistributable for Windows and a lot of Linux distros ship MKL in their repositories, for example, on *Arch Linux*, it's simply `pacman -S intel-oneapi-mkl` and you're done. Of course, you can also download the [official](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) installer from Intel for a full installation which has a **huge** on-disk footprint (~ 8 GiB). Finally, the MKL binaries included in Python distributions like [Anaconda](https://www.anaconda.com/products/individual) or [WinPython](https://winpython.github.io/) should also be usable.

By default, dedekind-MKL searches for its own *dedekind-mkl.so/.dll* native library on the *java.library.path* or in a directory that can be specified either through the environment variable `DEDEKIND_SHAREDLIB_DIR` or via the Java system property `dedekind.sharedlib.dir`. If that doesn't succeed, the shared library gets automatically unpacked from the jar file into *java.io.tmpdir* and loaded from there. The *dedekind-mkl.so/.dll* can <ins>only</ins> be loaded if the OS can find the dependent `MKL` libraries somewhere on the `PATH`.

You'll see an exception stack trace in the console if that doesn't work but *dedekind-MKL* will still be (partially) usable as it will fallback to the [F2J](https://repo1.maven.org/maven2/net/sourceforge/f2j/arpack_combined_all/0.1/) pure Java implementation. However, none of the BLAS / LAPACK routines for complex datatypes nor the BLAS extensions will be usable in that case and you'll surrender major performance improvements especially for larger matrices.

 