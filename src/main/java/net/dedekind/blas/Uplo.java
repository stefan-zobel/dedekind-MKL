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
package net.dedekind.blas;

public enum Uplo {

    U(121), //
    L(122); //

    public int code() {
        return uplo;
    }

    private final int uplo;

    private Uplo(int uplo) {
        this.uplo = uplo;
    }

    public static Uplo of(String uplo) {
        switch (Character.toUpperCase(uplo.charAt(0))) {
        case 'U':
            return U;
        case 'L':
            return L;
        default:
            throw new IllegalArgumentException("Invalid UPLO: " + uplo);
        }
    }
}
