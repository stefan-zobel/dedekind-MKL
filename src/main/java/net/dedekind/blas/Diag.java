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

public enum Diag {

    N(131), //
    U(132); //

    public int code() {
        return diag;
    }

    private final int diag;

    private Diag(int diag) {
        this.diag = diag;
    }

    public static Diag of(String diag) {
        switch (Character.toUpperCase(diag.charAt(0))) {
        case 'N':
            return N;
        case 'U':
            return U;
        default:
            throw new IllegalArgumentException("Invalid DIAG: " + diag);
        }
    }
}
