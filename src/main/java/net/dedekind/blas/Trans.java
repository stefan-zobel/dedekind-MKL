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

public enum Trans {

    N(111), //
    T(112), //
    C(113); //

    public int code() {
        return trans;
    }

    private final int trans;

    private Trans(int trans) {
        this.trans = trans;
    }

    public static Trans of(String trans) {
        switch (Character.toUpperCase(trans.charAt(0))) {
        case 'N':
            return N;
        case 'T':
            return T;
        case 'C':
            return C;
        default:
            throw new IllegalArgumentException("Invalid transposition: " + trans);
        }
    }
}
