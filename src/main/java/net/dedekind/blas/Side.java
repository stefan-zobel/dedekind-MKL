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

public enum Side {

    L(141), //
    R(142); //

    public int code() {
        return side;
    }

    private final int side;

    private Side(int side) {
        this.side = side;
    }

    public static Side of(String side) {
        switch (Character.toUpperCase(side.charAt(0))) {
        case 'L':
            return L;
        case 'R':
            return R;
        default:
            throw new IllegalArgumentException("Invalid side: " + side);
        }
    }
}
