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
package net.incubator.banach.matrix;

public abstract class DimensionsBase implements Dimensions {

    protected final int rows;
    protected final int cols;

    public DimensionsBase(int rows, int cols) {
        checkRows(rows);
        checkCols(cols);
        this.rows = rows;
        this.cols = cols;
    }

    @Override
    public boolean isScalar() {
        return rows == 1 && cols == 1;
    }

    @Override
    public boolean isColumnVector() {
        return cols == 1;
    }

    @Override
    public boolean isRowVector() {
        return rows == 1;
    }

    @Override
    public boolean isSquareMatrix() {
        return rows == cols;
    }

    @Override
    public int numColumns() {
        return cols;
    }

    @Override
    public int numRows() {
        return rows;
    }

    @Override
    public void checkIndex(int i, int j) {
        checkIJ(i, j);
    }

    @Override
    public void checkSubmatrixIndexes(int rFrom, int cFrom, int rTo, int cTo) {
        checkSubmatrixIndices(rFrom, cFrom, rTo, cTo);
    }

    protected int idx(int row, int col) {
        return col * rows + row;
    }

    protected void checkIJ(int i, int j) {
        if (i < 0 || i >= rows) {
            throw new IllegalArgumentException("Illegal row index " + i + " in (" + rows + " x " + cols + ") matrix");
        }
        if (j < 0 || j >= cols) {
            throw new IllegalArgumentException(
                    "Illegal column index " + j + " in (" + rows + " x " + cols + ") matrix");
        }
    }

    protected static int checkRows(int rows) {
        if (rows <= 0) {
            throw new IllegalArgumentException("number of rows must be strictly positive : " + rows);
        }
        return rows;
    }

    protected static int checkCols(int cols) {
        if (cols <= 0) {
            throw new IllegalArgumentException("number of columns must be strictly positive : " + cols);
        }
        return cols;
    }

    protected void checkSubmatrixIndices(int r0, int c0, int r1, int c1) {
        checkIJ(r0, c0);
        checkIJ(r1, c1);
        int _rows = r1 - r0 + 1;
        int _cols = c1 - c0 + 1;
        if (_rows <= 0 || _cols <= 0) {
            throw new IllegalArgumentException(
                    "Illegal submatrix indices : [" + r0 + ", " + c0 + ", " + r1 + ", " + c1 + "]");
        }
    }
}
