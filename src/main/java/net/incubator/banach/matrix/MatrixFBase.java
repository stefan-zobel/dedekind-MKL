package net.incubator.banach.matrix;

import java.util.Arrays;

/**
 * A {@code MatrixFBase} is a partial implementation of a dense matrix of
 * primitive floats with column-major storage layout. The addressing is zero
 * based. All operations throw a {@code NullPointerException} if any of the
 * method arguments is {@code null}.
 * <p>
 * <b>Note: this is experimental, unfinished and completely untested code!</b>
 */
public abstract class MatrixFBase extends DimensionsBase implements MatrixF {

    protected final float[] a;

    public MatrixFBase(int rows, int cols, float[] array, boolean doArrayCopy) {
        super(rows, cols);
        checkArrayLength(array, rows, cols);
        if (doArrayCopy) {
            float[] copy = new float[array.length];
            System.arraycopy(array, 0, copy, 0, copy.length);
            a = copy;
        } else {
            a = array;
        }
    }

    @Override
    public float toScalar() {
        if (!isScalar()) {
            throw new IllegalStateException("(" + rows + " x " + cols + ") matrix is not a scalar");
        }
        return a[0];
    }

    /**
     * {@code A = alpha * A}
     * 
     * @param alpha
     * @return {@code A}
     */
    @Override
    public MatrixF scaleInplace(float alpha) {
        if (alpha == 0.0f) {
            return this.zeroInplace();
        }
        float[] _a = a;
        for (int i = 0; i < _a.length; ++i) {
            _a[i] *= alpha;
        }
        return this;
    }

    /**
     * {@code B = alpha * A}
     * 
     * @param alpha
     * @param B
     * @return {@code B}
     */
    @Override
    public MatrixF scale(float alpha, MatrixF B) {
        Checks.checkEqualDimension(this, B);
        if (alpha == 0.0f) {
            Arrays.fill(B.getArrayUnsafe(), 0.0f);
            return B;
        }
        float[] _a = a;
        float[] _b = B.getArrayUnsafe();
        for (int i = 0; i < _b.length; ++i) {
            _b[i] = alpha * _a[i];
        }
        return B;
    }

    /**
     * <code>AT = A<sup>T</sup></code>
     * 
     * @param AT
     * @return {@code AT}
     */
    @Override
    public MatrixF trans(MatrixF AT) {
        Checks.checkTrans(this, AT);
        int cols_ = cols;
        int rows_ = rows;
        for (int col = 0; col < cols_; ++col) {
            for (int row = 0; row < rows_; ++row) {
                AT.setUnsafe(col, row, getUnsafe(row, col));
            }
        }
        return AT;
    }

    /**
     * {@code A = A + B}
     * 
     * @param B
     * @return {@code A}
     */
    @Override
    public MatrixF addInplace(MatrixF B) {
        return addInplace(1.0f, B);
    }

    /**
     * {@code A = A + alpha * B}
     * 
     * @param alpha
     * @param B
     * @return {@code A}
     */
    @Override
    public MatrixF addInplace(float alpha, MatrixF B) {
        Checks.checkEqualDimension(this, B);
        if (alpha != 0.0f) {
            int cols_ = cols;
            int rows_ = rows;
            for (int col = 0; col < cols_; ++col) {
                for (int row = 0; row < rows_; ++row) {
                    addUnsafe(row, col, alpha * B.getUnsafe(row, col));
                }
            }
        }
        return this;
    }

    /**
     * {@code C = A + B}
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF add(MatrixF B, MatrixF C) {
        return add(1.0f, B, C);
    }

    /**
     * {@code C = A + alpha * B}
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF add(float alpha, MatrixF B, MatrixF C) {
        Checks.checkAdd(this, B, C);
        if (alpha == 0.0f) {
            System.arraycopy(a, 0, C.getArrayUnsafe(), 0, a.length);
        } else {
            int cols_ = cols;
            int rows_ = rows;
            for (int col = 0; col < cols_; ++col) {
                for (int row = 0; row < rows_; ++row) {
                    C.setUnsafe(row, col, getUnsafe(row, col) + alpha * B.getUnsafe(row, col));
                }
            }
        }
        return C;
    }

    /**
     * {@code C = A * B}
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF mult(MatrixF B, MatrixF C) {
        return mult(1.0f, B, C);
    }

    /**
     * {@code C = alpha * A * B}
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF mult(float alpha, MatrixF B, MatrixF C) {
        return multAdd(alpha, B, C.zeroInplace());
    }

    /**
     * {@code C = A * B + C}
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF multAdd(MatrixF B, MatrixF C) {
        return multAdd(1.0f, B, C);
    }

    /**
     * {@code C = alpha * A * B + C}
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public abstract MatrixF multAdd(float alpha, MatrixF B, MatrixF C);

    /**
     * <code>C = A<sup>T</sup> * B<sup>T</sup></code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF transABmult(MatrixF B, MatrixF C) {
        return transABmult(1.0f, B, C);
    }

    /**
     * <code>C = alpha * A<sup>T</sup> * B<sup>T</sup></code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF transABmult(float alpha, MatrixF B, MatrixF C) {
        return transABmultAdd(alpha, B, C.zeroInplace());
    }

    /**
     * <code>C = A<sup>T</sup> * B</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF transAmult(MatrixF B, MatrixF C) {
        return transAmult(1.0f, B, C);
    }

    /**
     * <code>C = alpha * A<sup>T</sup> * B</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF transAmult(float alpha, MatrixF B, MatrixF C) {
        return transAmultAdd(alpha, B, C.zeroInplace());
    }

    /**
     * <code>C = A * B<sup>T</sup></code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF transBmult(MatrixF B, MatrixF C) {
        return transBmult(1.0f, B, C);
    }

    /**
     * <code>C = alpha * A * B<sup>T</sup></code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF transBmult(float alpha, MatrixF B, MatrixF C) {
        return transBmultAdd(alpha, B, C.zeroInplace());
    }

    /**
     * <code>C = A<sup>T</sup> * B<sup>T</sup> + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF transABmultAdd(MatrixF B, MatrixF C) {
        return transABmultAdd(1.0f, B, C);
    }

    /**
     * <code>C = alpha * A<sup>T</sup> * B<sup>T</sup> + C</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public abstract MatrixF transABmultAdd(float alpha, MatrixF B, MatrixF C);

    /**
     * <code>C = A<sup>T</sup> * B + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF transAmultAdd(MatrixF B, MatrixF C) {
        return transAmultAdd(1.0f, B, C);
    }

    /**
     * <code>C = alpha * A<sup>T</sup> * B + C</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public abstract MatrixF transAmultAdd(float alpha, MatrixF B, MatrixF C);

    /**
     * <code>C = A * B<sup>T</sup> + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixF transBmultAdd(MatrixF B, MatrixF C) {
        return transBmultAdd(1.0f, B, C);
    }

    /**
     * <code>C = alpha * A * B<sup>T</sup> + C</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public abstract MatrixF transBmultAdd(float alpha, MatrixF B, MatrixF C);

    @Override
    public MatrixF zeroInplace() {
        Arrays.fill(a, 0.0f);
        return this;
    }

    @Override
    public float get(int row, int col) {
        checkIJ(row, col);
        return a[idx(row, col)];
    }

    public float getUnsafe(int row, int col) {
        return a[idx(row, col)];
    }

    @Override
    public MatrixF set(int row, int col, float val) {
        checkIJ(row, col);
        a[idx(row, col)] = val;
        return this;
    }

    public void setUnsafe(int row, int col, float val) {
        a[idx(row, col)] = val;
    }

    @Override
    public MatrixF add(int row, int col, float val) {
        checkIJ(row, col);
        a[idx(row, col)] += val;
        return this;
    }

    protected void addUnsafe(int row, int col, float val) {
        a[idx(row, col)] += val;
    }

    @Override
    public float[] getArrayUnsafe() {
        return a;
    }

    protected static void checkArrayLength(float[] array, int rows, int cols) {
        if (array.length != rows * cols) {
            throw new IllegalArgumentException(
                    "data array has wrong length. Needed : " + rows * cols + " , Is : " + array.length);
        }
    }
}