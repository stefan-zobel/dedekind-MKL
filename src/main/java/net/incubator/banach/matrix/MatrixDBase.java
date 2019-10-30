package net.incubator.banach.matrix;

import java.util.Arrays;

/**
 * A {@code MatrixDBase} is a partial implementation of a dense matrix of
 * primitive doubles with column-major storage layout. The addressing is zero
 * based. All operations throw a {@code NullPointerException} if any of the
 * method arguments is {@code null}.
 * <p>
 * <b>Note: this is experimental, unfinished and completely untested code!</b>
 */
public abstract class MatrixDBase extends DimensionsBase implements MatrixD {

    protected final double[] a;

    public MatrixDBase(int rows, int cols, double[] array, boolean doArrayCopy) {
        super(rows, cols);
        checkArrayLength(array, rows, cols);
        if (doArrayCopy) {
            double[] copy = new double[array.length];
            System.arraycopy(array, 0, copy, 0, copy.length);
            a = copy;
        } else {
            a = array;
        }
    }

    @Override
    public double toScalar() {
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
    public MatrixD scaleInplace(double alpha) {
        if (alpha == 0.0) {
            return this.zeroInplace();
        }
        if (alpha == 1.0) {
            return this;
        }
        double[] _a = a;
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
    public MatrixD scale(double alpha, MatrixD B) {
        Checks.checkEqualDimension(this, B);
        if (alpha == 0.0) {
            Arrays.fill(B.getArrayUnsafe(), 0.0);
            return B;
        }
        double[] _a = a;
        double[] _b = B.getArrayUnsafe();
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
    public MatrixD trans(MatrixD AT) {
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
    public MatrixD addInplace(MatrixD B) {
        return addInplace(1.0, B);
    }

    /**
     * {@code A = A + alpha * B}
     * 
     * @param alpha
     * @param B
     * @return {@code A}
     */
    @Override
    public MatrixD addInplace(double alpha, MatrixD B) {
        Checks.checkEqualDimension(this, B);
        if (alpha != 0.0) {
            double[] _a = a;
            double[] _b = B.getArrayUnsafe();
            for (int i = 0; i < _b.length; ++i) {
                _a[i] += alpha * _b[i];
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
    public MatrixD add(MatrixD B, MatrixD C) {
        return add(1.0, B, C);
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
    public MatrixD add(double alpha, MatrixD B, MatrixD C) {
        Checks.checkAdd(this, B, C);
        if (alpha == 0.0) {
            System.arraycopy(a, 0, C.getArrayUnsafe(), 0, a.length);
        } else {
            double[] _a = a;
            double[] _b = B.getArrayUnsafe();
            double[] _c = C.getArrayUnsafe();
            for (int i = 0; i < _a.length; ++i) {
                _c[i] = _a[i] + alpha * _b[i];
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
    public MatrixD mult(MatrixD B, MatrixD C) {
        return mult(1.0, B, C);
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
    public MatrixD mult(double alpha, MatrixD B, MatrixD C) {
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
    public MatrixD multAdd(MatrixD B, MatrixD C) {
        return multAdd(1.0, B, C);
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
    public abstract MatrixD multAdd(double alpha, MatrixD B, MatrixD C);

    /**
     * <code>C = A<sup>T</sup> * B<sup>T</sup></code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixD transABmult(MatrixD B, MatrixD C) {
        return transABmult(1.0, B, C);
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
    public MatrixD transABmult(double alpha, MatrixD B, MatrixD C) {
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
    public MatrixD transAmult(MatrixD B, MatrixD C) {
        return transAmult(1.0, B, C);
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
    public MatrixD transAmult(double alpha, MatrixD B, MatrixD C) {
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
    public MatrixD transBmult(MatrixD B, MatrixD C) {
        return transBmult(1.0, B, C);
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
    public MatrixD transBmult(double alpha, MatrixD B, MatrixD C) {
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
    public MatrixD transABmultAdd(MatrixD B, MatrixD C) {
        return transABmultAdd(1.0, B, C);
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
    public abstract MatrixD transABmultAdd(double alpha, MatrixD B, MatrixD C);

    /**
     * <code>C = A<sup>T</sup> * B + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixD transAmultAdd(MatrixD B, MatrixD C) {
        return transAmultAdd(1.0, B, C);
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
    public abstract MatrixD transAmultAdd(double alpha, MatrixD B, MatrixD C);

    /**
     * <code>C = A * B<sup>T</sup> + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    @Override
    public MatrixD transBmultAdd(MatrixD B, MatrixD C) {
        return transBmultAdd(1.0, B, C);
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
    public abstract MatrixD transBmultAdd(double alpha, MatrixD B, MatrixD C);

    @Override
    public MatrixD zeroInplace() {
        Arrays.fill(a, 0.0);
        return this;
    }

    @Override
    public MatrixD setInplace(MatrixD other) {
        return setInplace(1.0, other);
    }

    @Override
    public MatrixD setInplace(double alpha, MatrixD other) {
        Checks.checkEqualDimension(this, other);
        if (alpha == 0.0) {
            return zeroInplace();
        }
        if (other == this) {
            return scaleInplace(alpha);
        }
        double[] _a = a;
        double[] _b = other.getArrayUnsafe();
        for (int i = 0; i < _b.length; ++i) {
            _a[i] = alpha * _b[i];
        }
        return this; 
    }

    @Override
    public double get(int row, int col) {
        checkIJ(row, col);
        return a[idx(row, col)];
    }

    public double getUnsafe(int row, int col) {
        return a[idx(row, col)];
    }

    @Override
    public MatrixD set(int row, int col, double val) {
        checkIJ(row, col);
        a[idx(row, col)] = val;
        return this;
    }

    public void setUnsafe(int row, int col, double val) {
        a[idx(row, col)] = val;
    }

    @Override
    public MatrixD add(int row, int col, double val) {
        checkIJ(row, col);
        a[idx(row, col)] += val;
        return this;
    }

    protected void addUnsafe(int row, int col, double val) {
        a[idx(row, col)] += val;
    }

    @Override
    public double[] getArrayUnsafe() {
        return a;
    }

    @Override
    public double normF() {
        // overflow resistant implementation
        double scale = 0.0;
        double sumsquared = 1.0;
        double[] _a = a;
        for (int i = 0; i < _a.length; ++i) {
            double xi = _a[i];
            if (xi != 0.0) {
                double absxi = Math.abs(xi);
                if (scale < absxi) {
                    double unsquared = scale / absxi;
                    sumsquared = 1.0 + sumsquared * (unsquared * unsquared);
                    scale = absxi;
                } else {
                    double unsquared = absxi / scale;
                    sumsquared = sumsquared + (unsquared * unsquared);
                }
            }
        }
        return scale * Math.sqrt(sumsquared);
    }

    protected static void checkArrayLength(double[] array, int rows, int cols) {
        if (array.length != rows * cols) {
            throw new IllegalArgumentException(
                    "data array has wrong length. Needed : " + rows * cols + " , Is : " + array.length);
        }
    }
}
