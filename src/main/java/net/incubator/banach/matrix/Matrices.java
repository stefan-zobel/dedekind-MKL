package net.incubator.banach.matrix;

import java.util.concurrent.ThreadLocalRandom;

/**
 * Static utility methods for matrices.
 */
public final class Matrices {

    public static MatrixD identityD(int n) {
        SimpleMatrixD m = new SimpleMatrixD(n, n);
        for (int i = 0; i < n; ++i) {
            m.set(i, i, 1.0);
        }
        return m;
    }

    public static MatrixF identityF(int n) {
        SimpleMatrixF m = new SimpleMatrixF(n, n);
        for (int i = 0; i < n; ++i) {
            m.set(i, i, 1.0f);
        }
        return m;
    }

    public static MatrixD randomUniformD(int rows, int cols) {
        SimpleMatrixD m = new SimpleMatrixD(rows, cols);
        ThreadLocalRandom rnd = ThreadLocalRandom.current();
        double[] _a = m.getArrayUnsafe();
        for (int i = 0; i < _a.length; ++i) {
            _a[i] = rnd.nextDouble();
        }
        return m;
    }

    public static MatrixF randomUniformF(int rows, int cols) {
        SimpleMatrixF m = new SimpleMatrixF(rows, cols);
        ThreadLocalRandom rnd = ThreadLocalRandom.current();
        float[] _a = m.getArrayUnsafe();
        for (int i = 0; i < _a.length; ++i) {
            _a[i] = rnd.nextFloat();
        }
        return m;
    }

    public static MatrixD randomNormalD(int rows, int cols) {
        SimpleMatrixD m = new SimpleMatrixD(rows, cols);
        ThreadLocalRandom rnd = ThreadLocalRandom.current();
        double[] _a = m.getArrayUnsafe();
        for (int i = 0; i < _a.length; ++i) {
            _a[i] = rnd.nextGaussian();
        }
        return m;
    }

    public static MatrixF randomNormalF(int rows, int cols) {
        SimpleMatrixF m = new SimpleMatrixF(rows, cols);
        ThreadLocalRandom rnd = ThreadLocalRandom.current();
        float[] _a = m.getArrayUnsafe();
        for (int i = 0; i < _a.length; ++i) {
            _a[i] = (float) rnd.nextGaussian();
        }
        return m;
    }

    /**
     * Useful for tests.
     * 
     * @param rows
     * @param cols
     * @return matrix filled with the natural numbers starting with 1 in
     *         row-major order
     */
    public static MatrixD naturalNumbersD(int rows, int cols) {
        SimpleMatrixD m = new SimpleMatrixD(rows, cols);
        int nat = 1;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                m.set(i, j, nat++);
            }
        }
        return m;
    }

    /**
     * Useful for tests.
     * 
     * @param rows
     * @param cols
     * @return matrix filled with the natural numbers starting with 1 in
     *         row-major order
     */
    public static MatrixF naturalNumbersF(int rows, int cols) {
        SimpleMatrixF m = new SimpleMatrixF(rows, cols);
        int nat = 1;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                m.set(i, j, nat++);
            }
        }
        return m;
    }

    private Matrices() {
        throw new AssertionError();
    }
}
