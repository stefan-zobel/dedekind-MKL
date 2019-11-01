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

    private Matrices() {
        throw new AssertionError();
    }
}
