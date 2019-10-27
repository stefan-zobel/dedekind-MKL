package net.incubator.banach.matrix;

import java.util.Arrays;

import net.dedekind.blas.Blas;

/**
 * <b>Note: this is experimental, unfinished and completely untested code!</b> 
 */
public class SimpleMatrixF extends MatrixFBase implements MatrixF {

    private static final float BETA = 1.0f;

    public SimpleMatrixF(int rows, int cols) {
        super(rows, cols, new float[rows * cols], false);
    }

    public SimpleMatrixF(int rows, int cols, float initialValue) {
        super(rows, cols, new float[rows * cols], false);
        Arrays.fill(a, initialValue);
    }

    @Override
    public MatrixF multAdd(float alpha, MatrixF B, MatrixF C) {
        // TODO Auto-generated method stub
        Checks.checkMultAdd(this, B, C);

        Blas blas = Blas.getInstance();
        throw new UnsupportedOperationException();
    }

    @Override
    public MatrixF transABmultAdd(float alpha, MatrixF B, MatrixF C) {
        // TODO Auto-generated method stub
        Checks.checkTransABmultAdd(this, B, C);

        Blas blas = Blas.getInstance();
        throw new UnsupportedOperationException();
    }

    @Override
    public MatrixF transAmultAdd(float alpha, MatrixF B, MatrixF C) {
        // TODO Auto-generated method stub
        Checks.checkTransAmultAdd(this, B, C);

        Blas blas = Blas.getInstance();
        throw new UnsupportedOperationException();
    }

    @Override
    public MatrixF transBmultAdd(float alpha, MatrixF B, MatrixF C) {
        // TODO Auto-generated method stub
        Checks.checkTransBmultAdd(this, B, C);

        Blas blas = Blas.getInstance();
        throw new UnsupportedOperationException();
    }
}
