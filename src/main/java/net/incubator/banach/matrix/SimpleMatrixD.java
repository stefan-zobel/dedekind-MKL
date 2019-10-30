package net.incubator.banach.matrix;

import java.util.Arrays;

import net.dedekind.blas.Blas;
import net.frobenius.TTrans;

/**
 * <b>Note: this is experimental, unfinished and completely untested code!</b>
 */
public class SimpleMatrixD extends MatrixDBase implements MatrixD {

    private static final double BETA = 1.0;

    public SimpleMatrixD(int rows, int cols) {
        super(rows, cols, new double[rows * cols], false);
    }

    public SimpleMatrixD(int rows, int cols, double initialValue) {
        super(rows, cols, new double[rows * cols], false);
        Arrays.fill(a, initialValue);
    }

    private SimpleMatrixD(SimpleMatrixD other) {
        super(other.rows, other.cols, other.a, true);
    }

    @Override
    public MatrixD multAdd(double alpha, MatrixD B, MatrixD C) {
        Checks.checkMultAdd(this, B, C);

        Blas blas = Blas.getInstance();
        blas.dgemm(TTrans.NO_TRANS.val(), TTrans.NO_TRANS.val(), C.numRows(), C.numColumns(), cols, alpha, a,
                Math.max(1, rows), B.getArrayUnsafe(), Math.max(1, B.numRows()), BETA, C.getArrayUnsafe(),
                Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public MatrixD transABmultAdd(double alpha, MatrixD B, MatrixD C) {
        Checks.checkTransABmultAdd(this, B, C);

        Blas blas = Blas.getInstance();
        blas.dgemm(TTrans.TRANS.val(), TTrans.TRANS.val(), C.numRows(), C.numColumns(), rows, alpha, a,
                Math.max(1, rows), B.getArrayUnsafe(), Math.max(1, B.numRows()), BETA, C.getArrayUnsafe(),
                Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public MatrixD transAmultAdd(double alpha, MatrixD B, MatrixD C) {
        Checks.checkTransAmultAdd(this, B, C);

        Blas blas = Blas.getInstance();
        blas.dgemm(TTrans.TRANS.val(), TTrans.NO_TRANS.val(), C.numRows(), C.numColumns(), rows, alpha, a,
                Math.max(1, rows), B.getArrayUnsafe(), Math.max(1, B.numRows()), BETA, C.getArrayUnsafe(),
                Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public MatrixD transBmultAdd(double alpha, MatrixD B, MatrixD C) {
        Checks.checkTransBmultAdd(this, B, C);

        Blas blas = Blas.getInstance();
        blas.dgemm(TTrans.NO_TRANS.val(), TTrans.TRANS.val(), C.numRows(), C.numColumns(), cols, alpha, a,
                Math.max(1, rows), B.getArrayUnsafe(), Math.max(1, B.numRows()), BETA, C.getArrayUnsafe(),
                Math.max(1, C.numRows()));

        return C;
    }

    @Override
    public MatrixD copy() {
        return new SimpleMatrixD(this); 
    }
}
