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
}
