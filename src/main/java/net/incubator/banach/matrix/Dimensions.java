package net.incubator.banach.matrix;

public interface Dimensions {

    boolean isScalar();

    boolean isColumnVector();

    boolean isRowVector();

    boolean isSquareMatrix();

    int numColumns();

    int numRows();

    void checkIndex(int i, int j);
}
