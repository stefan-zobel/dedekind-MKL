package net.incubator.banach.matrix;

public interface Dimensions {

    boolean isScalar();

    boolean isColumnVector();

    boolean isRowVector();

    boolean isSquareMatrix();

    int numColumns();

    int numRows();

    void checkIndex(int i, int j);

    /**
     * {@code (rFrom, cFrom)} upper left corner, {@code (rTo, cTo)} lower right
     * corner. All indexes must be valid and the submatrix must select at least
     * one element.
     * 
     * @param rFrom
     *            {@code <= rTo}
     * @param cFrom
     *            {@code <= cTo}
     * @param rTo
     *            {@code >= rFrom}
     * @param cTo
     *            {@code >= rTo}
     */
    void checkSubmatrixIndexes(int rFrom, int cFrom, int rTo, int cTo);
}
