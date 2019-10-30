package net.incubator.banach.matrix;

/**
 * A {@code MatrixF} is a dense matrix of primitive floats with column-major
 * storage layout. The addressing is zero based. All operations throw a
 * {@code NullPointerException} if any of the method arguments is {@code null}.
 */
public interface MatrixF extends Dimensions {

    float toScalar();

    /**
     * {@code A = alpha * A}
     * 
     * @param alpha
     * @return {@code A}
     */
    MatrixF scaleInplace(float alpha);

    /**
     * {@code B = alpha * A}
     * 
     * @param alpha
     * @param B
     * @return {@code B}
     */
    MatrixF scale(float alpha, MatrixF B);

    /**
     * <code>AT = A<sup>T</sup></code>
     * 
     * @param AT
     * @return {@code AT}
     */
    MatrixF trans(MatrixF AT);

    /**
     * {@code A = A + B}
     * 
     * @param B
     * @return {@code A}
     */
    MatrixF addInplace(MatrixF B);

    /**
     * {@code A = A + alpha * B}
     * 
     * @param alpha
     * @param B
     * @return {@code A}
     */
    MatrixF addInplace(float alpha, MatrixF B);

    /**
     * {@code C = A + B}
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF add(MatrixF B, MatrixF C);

    /**
     * {@code C = A + alpha * B}
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF add(float alpha, MatrixF B, MatrixF C);

    /**
     * {@code C = A * B}
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF mult(MatrixF B, MatrixF C);

    /**
     * {@code C = alpha * A * B}
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF mult(float alpha, MatrixF B, MatrixF C);

    /**
     * {@code C = A * B + C}
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF multAdd(MatrixF B, MatrixF C);

    /**
     * {@code C = alpha * A * B + C}
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF multAdd(float alpha, MatrixF B, MatrixF C);

    /**
     * <code>C = A<sup>T</sup> * B<sup>T</sup></code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transABmult(MatrixF B, MatrixF C);

    /**
     * <code>C = alpha * A<sup>T</sup> * B<sup>T</sup></code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transABmult(float alpha, MatrixF B, MatrixF C);

    /**
     * <code>C = A<sup>T</sup> * B</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transAmult(MatrixF B, MatrixF C);

    /**
     * <code>C = alpha * A<sup>T</sup> * B</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transAmult(float alpha, MatrixF B, MatrixF C);

    /**
     * <code>C = A * B<sup>T</sup></code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transBmult(MatrixF B, MatrixF C);

    /**
     * <code>C = alpha * A * B<sup>T</sup></code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transBmult(float alpha, MatrixF B, MatrixF C);

    /**
     * <code>C = A<sup>T</sup> * B<sup>T</sup> + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transABmultAdd(MatrixF B, MatrixF C);

    /**
     * <code>C = alpha * A<sup>T</sup> * B<sup>T</sup> + C</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transABmultAdd(float alpha, MatrixF B, MatrixF C);

    /**
     * <code>C = A<sup>T</sup> * B + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transAmultAdd(MatrixF B, MatrixF C);

    /**
     * <code>C = alpha * A<sup>T</sup> * B + C</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transAmultAdd(float alpha, MatrixF B, MatrixF C);

    /**
     * <code>C = A * B<sup>T</sup> + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transBmultAdd(MatrixF B, MatrixF C);

    /**
     * <code>C = alpha * A * B<sup>T</sup> + C</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixF transBmultAdd(float alpha, MatrixF B, MatrixF C);

    MatrixF copy();

    MatrixF zeroInplace();

    MatrixF setInplace(MatrixF other);

    /**
     * {@code A = alpha * B}
     * @param alpha
     * @param B
     * @return {@code A}
     */
    MatrixF setInplace(float alpha, MatrixF B);

    float get(int row, int col);

    MatrixF set(int row, int col, float val);

    MatrixF add(int row, int col, float val);

    float[] getArrayUnsafe();

    float getUnsafe(int row, int col);

    void setUnsafe(int row, int col, float val);
}
