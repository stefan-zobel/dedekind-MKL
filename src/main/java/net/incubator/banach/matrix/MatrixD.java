package net.incubator.banach.matrix;

/**
 * A {@code MatrixD} is a dense matrix of primitive doubles with column-major
 * storage layout. The addressing is zero based. All operations throw a
 * {@code NullPointerException} if any of the method arguments is {@code null}.
 */
public interface MatrixD extends Dimensions {

    double toScalar();

    /**
     * {@code A = alpha * A}
     * 
     * @param alpha
     * @return {@code A}
     */
    MatrixD scaleInplace(double alpha);

    /**
     * {@code B = alpha * A}
     * 
     * @param alpha
     * @param B
     * @return {@code B}
     */
    MatrixD scale(double alpha, MatrixD B);

    /**
     * <code>AT = A<sup>T</sup></code>
     * 
     * @param AT
     * @return {@code AT}
     */
    MatrixD trans(MatrixD AT);

    /**
     * {@code A = A + B}
     * 
     * @param B
     * @return {@code A}
     */
    MatrixD addInplace(MatrixD B);

    /**
     * {@code A = A + alpha * B}
     * 
     * @param alpha
     * @param B
     * @return {@code A}
     */
    MatrixD addInplace(double alpha, MatrixD B);

    /**
     * {@code C = A + B}
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD add(MatrixD B, MatrixD C);

    /**
     * {@code C = A + alpha * B}
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD add(double alpha, MatrixD B, MatrixD C);

    /**
     * {@code C = A * B}
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD mult(MatrixD B, MatrixD C);

    /**
     * {@code C = alpha * A * B}
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD mult(double alpha, MatrixD B, MatrixD C);

    /**
     * {@code C = A * B + C}
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD multAdd(MatrixD B, MatrixD C);

    /**
     * {@code C = alpha * A * B + C}
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD multAdd(double alpha, MatrixD B, MatrixD C);

    /**
     * <code>C = A<sup>T</sup> * B<sup>T</sup></code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transABmult(MatrixD B, MatrixD C);

    /**
     * <code>C = alpha * A<sup>T</sup> * B<sup>T</sup></code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transABmult(double alpha, MatrixD B, MatrixD C);

    /**
     * <code>C = A<sup>T</sup> * B</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transAmult(MatrixD B, MatrixD C);

    /**
     * <code>C = alpha * A<sup>T</sup> * B</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transAmult(double alpha, MatrixD B, MatrixD C);

    /**
     * <code>C = A * B<sup>T</sup></code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transBmult(MatrixD B, MatrixD C);

    /**
     * <code>C = alpha * A * B<sup>T</sup></code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transBmult(double alpha, MatrixD B, MatrixD C);

    /**
     * <code>C = A<sup>T</sup> * B<sup>T</sup> + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transABmultAdd(MatrixD B, MatrixD C);

    /**
     * <code>C = alpha * A<sup>T</sup> * B<sup>T</sup> + C</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transABmultAdd(double alpha, MatrixD B, MatrixD C);

    /**
     * <code>C = A<sup>T</sup> * B + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transAmultAdd(MatrixD B, MatrixD C);

    /**
     * <code>C = alpha * A<sup>T</sup> * B + C</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transAmultAdd(double alpha, MatrixD B, MatrixD C);

    /**
     * <code>C = A * B<sup>T</sup> + C</code>
     * 
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transBmultAdd(MatrixD B, MatrixD C);

    /**
     * <code>C = alpha * A * B<sup>T</sup> + C</code>
     * 
     * @param alpha
     * @param B
     * @param C
     * @return {@code C}
     */
    MatrixD transBmultAdd(double alpha, MatrixD B, MatrixD C);

    MatrixD copy();

    MatrixD zeroInplace();

    MatrixD setInplace(MatrixD other);

    /**
     * {@code A = alpha * B}
     * @param alpha
     * @param B
     * @return {@code A}
     */
    MatrixD setInplace(double alpha, MatrixD other);

    double get(int row, int col);

    MatrixD set(int row, int col, double val);

    MatrixD add(int row, int col, double val);

    double[] getArrayUnsafe();

    double getUnsafe(int row, int col);

    void setUnsafe(int row, int col, double val);
}
