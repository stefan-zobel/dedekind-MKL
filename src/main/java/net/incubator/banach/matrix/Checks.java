package net.incubator.banach.matrix;

final class Checks {

    static void checkMult(Dimensions A, Dimensions B) {
        if (A.numColumns() != B.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != B.numRows() (" + A.numColumns() + " != " + B.numRows() + ")");
        }
    }

    static void checkTrans(Dimensions A, Dimensions AT) {
        if (A.numRows() != AT.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "A.numRows() != B.numColumns() (" + A.numRows() + " != " + AT.numColumns() + ")");
        }
        if (A.numColumns() != AT.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != B.numRows() (" + A.numColumns() + " != " + AT.numRows() + ")");
        }
    }

    static void checkEqualDimension(Dimensions A, Dimensions B) {
        if (A.numRows() != B.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numRows() != B.numRows() (" + A.numRows() + " != " + B.numRows() + ")");
        }
        if (A.numColumns() != B.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != B.numColumns() (" + A.numColumns() + " != " + B.numColumns() + ")");
        }
    }

    static void checkAdd(Dimensions A, Dimensions B, Dimensions C) {
        if (A.numRows() != B.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numRows() != B.numRows() (" + A.numRows() + " != " + B.numRows() + ")");
        }
        if (A.numColumns() != B.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != B.numColumns() (" + A.numColumns() + " != " + B.numColumns() + ")");
        }
        if (B.numRows() != C.numRows()) {
            throw new IndexOutOfBoundsException(
                    "B.numRows() != C.numRows() (" + B.numRows() + " != " + C.numRows() + ")");
        }
        if (B.numColumns() != C.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "B.numColumns() != C.numColumns() (" + B.numColumns() + " != " + C.numColumns() + ")");
        }
    }

    static void checkMultAdd(Dimensions A, Dimensions B, Dimensions C) {
        if (A.numRows() != C.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numRows() != C.numRows() (" + A.numRows() + " != " + C.numRows() + ")");
        }
        if (A.numColumns() != B.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != B.numRows() (" + A.numColumns() + " != " + B.numRows() + ")");
        }
        if (B.numColumns() != C.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "B.numColumns() != C.numColumns() (" + B.numColumns() + " != " + C.numColumns() + ")");
        }
    }

    static void checkTransABmultAdd(Dimensions A, Dimensions B, Dimensions C) {
        if (A.numRows() != B.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "A.numRows() != B.numColumns() (" + A.numRows() + " != " + B.numColumns() + ")");
        }
        if (A.numColumns() != C.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != C.numRows() (" + A.numColumns() + " != " + C.numRows() + ")");
        }
        if (B.numRows() != C.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "B.numRows() != C.numColumns() (" + B.numRows() + " != " + C.numColumns() + ")");
        }
    }

    static void checkTransAmultAdd(Dimensions A, Dimensions B, Dimensions C) {
        if (A.numRows() != B.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numRows() != B.numRows() (" + A.numRows() + " != " + B.numRows() + ")");
        }
        if (A.numColumns() != C.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != C.numRows() (" + A.numColumns() + " != " + C.numRows() + ")");
        }
        if (B.numColumns() != C.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "B.numColumns() != C.numColumns() (" + B.numColumns() + " != " + C.numColumns() + ")");
        }
    }

    static void checkTransBmultAdd(Dimensions A, Dimensions B, Dimensions C) {
        if (A.numColumns() != B.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != B.numColumns() (" + A.numColumns() + " != " + B.numColumns() + ")");
        }
        if (A.numRows() != C.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numRows() != C.numRows() (" + A.numRows() + " != " + C.numRows() + ")");
        }
        if (B.numRows() != C.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "B.numRows() != C.numColumns() (" + B.numRows() + " != " + C.numColumns() + ")");
        }
    }

    private Checks() {
        throw new AssertionError();
    }
}
