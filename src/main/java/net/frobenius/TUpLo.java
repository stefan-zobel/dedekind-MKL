package net.frobenius;

public enum TUpLo {

    UPPER("U"), //
    LOWER("L"); //

    private final String val;

    private TUpLo(String s) {
        val = s;
    }

    public String val() {
        return val;
    }
}
