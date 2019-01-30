package net.frobenius;

public enum TSvdJob {

    ALL("A"), //
    PART("S"), //
    OVERWRITE("O"), //
    NONE("N"); //

    private final String val;

    private TSvdJob(String s) {
        val = s;
    }

    public String val() {
        return val;
    }
}
