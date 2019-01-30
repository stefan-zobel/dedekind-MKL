package net.frobenius;

public enum TDiag {

	UNIT("U"), //
	NON_UNIT("N"); //

	private final String val;

	private TDiag(String s) {
		val = s;
	}

	public String val() {
		return val;
	}
}
