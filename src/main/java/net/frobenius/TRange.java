package net.frobenius;

public enum TRange {

	ALL("A"), //
	INDICES("I"), //
	INTERVAL("V"); //

	private final String val;

	private TRange(String s) {
		val = s;
	}

	public String val() {
		return val;
	}
}
