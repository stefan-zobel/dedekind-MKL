package net.frobenius;

public enum TSide {

	LEFT("L"), //
	RIGHT("R"); //

	private final String val;

	private TSide(String s) {
		val = s;
	}

	public String val() {
		return val;
	}
}
