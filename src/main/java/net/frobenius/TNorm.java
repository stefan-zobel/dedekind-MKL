package net.frobenius;

public enum TNorm {

	ONE("1"), //
	INF("I"); //

	private final String val;

	private TNorm(String s) {
		val = s;
	}

	public String val() {
		return val;
	}
}
