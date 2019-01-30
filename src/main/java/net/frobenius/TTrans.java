package net.frobenius;

public enum TTrans {

	NO_TRANS("N"), //
	TRANS("T"), //
	CONJ_TRANS("C"); //

	private final String val;

	private TTrans(String s) {
		val = s;
	}

	public String val() {
		return val;
	}
}
