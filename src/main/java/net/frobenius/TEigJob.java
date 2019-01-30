package net.frobenius;

public enum TEigJob {

	ALL("V"), //
	VALUES_ONLY("N"); //

	private final String val;

	private TEigJob(String s) {
		val = s;
	}

	public String val() {
		return val;
	}
}
