package net.frobenius;

@SuppressWarnings("serial")
public class ComputationTruncatedException extends RuntimeException {

	public ComputationTruncatedException() {
	}

	public ComputationTruncatedException(String message) {
		super(message);
	}

	public ComputationTruncatedException(Throwable cause) {
		super(cause);
	}

	public ComputationTruncatedException(String message, Throwable cause) {
		super(message, cause);
	}

	public ComputationTruncatedException(String message, Throwable cause, boolean enableSuppression,
			boolean writableStackTrace) {
		super(message, cause, enableSuppression, writableStackTrace);
	}
}
