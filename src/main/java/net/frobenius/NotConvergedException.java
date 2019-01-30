package net.frobenius;

@SuppressWarnings("serial")
public class NotConvergedException extends RuntimeException {

	public NotConvergedException() {
	}

	public NotConvergedException(String message) {
		super(message);
	}

	public NotConvergedException(Throwable cause) {
		super(cause);
	}

	public NotConvergedException(String message, Throwable cause) {
		super(message, cause);
	}

	public NotConvergedException(String message, Throwable cause, boolean enableSuppression,
			boolean writableStackTrace) {
		super(message, cause, enableSuppression, writableStackTrace);
	}
}
