package painting;


/**
 * Generic Class to hold max value and position in an array
 * @author Chris Gamble, DPhil Candidate in
 *  Statistical Genetics, University Of Oxford,
 *  Copyright 2012.
 *
 */
public class Max {
	private int maxIndex;
	private double max;

	public Max (int maxIndex, double max) {
		this.max = max;
		this.maxIndex = maxIndex;
	}

	public int getMaxIndex() {
		return maxIndex;
	}

	public double getMax() {
		return max;
	}
}

