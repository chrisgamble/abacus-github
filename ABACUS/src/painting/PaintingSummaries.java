package painting;

/**
 * Class to hold the values of the chromosome wide
 *  painting algorithm.
 * @author chrisgamble
 *
 */
public class PaintingSummaries {
	private float[] totalPairLengths;
	private float[] numberPairLengths;
	
	public PaintingSummaries (float[] totalPairLengths, float[] numberPairLengths) {
		this.totalPairLengths = totalPairLengths;
		this.numberPairLengths = numberPairLengths;
	}
	
	public float[] getTotal() {
		return totalPairLengths;
	}
	
	public float[] getNumber() {
		return numberPairLengths;
	}
	
}
