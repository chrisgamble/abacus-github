package bayesfactor;

public class LaplacePair {
	private float[][] array;
	private double[] genomeWideMean;

	public LaplacePair (float[][] array, double[] genomeWideMean) {
		this.array = array;
		this.genomeWideMean = genomeWideMean;
	}

	public float[][] getArray () {
		return array;
	}

	public double[] getGenomeWideMean () {
		return genomeWideMean;
	}

}
