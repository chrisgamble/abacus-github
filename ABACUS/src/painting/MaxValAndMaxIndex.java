package painting;
import java.util.ArrayList;



public class MaxValAndMaxIndex {
	private double[] maxValues;
	private int[] maxIndexes;
	
	public MaxValAndMaxIndex (double[] maxValues, int[] maxIndexes) {
		this.maxIndexes = maxIndexes;
		this.maxValues = maxValues;
	}
	
	public double[] getMaxValues() {
		return maxValues;
	}
	
	public int[] getMaxIndexes() {
		return maxIndexes;
	}
	
	public static MaxValAndMaxIndex create (double[] logProbability, double probabilityOfTransition,
			double probabilityOfNotTransitioning, Max referenceMaximum, double[] localEmmision,
			ArrayList<Integer> haplotypeIndexes, ArrayList<Integer> notHaplotypeIndexes) {
		
		double currentMaximum = referenceMaximum.getMax() + Math.log(probabilityOfTransition);
		
		double[] maxVal = new double[logProbability.length];
		int[] maxIndex = new int[logProbability.length];
		if (logProbability.length > 0) {
			for (int i : haplotypeIndexes) {
				double maxValPreLogic = logProbability[i] + 
						Math.log(probabilityOfNotTransitioning) - currentMaximum;
				if (maxValPreLogic > 0) {
					maxVal[i] = logProbability[i] + 
							Math.log(probabilityOfNotTransitioning) +
							localEmmision[i];
					maxIndex[i] = i + 1;
				} else {
					maxVal[i] = currentMaximum + localEmmision[i];
					maxIndex[i] = referenceMaximum.getMaxIndex();
				}
			}
			for (int i : notHaplotypeIndexes) {
				maxVal[i] = Double.NaN;
				maxIndex[i] = 0;
			}
		}
		return new MaxValAndMaxIndex(maxVal, maxIndex);
	}

}
