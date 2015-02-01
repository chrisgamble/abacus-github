package painting;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * This is the main functional class which does much of the
 * heavy lifting to determine the Viterbi Paths. Depending on
 * performance, a simple multi-threaded application could be written.
 * @author Chris Gamble, DPhil Candidate in
 *  Statistical Genetics, University Of Oxford,
 *  Copyright 2012.
 *
 */
public class PaintingDonorAlgorithm implements Callable<PaintingSummaries> {
	private int haplotypeIndex;
	private ArrayList<Double> recombinationMap;
	private ArrayList<Integer> snpMap;
	private ArrayList<String> rsid;
	private ArrayList<Integer> donor;
	private Byte[][] haplotype;
	private double effectivePopulation;
	private double mutationParameter;
	private int numberOfSnps;
	
	private String prefix;

	/**
	 * Class Constructor
	 * @param logProbability
	 * @param pathTraceArray
	 */
	public PaintingDonorAlgorithm (int haplotypeIndex, ArrayList<Double> recommbinationMap,
			ArrayList<Integer> snpMap, ArrayList<String> rsid,
			ArrayList<Integer> donor, Byte[][] haplotype, double effectivePopulation, double mutationParameter,
			int numberOfSnps, String prefix) {
		this.haplotypeIndex = haplotypeIndex;
		this.recombinationMap = recommbinationMap;
		this.snpMap = snpMap;
		this.rsid = rsid;
		this.donor = donor;
		this.haplotype = haplotype;
		this.effectivePopulation = effectivePopulation;
		this.mutationParameter = mutationParameter;
		this.numberOfSnps = numberOfSnps;
		this.prefix = prefix;

	}

	public PaintingSummaries call () throws Exception{
		double sizeOfStateSpace = (double) donor.size();
		byte[] haplotypeRow = new byte[numberOfSnps];
		ArrayList<Integer> haplotypeIndexes = new ArrayList<Integer>();
		ArrayList<Integer> notHaplotypeIndexes = new ArrayList<Integer>();
		for (int i = 0; i < donor.size(); i++) haplotypeIndexes.add(i);

		haplotypeRow[0] = haplotype[0][haplotypeIndex];
		Byte[] haplotypeColumn = new Byte[donor.size()];
		for (int i = 0; i < donor.size(); i++) haplotypeColumn[i] = haplotype[0][donor.get(i)]; 
		double[] logProbability = getEmmisionValues(haplotypeRow[0], haplotypeColumn,
				mutationParameter, haplotypeIndexes);
		for (int i = 0; i < donor.size(); i++) {
			logProbability[i] = logProbability[i] - Math.log(sizeOfStateSpace);
		}
		
		short[][] pathTraceArray = new short[numberOfSnps][donor.size()];
		for (int snp = 1; snp < numberOfSnps; snp++) {
			
			haplotypeRow[snp] = haplotype[snp][haplotypeIndex];
			haplotypeColumn = new Byte[donor.size()];
			for (int i = 0; i < donor.size(); i++) haplotypeColumn[i] = haplotype[0][donor.get(i)]; 
			
			double localRecombination = (recombinationMap.get(snp) -
					recombinationMap.get(snp - 1)) / 100.0;
			double p = -1 * localRecombination * 4 * effectivePopulation / sizeOfStateSpace;
			double probabilityOfTransition = (1 - Math.exp(p)) / sizeOfStateSpace;
			double probabilityOfNotTransitioning = probabilityOfTransition + Math.exp(p);
			Max referenceMaximum = oneShotMax(logProbability, haplotypeIndexes);
			double[] localEmmision = getEmmisionValues(haplotypeRow[snp], haplotypeColumn,
					mutationParameter, haplotypeIndexes);
			MaxValAndMaxIndex betterName = MaxValAndMaxIndex.create(logProbability,
					probabilityOfTransition, probabilityOfNotTransitioning,
					referenceMaximum, localEmmision, haplotypeIndexes, notHaplotypeIndexes);
			logProbability = betterName.getMaxValues();
			int[] next = betterName.getMaxIndexes();
			for (int i = 0; i < next.length; i++) {
				pathTraceArray[snp - 1][i] = (short) next[i];
			}
		}
		short[] vPath = new short[numberOfSnps];
		vPath[numberOfSnps - 1] = (short) oneShotMax(logProbability, haplotypeIndexes).getMaxIndex();
		for (int l = numberOfSnps - 1; l > 0; l--) {
			vPath[l - 1] = pathTraceArray[l - 1][vPath[l] - 1];
		}
		pathTraceArray = null;
		BufferedWriter tmpOut = new BufferedWriter(new  FileWriter(prefix + "_" + haplotypeIndex + ".viterbi"));
		LengthsFromPainting output = getLengths(vPath);
		for (int i = 0; i < numberOfSnps; i++) {
			tmpOut.write(snpMap.get(i) + " " + rsid.get(i) + " " + vPath[i] + " " + Double.toString(output.getBlocks()[i]));
			tmpOut.newLine();
		}
		tmpOut.close();
		System.out.println("Viterbi " + haplotypeIndex + " complete!");
		return new PaintingSummaries(output.getLengths(), output.getCounts());
	}

	/**
	 * Determines the lengths of each contiguous block using the painted
	 * path (vPath).
	 * @param vPath {@code short[]}
	 * @return {@link LengthsFromPainting}
	 * 
	 */
	public LengthsFromPainting getLengths (short[] vPath) throws IOException {
		float[] totalPairLengths = new float[donor.size()];
		float[] numberPairLengths = new float[donor.size()];
		float[] blockLength = new float[numberOfSnps];
		int beginningOfCurrentBlock = 0;
		for (int i = 0; i < numberOfSnps - 1; i++) {
			if (vPath[i] != vPath[i + 1]) {
				double recombination = recombinationMap.get(i) - recombinationMap.get(beginningOfCurrentBlock);
				for (int j = beginningOfCurrentBlock; j <= i; j++) {
					blockLength[j] = (float) recombination;
				}
				if (i == beginningOfCurrentBlock) {
					recombination = (recombinationMap.get(i + 1) - recombinationMap.get(beginningOfCurrentBlock)) / 2.0;
					blockLength[i] = (float) recombination;
				}
				totalPairLengths[vPath[beginningOfCurrentBlock] - 1] += recombination;
				numberPairLengths[vPath[beginningOfCurrentBlock] - 1] ++;
				beginningOfCurrentBlock = i + 1;
			}
		}

		for (int j = beginningOfCurrentBlock; j < numberOfSnps; j++) {
			blockLength[j] = (float) (recombinationMap.get(numberOfSnps - 1) - recombinationMap.get(beginningOfCurrentBlock));
		}
		totalPairLengths[vPath[beginningOfCurrentBlock] - 1] += recombinationMap.get(numberOfSnps - 1) - recombinationMap.get(beginningOfCurrentBlock);
		numberPairLengths[vPath[beginningOfCurrentBlock] - 1] ++;
		return new LengthsFromPainting(haplotypeIndex, blockLength, vPath, totalPairLengths, numberPairLengths);
	}
	
	
	/**
	 * Takes the values to determine the local emmission probabilities
	 * and returns a double array where of length equal to the number
	 * of haplotypes.  Emmission probabilities calculated using 
	 * Li and Stephens [2003].
	 * @param array {@code byte[][]}
	 * @param haplotypeIndex {@code int}
	 * @param columnIndex {@code int}
	 * @param numberOfHaplotypes {@code int}
	 * @param mutationParameter {@code double}
	 * @param diploid {@code boolean}
	 * @return {@code new double[2]}
	 */
	public static double[] getEmmisionValues (Byte haplotypeRow, Byte[] haplotypeColumn,
			double mutationParameter, ArrayList<Integer> hapsIn) {
		double[] emmisionVector = new double[2];
		emmisionVector[0] = (hapsIn.size() + 0.5 * mutationParameter) / 
				(hapsIn.size() + mutationParameter);
		emmisionVector[1] = (0.5 * mutationParameter) /
				(hapsIn.size() - 1.0 + mutationParameter);
		double[] logProbability = new double[haplotypeColumn.length];
		for (int k : hapsIn) {
			if (haplotypeRow.equals(haplotypeColumn[k])) {
				logProbability[k] = Math.log1p(emmisionVector[0] - 1);
			} else {
				logProbability[k] = Math.log1p(emmisionVector[1] - 1);
			}
		}
		return logProbability;
	}

	/**
	 * Finds the position and values of the max.
	 * element in an array.
	 * @param array {@code double[]}
	 * @return {@link Max}
	 */
	public static Max oneShotMax (double[] array, ArrayList<Integer> haplotypeIndexes) {
		double max = Double.NEGATIVE_INFINITY;
		ArrayList<Integer> maxIndex = new ArrayList<Integer>();
		for (int i : haplotypeIndexes) {
			if (max < array[i]) {
				max = array[i];
			}
		}
		for (int i : haplotypeIndexes) {
			if (array[i] == max) {
				maxIndex.add(i + 1);
			}
		}
		int sIndex = maxIndex.get((int) (Math.random() * maxIndex.size()));
		return new Max(sIndex, max);
	}

	/**
	 * Finds the individual, and the joint haplotypes 
	 * for any haplotype index given a diploid dataset.
	 * @param haplotypeIndex {@code int}
	 * @return {@link IndHap} 
	 */
	public static IndHap haplotypeToIndividual (int haplotypeIndex) {
		int individual = (int) Math.floor(haplotypeIndex / 2.0);
		int haplotypeOne = 2 * individual;
		int haplotypeTwo = 2 * individual + 1;
		return new IndHap (individual, haplotypeOne, haplotypeTwo);
	}

}


