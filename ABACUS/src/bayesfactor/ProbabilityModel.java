package bayesfactor;
import java.util.concurrent.Callable;

import jsc.distributions.Gamma;

/**
 * Functional class to calculate the posterior probabilities
 * of being long, given the pairwise genome wide average.  Using these
 * probabilities we perform two tests of association; the ratio of Monte Carlo
 * approximated marginal likelihoods and the average of the Monte Carlo Bayes Factors.
 * 
 * @author Chris Gamble, DPhil Candidate in
 *  Statistical Genetics, University Of Oxford,
 *  Copyright 2012.
 *
 */
public class ProbabilityModel implements Callable<BayesfactorOut> {
	private final double[] lengths;
	private final int[] viterbi;
	private final boolean[] alleles;
	private final float[][] average;
	private final int snpIndex;
	private final double sigmaBeta;
	private final int K;

	public ProbabilityModel(double[] lengths, int[] viterbi, boolean[] alleles, 
			float[][] average, int snpIndex, double sigmaBeta, int K) {
		this.lengths = lengths;
		this.viterbi = viterbi;
		this.alleles = alleles;
		this.average = average;
		this.snpIndex = snpIndex;
		this.sigmaBeta = sigmaBeta;
		this.K = K;
	}

	public BayesfactorOut call() {
		LaplaceApproximation laplaceApproximation = new LaplaceApproximation(snpIndex,
				getProb(getMean(average)), alleles, sigmaBeta, K);
		return laplaceApproximation.getBF();
	}
	
	/**
	 * Takes Genome wide average lengths, and finds the relevant means
	 * for the current SNP  
	 * @param {@link float[][] average}
	 * @return {@link double[] mean}
	 */
	private double[] getMean (float[][] average) { 
		double[] mean = new double[lengths.length];
		for (int haplotype = 0; haplotype < lengths.length; haplotype++) {
			mean[haplotype] = average[haplotype][viterbi[haplotype] - 1];
		}
		return mean;
	}
	
	/**
	 * Takes the local mean as a parameter, this function calculates the 
	 * posterior probability of begin long, and the two tests of
	 * association 
	 * @param {@link double[]} local mean
	 * @return {@link BayesFactorPair} Object containing the two tests and snp location
	 */
	private double[] getProb (double[] mean) {
		double eTheta = 0.5;
		double[] probLong = new double[lengths.length];
		for (int i = 0; i < lengths.length; i++) {
			if (lengths[i] > 0 & mean[i] > 0) {
				Gamma localDist = new Gamma(2, lengths[i] / mean[i]);
				double logA = Math.log(mean[i]) - 2.0 * Math.log(lengths[i]) + Math.log(localDist.cdf(1.0)) + Math.log(eTheta);
				double logB = -1.0 * Math.log(mean[i]) - lengths[i] / mean[i] + Math.log(1 - eTheta);
				probLong[i] = Math.exp(-1.0 * Math.log(1 + Math.exp(logB - logA)));
			} else {
				probLong[i] = 0.0;		
			}
		}
		return probLong;
		
	}
}
