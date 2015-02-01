package painting;


/**
 * Generic Class to hold Individual to Haplotypes mapping
 * @author Chris Gamble, DPhil Candidate in
 *  Statistical Genetics, University Of Oxford,
 *  Copyright 2012.
 *
 */
public class IndHap {
	private int individual;
	private int haplotypeOne;
	private int haplotypeTwo;
	
	public IndHap (int individual, int haplotypeOne, int haplotypeTwo) {
		this.haplotypeOne = haplotypeOne;
		this.haplotypeTwo = haplotypeTwo;
		this.individual = individual;
	}
	
	public int getHaplotypeOne () {
		return haplotypeOne;
	}
	
	public int getHaplotypeTwo () {
		return haplotypeTwo;
	}
	
	public int getIndividual () {
		return individual;
	}
	
	public static IndHap create (int haplotypeIndex) {
		int individual = (int) Math.floor(haplotypeIndex / 2.0);
		int haplotypeOne = 2 * individual;
		int haplotypeTwo = 2 * individual + 1;
		return new IndHap (individual, haplotypeOne, haplotypeTwo);
	}
}
