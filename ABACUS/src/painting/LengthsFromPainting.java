package painting;


/**
 * Generic Class to hold max value and position in an array
 * @author Chris Gamble, DPhil Candidate in
 *  Statistical Genetics, University Of Oxford,
 *  Copyright 2012.
 *
 */
public class LengthsFromPainting {
	private int index;
	private float[] blocks;
	private short[] path;
	private float[] lengths;
	private float[] counts;
 
	public LengthsFromPainting (int index, float[] blocks, short[] path, float[] lengths, float[] counts) {
		this.index = index;
		this.blocks = blocks;
		this.path = path;
		this.lengths = lengths;
		this.counts = counts;
	}

	public int getIndex () {
		return index;
	}
	public float[] getBlocks () {
		return blocks;
	}
	
	public short[] getPath () {
		return path;
	}

	public float[] getLengths () {
		return lengths;
	}
	
	public float[] getCounts () {
		return counts;
	}
}

