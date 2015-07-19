package bayesfactor;

import java.util.Set;

import utilities.IntegerFlag;
import utilities.StringFlag;

import com.google.common.collect.ImmutableSet;

public class ABACUSFlags {
	
	private static final String HAPLOTYPE_FLAG_NAME =  "haplotype_map";
	private StringFlag haplotypeFlag = new StringFlag(HAPLOTYPE_FLAG_NAME,
			"./genotypes_test.haplotype",
			"Location of haplotype file.");
	
	private static final String GENETIC_MAP_FLAG_NAME = "genetic_map";
	private StringFlag geneticMapFlag = new StringFlag(GENETIC_MAP_FLAG_NAME,
			"./genotypes_test.map",
			"Location of recombination file.");
	
	private static final String COUNTS_FLAG_NAME = "counts";
	private StringFlag countsFlag = new StringFlag(COUNTS_FLAG_NAME,
			"./genotypes_test.viterbiCounts",
			"Location of counts file for population structure adjustment.");
	
	private static final String LENGTHS_FLAG_NAME = "lengths";
	private StringFlag lengthsFlag = new StringFlag(LENGTHS_FLAG_NAME,
			"./genotypes_test.viterbiLengths",
			"Location of lengths file for population structure adjustment.");
	
	private static final String NUMBER_OF_THREADS_FLAG_NAME = "number_of_threads";
	private IntegerFlag numberOfThreadsFlag = new IntegerFlag(NUMBER_OF_THREADS_FLAG_NAME,
			2,
			"Number of threads to run ABACUS.");
	
	private static final Set<String> FLAG_NAMES = ImmutableSet.of(HAPLOTYPE_FLAG_NAME,
			GENETIC_MAP_FLAG_NAME, COUNTS_FLAG_NAME, LENGTHS_FLAG_NAME,
			NUMBER_OF_THREADS_FLAG_NAME);

	public static final double SIGMA_BETA = 1.0;
	public static final int K = 10;
	public static final int NUMBER_OF_LOOPS = 100;
	
	public static ABACUSFlags makeFlags(String[] strings) throws Exception {
		ABACUSFlags flags = new ABACUSFlags();
		for (String s : strings) {
			// --<flag_name>:value
			String[] flag = s.split("--")[1].split(":");
			if (FLAG_NAMES.contains(flag[0])) {
				switch(flag[0]) {
				case HAPLOTYPE_FLAG_NAME : flags.setHaplotype(flag[1]);
				continue;
				case GENETIC_MAP_FLAG_NAME : flags.setGeneticMap(flag[1]);
				continue;
				case COUNTS_FLAG_NAME : flags.setCounts(flag[1]);
				continue;
				case LENGTHS_FLAG_NAME : flags.setLengths(flag[1]);
				continue;
				case NUMBER_OF_THREADS_FLAG_NAME : flags.setTheads(Integer.parseInt(flag[1]));
				continue;
				default : throw new Exception("Flag " + flag[0] + " is not a valid name.\n");
				}
			}
		}
		flags.printDialogue();
		return flags;
	}
	
	private void setHaplotype(String haplotype) {
		haplotypeFlag.setValue(haplotype);
	}

	private void setGeneticMap(String map) {
		geneticMapFlag.setValue(map);
	}
	
	private void setCounts(String counts) {
		countsFlag.setValue(counts);
	}
	
	private void setLengths(String lengths) {
		lengthsFlag.setValue(lengths);
	}
	
	private void setTheads(int threads) {
		numberOfThreadsFlag.setValue(threads);
	}
	
	public String getHaplotype() {
		return haplotypeFlag.getValue();
	}
	public String getGeneticMap() {
		return geneticMapFlag.getValue();
	}
	public String getCountsPath() {
		return countsFlag.getValue();
	}
	public String getLengthsPath() {
		return lengthsFlag.getValue();
	}
	public Integer getNumberOfThreads() {
		return numberOfThreadsFlag.getValue();
	}
	
	void printDialogue() {
		System.out.println("Commencing ABACUS on:");
		System.out.println("\t haplotype: " + haplotypeFlag.getValue());
		System.out.println("***USING " + numberOfThreadsFlag.getValue() + " CORES***");
		System.out.println("Starting data import");
	}
}
