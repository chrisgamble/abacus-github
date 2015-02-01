package painting;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.google.common.collect.ImmutableSet;

import utilities.Data;

/**
 * Creates the Viterbi Path for each haplotype in the sample.
 * It writes the output files to the same directory as .hap files. 
 * @param -hap:<./path/to/hap>
 * @param -map:<./path/to/map>
 * @param -e:10000 effective population size
 * @param -cores:12 number of cores
 * @param -diploid default haploid
 * @param -ancestry default to output maps, add if require just summary matrices
 * 
 * @author Chris Gamble, DPhil Candidate in
 *  Statistical Genetics, University Of Oxford,
 *  Copyright 2012.
 *
 */
public class PaintingMain {
	
	public static void main (String[] args) {
		try {
			// Set flags.
			Flags flags = new Flags();
			for (String s : args) {
				flags.setFlags(s);
			}
			flags.printDialogue();
			String prefix = flags.getHaplotype().split(".hap")[0];
			BufferedWriter lengths = new BufferedWriter(new 
					FileWriter(prefix + ".viterbiLengths"));
			BufferedWriter counts = new BufferedWriter(new 
					FileWriter(prefix + ".viterbiCounts"));
			
			Data data = Data.read(flags.getHaplotype(), flags.getGeneticMap());
			int numberOfHaplotypes = data.numberOfHaplotypes();
			int numberOfSnps = data.numberOfSnps();
			System.out.println("Number of SNPS " + numberOfSnps);
			System.out.println("Number of Haplotypes " + numberOfHaplotypes);
			double mutationParameter = wattersonsEstimate(numberOfHaplotypes);
			System.out.println("Finished data import, starting painting!");

			ExecutorService executor = Executors.newFixedThreadPool(flags.getNumberOfThreads());
			List<Future<PaintingSummaries>> list = new ArrayList<Future<PaintingSummaries>>();
			for (int haplotypeIndex = 0; haplotypeIndex < numberOfHaplotypes; haplotypeIndex++) {
				PaintingAlgorithm worker = new PaintingAlgorithm(haplotypeIndex, data.getMap(), 
						data.getSnps(), data.getRsid(), data.getHaplotype(), 
						flags.getEffectivePopulation(), mutationParameter, numberOfHaplotypes,
						numberOfSnps, flags.getDiploid(), prefix);
				Future<PaintingSummaries> submit = executor.submit(worker);
				list.add(submit);
			}
			int outputLines = 0;
			for (Future<PaintingSummaries> future : list) {
				try {
					outputLines++;
					System.out.print(outputLines + " ");
					PaintingSummaries threadOutput = future.get();
					lengths.write(arrayToString(threadOutput.getTotal()));
					lengths.newLine();
					counts.write(arrayToString(threadOutput.getNumber()));
					counts.newLine();
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				} catch (NullPointerException e) {
					e.printStackTrace();
				}
			}
			System.out.println();
			System.out.println("Painting complete!");
			executor.shutdownNow();
			lengths.close();
			counts.close();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static class Flags {
		private static final String HAPLOTYPE_FLAG = "haplotype_map";
		private static final String GENETIC_MAP_FLAG = "genetic_map";
		private static final String EFFECTIVE_POPULATION_FLAG = "effective_population";
		private static final String NUMBER_OF_THREADS_FLAG = "number_of_threads";
		private static final String DIPLOID_FLAG = "is_diploid";
		
		private static final Set<String> FLAG_NAMES = ImmutableSet.of(HAPLOTYPE_FLAG,
				GENETIC_MAP_FLAG, EFFECTIVE_POPULATION_FLAG, NUMBER_OF_THREADS_FLAG, DIPLOID_FLAG);
		
		private String haplotype = "./genotypes_test.haplotype";
		private String geneticMap = "./genotypes_test.map";
		private int effectivePopulation = 10000;
		private int numberOfThreads = 2;
		private boolean diploid = false;

		public void setFlags(String s) throws Exception {
			// -<flag_name>=value
			String[] flag = s.split("-")[1].split("=");
			if (FLAG_NAMES.contains(flag[0])) {
				switch(flag[0]) {
				case HAPLOTYPE_FLAG : haplotype = flag[1];
				return;
				case GENETIC_MAP_FLAG : geneticMap = flag[1];
				return;
				case EFFECTIVE_POPULATION_FLAG : effectivePopulation = Integer.parseInt(flag[1]);
				return;
				case NUMBER_OF_THREADS_FLAG : numberOfThreads = Integer.parseInt(flag[1]);
				return;
				case DIPLOID_FLAG : diploid = Boolean.parseBoolean(flag[1]);
				return;
				default : throw new Exception("Flag " + flag[0] + " is not a valid name.\n");
				}
			}
		}

		public String getHaplotype() {
			return haplotype;
		}
		public String getGeneticMap() {
			return geneticMap;
		}
		public int getEffectivePopulation() {
			return effectivePopulation;
		}
		public int getNumberOfThreads() {
			return numberOfThreads;
		}
		public boolean getDiploid() {
			return diploid;
		}
		
		public void printDialogue() {
			System.out.println("Commencing most likely painting algorithm for:");
			System.out.println("\t haplotype: " + haplotype);
			System.out.println("\t recombination map:" + geneticMap);
			if (diploid) {
				System.out.println("\t Treating haplotypes as diploid");
			} else {
				System.out.println("\t Treating haplotypes as haploid");
			}
			System.out.println("\t Effective population size: " + effectivePopulation);
			System.out.println("***USING " + numberOfThreads + " CORES***");
			System.out.println("Starting data import");
		}
	}
	
	private static String arrayToString (float[] values) {
		String str = "";
		for (float d : values) {
			str += d + " ";
		}
		return str;
	}

	private static double wattersonsEstimate (int numberOfHaplotypes) {
		double t = 0;
		for (int i = 1; i < numberOfHaplotypes; i++) {
			t += 1.0 / ((double) i);
		}
		return 1.0 / t;
	}
}