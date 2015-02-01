package utilities;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;

/**
 * Functional Class to Process Recombination Map for each chromosome.
 * Holds the processed recombination map and snp map. 
 * @author Chris Gamble, DPhil Candidate in
 *  Statistical Genetics, University Of Oxford,
 *  Copyright 2014.
 *
 */
public class GetRecombinationMap {
	private ArrayList<Double> processedMap;
	private ArrayList<Integer> snpMap;

	public GetRecombinationMap (ArrayList<Double> processedMap, ArrayList<Integer> snpMap) {
		this.processedMap = processedMap;
		this.snpMap = snpMap;
	}

	public ArrayList<Double> getRecombinationMap() {
		return processedMap;
	}

	public ArrayList<Integer> getSnpMap () {
		return snpMap;
	}

	private static class MappedRecombination {
		private int snp;
		private double interSnpRecombination;
		private double cumulativeRecombination;

		public MappedRecombination (int snp, double interSnpRecombination, double cumulativeRecombination) {
			this.snp = snp;
			this.interSnpRecombination = interSnpRecombination;
			this.cumulativeRecombination = cumulativeRecombination;
		}

		public int getSnp () {
			return snp;
		}

		public double getInterSnpRecombination () {
			return interSnpRecombination;
		}

		public double getCumulativeRecombination () {
			return cumulativeRecombination;
		}
	}

	public static class Interpolate {
		private int lastInt;
		private double value;

		public Interpolate(int lastInt, double value) {
			this.lastInt = lastInt;
			this.value = value;
		}

		public int getLastInt () {
			return lastInt;
		}

		public double getValue () {
			return value;
		}
	}

	public static void main (String args[]) {
		String hap = "/Users/chrisgamble/Downloads/chr19.sub.dom.haps";
		String map = "/Users/chrisgamble/Downloads/dom_chr19_cleaned_cM.hapmapFormat.txt";
		for (String s : args) {
			if (s.contains("-hap")) {
				hap = s.split(":")[1];
			}
			if (s.contains("-map")) {
				map = s.split(":")[1];
			}
		}
		try {
			String prefix = hap.split(".hap")[0];
			ArrayList<Integer> snps = readSnpMap(prefix + ".snp");
			ArrayList<String> recombinationMap = read(map);
			GetRecombinationMap processedRecombination = fixRecombinationMap(create(recombinationMap, snps));
			ArrayList<Double> processedMap = processedRecombination.getRecombinationMap();
			BufferedWriter recombPath = new BufferedWriter(new FileWriter(prefix + ".map"));
			for (int i = 0; i < processedMap.size(); i++) {
				recombPath.write(Double.toString(processedMap.get(i)));
				recombPath.newLine();
			}
			recombPath.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static ArrayList<MappedRecombination> make (ArrayList<Integer> snps, Map<Integer, Double> map) {
		ArrayList<MappedRecombination> out = new ArrayList<MappedRecombination>();
		ArrayList<Integer> mappedSnps = new ArrayList<Integer>(map.keySet());
		int counter = 0;
		for (int snp : snps) {
			System.out.println(snp);
			int s = 0;
			while (mappedSnps.get(s) < snp && s < mappedSnps.size() - 1) {
				s++;
			}
			out.add(counter, new MappedRecombination(snp, map.get(mappedSnps.get(s)), snp * map.get(mappedSnps.get(s)) / 1E6));
			counter++;
		}
		return out;
	}

	/**
	 * Fucntion creates the recombination map, and checks there
	 *  are no leading zeros.  If there are it interpolates these
	 *   values according to the first recombination distance and
	 *    the relative distance of each SNP to the position of the
	 *     first recombination distance. 
	 * @param recombinationMap {@code ArrayList<String>}
	 * @param haplotypeLegend {@code ArrayList<String>}
	 * @return {@link GetRecombinationMap}
	 */
	public static GetRecombinationMap create (ArrayList<String> recombinationMap,
			ArrayList<Integer> haplotypeLegend) {

		//		int numberOfSnps = haplotypeLegend.size() - 1;
		int numberOfSnps = haplotypeLegend.size();
		int numberOfMappedSnps = recombinationMap.size();
		ArrayList<Integer> legend = new ArrayList<Integer>();
		ArrayList<MappedRecombination> mappedRecombination = new ArrayList<MappedRecombination>();
		ArrayList<Double> processedMap = new ArrayList<Double>();
		for (int i = 1; i < numberOfMappedSnps - 1; i++) {
			String[] mappedSnpIndexValue = recombinationMap.get(i).split(" ");
			mappedRecombination.add(i - 1, new MappedRecombination(Integer.parseInt(mappedSnpIndexValue[0]),
					Double.parseDouble(mappedSnpIndexValue[1]), Double.parseDouble(mappedSnpIndexValue[2])));
		}
		int lastInt = 0;
		for (int i = 0; i < numberOfSnps; i++) {
			int haplotypeKey = haplotypeLegend.get(i);
			Interpolate values = getRecombinationValuesFixed(mappedRecombination, haplotypeKey, lastInt);
			lastInt = values.getLastInt();
			processedMap.add(values.getValue());
			legend.add(haplotypeKey);
		}
		System.out.println("Number of snps in legend : " + legend.size());
		System.out.println("Number of snps in map : " + processedMap.size());
		return new GetRecombinationMap(processedMap, legend);

	}

	/**
	 * Generic holding class to help speed up interpolation.
	 * @author Chris Gamble, DPhil Candidate in
	 *  Statistical Genetics, University Of Oxford,
	 *  Copyright 2012.
	 *
	 */


	/**
	 * Method to interpolate leading zeros.
	 * @param originalValues {@code ProbabilityData}
	 * @return return {@code ProbabilityData}
	 */
	public static GetRecombinationMap fixRecombinationMap (GetRecombinationMap originalValues) {
		int i;
		for (i = 0; i < originalValues.getRecombinationMap().size(); i++) {
			if (originalValues.getRecombinationMap().get(i) > 0.0) {
				break;
			}
		}
		if (i > 0) {
			int firstNonZeroSnp = originalValues.getSnpMap().get(i);
			double firstNonZeroValue = originalValues.getRecombinationMap().get(i);
			for (int j = 0; j < i; j++) {
				double interpolatedValue = originalValues.getSnpMap().get(j) / (double) firstNonZeroSnp * firstNonZeroValue;
				originalValues.getRecombinationMap().set(j, interpolatedValue);
			}
		}
		return originalValues;
	}

	/**
	 * Method to interpolate between mapped SNPs
	 * @param mappedSnps {@code ArrayList<Integer>}
	 * @param mappedRecombination {@code ArrayList<Double>}
	 * @param currentSnpOfInterest {@code int}
	 * @param lastIndex {@code int}
	 * @return {@link Interpolate}
	 */
	public static Interpolate getRecombinationValuesFixed (ArrayList<MappedRecombination> mappedRecombination,
			int currentSnpOfInterest, int lastIndex) {
		double value = 0.0;
		int lastInt = 0;
		if (currentSnpOfInterest <= mappedRecombination.get(0).getSnp()) {
			lastInt = 0;
			value = 0;
			return new Interpolate(lastInt, value);
		} else {
			if (currentSnpOfInterest >= mappedRecombination.get(mappedRecombination.size() - 1).getSnp()) {
				lastInt = mappedRecombination.size() - 1;
				value = (currentSnpOfInterest - mappedRecombination.get(lastInt).getSnp()) / 1000000.0 *
						mappedRecombination.get(lastInt).getInterSnpRecombination() + mappedRecombination.get(lastInt).getCumulativeRecombination();
				return new Interpolate(lastInt, value);
			} else {
				int start = lastIndex - 10;
				if (start <= 0) {
					start = 0;
				}
				for (int i = start; i < mappedRecombination.size() - 1; i++) { 
					if (mappedRecombination.get(i).getSnp() <= currentSnpOfInterest &&
							mappedRecombination.get(i + 1).getSnp() > currentSnpOfInterest){
						lastInt = i;
						int recStart = i - 1;
						if (recStart <= 0) {
							recStart = 0;
						}
						value = (currentSnpOfInterest - mappedRecombination.get(i).getSnp()) / 1000000.0 *
								mappedRecombination.get(i + 1).getInterSnpRecombination() + mappedRecombination.get(i).getCumulativeRecombination();
						return new Interpolate(lastInt, value);
					}
				}
			}
		}
		return new Interpolate(lastInt, value);
	}

	public static ArrayList<Double> readProcessedMap (String path) throws IOException {
		BufferedReader recombPath = new BufferedReader(new FileReader(path));
		ArrayList<Double> map = new ArrayList<Double>();
		String strLine;
		while ((strLine = recombPath.readLine()) != null) {
			map.add(Double.parseDouble(strLine));
		}
		recombPath.close();
		return map;
	}

	public static ArrayList<Integer> readSnpMap (String path) throws IOException {
		BufferedReader recombPath = new BufferedReader(new FileReader(path));
		ArrayList<Integer> map = new ArrayList<Integer>();
		String strLine;
		while ((strLine = recombPath.readLine()) != null) {
			map.add(Integer.parseInt(strLine));
		}
		recombPath.close();
		return map;
	}

	public static ArrayList<String> read(String path) throws IOException {

		BufferedReader textReader = new BufferedReader(new FileReader(path));

		String strLine;
		ArrayList<String> output = new ArrayList<String>();
		while ((strLine = textReader.readLine()) != null) {
			output.add(strLine);
		}
		textReader.close();
		return output;
	}
}