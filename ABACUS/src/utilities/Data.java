package utilities;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;



/**
 * Process and hold raw input data (haplotype file, rsids, 
 * snp positions, recombination map, number of snps, 
 * number of haplotypes)
 * @author Chris Gamble, DPhil Candidate in
 *  Statistical Genetics, University Of Oxford,
 *  Copyright 2014.
 *
 */
public class Data {
	private final boolean[][] haplotype;
	private final List<Integer> snpMap;
	private final List<Double> recombinationMap;
	private final List<String> rsid;
	private final int numberOfSnps;
	private final int numberOfHaplotypes;

	public Data (boolean[][] haplotype, List<Integer> snpMap, List<Double> recombinationMap,
			List<String> rsid, int numberOfSnps, int numberOfHaplotypes) {
		this.haplotype = haplotype;
		this.snpMap = snpMap;
		this.recombinationMap = recombinationMap;
		this.rsid = rsid;
		this.numberOfSnps = numberOfSnps;
		this.numberOfHaplotypes = numberOfHaplotypes;
	}

	public boolean[][] getHaplotype () {
		return haplotype;
	}

	public List<Integer> getSnps () {
		return snpMap;
	}

	public List<Double> getMap () {
		return recombinationMap;
	}

	public List<String> getRsid () {
		return rsid;
	}

	public int numberOfSnps () {
		return numberOfSnps;
	}

	public int numberOfHaplotypes () {
		return numberOfHaplotypes;
	}

	public static Data read(String haplotypePath, String recombinationPath) {
		Size size = determineSize(haplotypePath);
		List<Double> recombinationMap = new ArrayList<>(size.getNumberOfSnps());
		List<Integer> snpMap = new ArrayList<>(size.getNumberOfSnps());
		List<String> rsid = new ArrayList<>(size.getNumberOfSnps());
		boolean[][] haplotype = new boolean[size.getNumberOfSnps()][size.getNumberOfHaplotypes()];
		String rawRow = "";
		try {
			BufferedReader recombinationReader = new BufferedReader(new FileReader(recombinationPath));
			for (int snp = 0; snp < size.getNumberOfSnps(); snp++) {
				rawRow = recombinationReader.readLine();
				if (rawRow == null) {
					throw new IOException("Miss-match between number of SNPs in haplotype file and genetic map.");
				}
				recombinationMap.add(snp, Double.parseDouble(rawRow));
			}
			recombinationReader.close();
		} catch (IOException e) {
			System.out.println("Cannot locate recombination map file. Check input flags -map:./path/to/recombination/file.map");
		}
		try {
			BufferedReader haplotypeReader = new BufferedReader(new FileReader(haplotypePath));
			for (int snp = 0; snp < size.getNumberOfSnps(); snp++) {
				Row row = processRawRowInput(haplotypeReader.readLine());
				haplotype[snp] = row.getHaplotypes();
				snpMap.add(snp, row.getSnpPosition());
				rsid.add(snp, row.getRsid());
			}
			haplotypeReader.close();
		} catch (IOException e) {
			System.out.println("Cannot locate haplotype file. Check input flags -hap:./path/to/haplotype/file.hap");
		}

		return new Data(haplotype, snpMap, recombinationMap, rsid, size.getNumberOfSnps(), size.getNumberOfHaplotypes());
	}

	private static Size determineSize(String haplotypePath) {
		int numberOfHaplotypes = 0;
		int numberOfSnps = 0;
		try {
			BufferedReader haplotypeReader = new BufferedReader(new FileReader(haplotypePath));
			String s = haplotypeReader.readLine();
			Row row = processRawRowInput(s);
			numberOfHaplotypes = row.getHaplotypes().length; 
			numberOfSnps = 1;
			while ((s = haplotypeReader.readLine()) != null) {
				numberOfSnps++;
			}
			haplotypeReader.close();
		} catch (IOException e) {
			System.out.println("Cannot locate haplotype file. Check input flags -hap:./path/to/haplotype/file.hap");
		}
		return new Size(numberOfHaplotypes, numberOfSnps);
	}
	
	private static Row processRawRowInput(String rawRow) {
		String[] row = rawRow.split(" ");
		int numberOfHaplotypes = row.length - 5;
		int snpPosition = Integer.parseInt(row[2]);
		String rsid = row[0];
		boolean[] haplotypes = new boolean[numberOfHaplotypes];
		Arrays.fill(haplotypes, false);
		for (int hap = 0; hap < numberOfHaplotypes; hap++) {
			if (row[hap + 5].equals("1")) {
				haplotypes[hap] = true;
			}
		}
		return new Row(snpPosition, rsid, haplotypes);
	}

	
	private static class Row {
			private final int snpPosition;
			private final String rsid;
			private final boolean[] haplotypes;
			
			public Row(int snpPosition, String rsid, boolean[] haplotypes) {
				this.snpPosition = snpPosition;
				this.rsid = rsid;
				this.haplotypes = haplotypes;
			}
			
			public int getSnpPosition() {
				return snpPosition;
			}
			public String getRsid() {
				return rsid;
			}
			public boolean[] getHaplotypes() {
				return haplotypes;
			}
	}
	
	private static class Size {
		private final int numberOfHaplotypes;
		private final int numberOfSnps;
		
		public Size(int numberOfHaplotypes, int numberOfSnps) {
			this.numberOfHaplotypes = numberOfHaplotypes;
			this.numberOfSnps = numberOfSnps;
		}
		
		public int getNumberOfHaplotypes() {
			return numberOfHaplotypes;
		}
		
		public int getNumberOfSnps() {
			return numberOfSnps;
		}
	}
}
