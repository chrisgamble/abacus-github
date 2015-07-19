package bayesfactor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import utilities.Data;

public class ABACUSMain {
	/**
	 * @author Chris Gamble, DPhil Student, University of Oxford
	 * @param -hap:pathtohapfile
	 * @param -lengths:pathtolengthfile
	 * @param -counts:pathtocountfile
	 * @param -cores:numberOfThreads
	 * @param -debug
	 */
	public static void main (String[] args) {
		try {
			ABACUSFlags flags = ABACUSFlags.makeFlags(args);
			flags.printDialogue();
			String prefix = flags.getHaplotype().split(".hap")[0];
			Data data = Data.read(flags.getHaplotype(), flags.getGeneticMap());
			System.out.println("Finished data import, starting ABACUS!");
			ArrayList<Integer> start = new ArrayList<Integer>();
			ArrayList<Integer> end = new ArrayList<Integer>();
			int numberOfHaplotypes = data.numberOfHaplotypes();
			int numberOfSnps = data.numberOfSnps();
			int size = numberOfSnps / ABACUSFlags.NUMBER_OF_LOOPS; 
			for (int loop = 0; loop < ABACUSFlags.NUMBER_OF_LOOPS - 1; loop++) {
				start.add(loop, loop * size);
				end.add(loop, (loop + 1) * size);
			}
			start.add(ABACUSFlags.NUMBER_OF_LOOPS - 1, (ABACUSFlags.NUMBER_OF_LOOPS - 1) * size);
			end.add(ABACUSFlags.NUMBER_OF_LOOPS - 1, numberOfSnps);
			
			float[][] counts = readArray(flags.getCountsPath(), numberOfHaplotypes);
			float[][] lengths = readArray(flags.getLengthsPath(), numberOfHaplotypes);
			
			LaplacePair arrayAndMean = getPairArray(lengths, counts);
			
			float[][] average = arrayAndMean.getArray();
			
			Map<Integer, BufferedReader> readers = new HashMap<Integer, BufferedReader>();
			for (int haplotypeIndex = 0; haplotypeIndex < numberOfHaplotypes; haplotypeIndex++) {
				String path = prefix + "_" + haplotypeIndex + ".viterbi";
				readers.put(haplotypeIndex, new BufferedReader(new FileReader(path)));
			}
			BufferedWriter stochasticPair = new BufferedWriter(new FileWriter(prefix + ".bayesfactor"));
			for (int loop = 0; loop < ABACUSFlags.NUMBER_OF_LOOPS; loop++) {
				ExecutorService executor = Executors.newFixedThreadPool(flags.getNumberOfThreads());
				List<Future<BayesfactorOut>> list = new ArrayList<Future<BayesfactorOut>>();
				for (int snp = start.get(loop); snp < end.get(loop); snp++) {
					double[] snpLength = new double[numberOfHaplotypes];
					int[] viterbi = new int[numberOfHaplotypes];
					for (int haplotypeIndex = 0; haplotypeIndex < numberOfHaplotypes; haplotypeIndex++) {
						String str = readers.get(haplotypeIndex).readLine();
						String[] tmp = str.split(" ");
						snpLength[haplotypeIndex] = Double.parseDouble(tmp[3]);
						viterbi[haplotypeIndex] = Integer.parseInt(tmp[2]);
					}
					ProbabilityModel worker = new ProbabilityModel(snpLength, viterbi,
							data.getHaplotype()[snp], average, snp, ABACUSFlags.SIGMA_BETA,
							ABACUSFlags.K); 
					Future<BayesfactorOut> submit = executor.submit(worker);
					list.add(submit);
				}
				for (Future<BayesfactorOut> future : list) {
					try {
						BayesfactorOut threadOutput = future.get();
						int snp = threadOutput.getSnp();
						stochasticPair.write(data.getRsid().get(snp) +
								" " + data.getSnps().get(snp) +
								" " + threadOutput.getBFAdj() +
								" " + threadOutput.getBetaAdj());
						stochasticPair.newLine();
					} catch (InterruptedException e) {
						e.printStackTrace();
					} catch (ExecutionException e) {
						e.printStackTrace();
					}
				}
				System.out.println(((int) (loop * 100.0 / ABACUSFlags.NUMBER_OF_LOOPS) + 1) + "%");
				executor.shutdownNow();
			}
			System.out.println("ABACUS Complete!");
			stochasticPair.close();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static float[][] readArray(String path, int numberOfHaplotypes) throws IOException {
		BufferedReader textReader = new BufferedReader(new FileReader(path));
		int haplotype = 0;
		String strLine;
		float[][] output = new float[numberOfHaplotypes][numberOfHaplotypes];
		while ((strLine = textReader.readLine()) != null) {
			String[] local = strLine.split(" ");
			for (int i = 0; i < local.length; i++) {
				output[haplotype][i] = Float.parseFloat(local[i]);
			}
			haplotype++;
		}
		textReader.close();
		return output;
	}
		
	public static LaplacePair getPairArray(float[][] lengths, float[][] counts) {
		float[][] average = new float[lengths.length][lengths.length];
		double[] genomeWideMean = new double[lengths.length];
		float mean = 0;
		for (int i = 0; i < lengths.length; i++) {
			for (int j = 0; j < lengths.length; j++) {
				if (!Float.isNaN(lengths[i][j]) && !Float.isNaN(counts[i][j]) && counts[i][j] > 0) {
					average[i][j] = lengths[i][j] / counts[i][j];
					mean += average[i][j];
				}
			}
		}
		for (int i = 0; i < lengths.length; i++) {
			genomeWideMean[i] = mean;
		}
		return new LaplacePair(average, genomeWideMean);
	}
}