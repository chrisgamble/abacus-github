package utilities;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;


public class CombineSummaries {
	public static void main (String[] args) {
		String first = "";
		String second = "";
		String output = "";
		for (String s : args) {
			if (s.contains("-first")) {
				first = s.split(":")[1];
			}
			if (s.contains("-second")) {
				second = s.split(":")[1];
			}
			if (s.contains("-output")) {
				output = s.split(":")[1];
			}
		
			
		}
		
		try {
			writeArray(output, getTArray(readArray(first), readArray(second)));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static float[][] readArray(String path) throws IOException {
		BufferedReader textReader = new BufferedReader(new FileReader(path));
		ArrayList<String[]> storage = new ArrayList<String[]>();
		int numberOfHaplotypes = 0;
		String strLine = "";
		while ((strLine = textReader.readLine()) != null) {
			storage.add(numberOfHaplotypes, strLine.split(" "));
			numberOfHaplotypes++;
		}
		float[][] output = new float[numberOfHaplotypes][numberOfHaplotypes];
		for (int i = 0; i < numberOfHaplotypes; i++) {
			for (int j = 0; j < numberOfHaplotypes; j++) {
				output[i][j] = Float.parseFloat(storage.get(i)[j]);
			}
		}
		textReader.close();
		return output;
	}
	
	public static void writeArray(String path, float[][] array) throws IOException {
		BufferedWriter textReader = new BufferedWriter(new FileWriter(path));
		for (float[] f : array) {
			String str = "";
			for (float fd : f) {
				str += fd + " "; 
			}
			textReader.write(str);
			textReader.newLine();
		}
		textReader.close();
	}
	
	public static float[][] getTArray (float[][] A, float[][] B) {
		float[][] T = new float[A.length][A.length];
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A.length; j++) {
				if (!Float.isNaN(A[i][j]) && !Float.isNaN(B[i][j])) {
					T[i][j] = A[i][j] + B[i][j];
				}
			}
		}
		return T;
	}
	
}