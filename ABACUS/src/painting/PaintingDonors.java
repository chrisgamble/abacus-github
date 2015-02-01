package painting;
import java.util.ArrayList;


public class PaintingDonors {
	private ArrayList<Integer> haplotypeIndexes;
	private ArrayList<Integer> notHaplotypeIndexes;

	public PaintingDonors (ArrayList<Integer> haplotypeIndexes, ArrayList<Integer> notHaplotypeIndexes) {
		this.haplotypeIndexes = haplotypeIndexes;
		this.notHaplotypeIndexes = notHaplotypeIndexes;
	}

	public ArrayList<Integer> getHaplotypeIndexes () {
		return haplotypeIndexes; // recipients
	}

	public ArrayList<Integer> getNotHaplotypeIndexes () {
		return notHaplotypeIndexes; //donors
	}

	public static PaintingDonors get (boolean diploid, boolean donor,
			ArrayList<Integer> donorIndexes, int numberOfHaplotypes, int haplotypeIndex) {
		ArrayList<Integer> haplotypeIndexes = new ArrayList<Integer>();
		ArrayList<Integer> notHaplotypeIndexes = new ArrayList<Integer>();
		if (donor) {
			for (int k = 0; k < numberOfHaplotypes; k++) {
				if (donorIndexes.contains(k) | k == haplotypeIndex) {
					notHaplotypeIndexes.add(k);
				} else {
					haplotypeIndexes.add(k);
				}
			}
		} else {
			if (diploid) {
				IndHap ind = haplotypeToIndividual (haplotypeIndex);
				notHaplotypeIndexes.add(ind.getHaplotypeOne());
				notHaplotypeIndexes.add(ind.getHaplotypeTwo());
				for (int k = 0; k < numberOfHaplotypes; k++) {
					if (k != ind.getHaplotypeOne() && k != ind.getHaplotypeTwo()) {
						haplotypeIndexes.add(k);
					}
				}
			} else {
				for (int k = 0; k < numberOfHaplotypes; k++) {
					if (k != haplotypeIndex) {
						haplotypeIndexes.add(k);
					} else {
						notHaplotypeIndexes.add(haplotypeIndex);
					}
				}
			}
		}
		return new PaintingDonors(haplotypeIndexes, notHaplotypeIndexes);
	}

	public static PaintingDonors getRecipients (boolean diploid, boolean donor,
			ArrayList<Integer> donorIndexes, int numberOfHaplotypes, int haplotypeIndex) {
		ArrayList<Integer> haplotypeIndexes = new ArrayList<Integer>();
		ArrayList<Integer> notHaplotypeIndexes = new ArrayList<Integer>();
		for (int k = 0; k < numberOfHaplotypes; k++) {
			if (donorIndexes.contains(k)) {
				notHaplotypeIndexes.add(k);
			} else {
				haplotypeIndexes.add(k);
			}
		}
		return new PaintingDonors(haplotypeIndexes, notHaplotypeIndexes);
	}
	
	public static IndHap haplotypeToIndividual (int haplotypeIndex) {
		int individual = (int) Math.floor(haplotypeIndex / 2.0);
		int haplotypeOne = 2 * individual;
		int haplotypeTwo = 2 * individual + 1;
		return new IndHap (individual, haplotypeOne, haplotypeTwo);
	}
}
