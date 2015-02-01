package bayesfactor;

public class LaplaceApproximation {
	private final int snpIndex;
	private final double[] probabilityLong; 
	private final boolean[] alleles;
	private final double sigmaBeta;
	private final int K;
	
	public LaplaceApproximation(int snpIndex, double[] probabilityLong, boolean[] alleles,
			double sigmaBeta, int K) {
		this.snpIndex = snpIndex;
		this.probabilityLong = probabilityLong;
		this.alleles = alleles;
		this.sigmaBeta = sigmaBeta;
		this.K = K;
	}
	
	public BayesfactorOut getBF () {
		double[] thetaH1 = mleH1();
		double thetaH0 = mleH0();
		double marginalH1 = likelihood(thetaH1[0], thetaH1[1]) +
				-0.5 * Math.log (2 * Math.PI) - 0.5 * Math.pow(thetaH1[0], 2) +
				-0.5 * Math.log (2 * Math.PI * Math.pow(sigmaBeta, 2)) - 0.5 * Math.pow(thetaH1[1] / sigmaBeta, 2) +
				-0.5 * Math.log(det(IH1(thetaH1[0], thetaH1[1]))); 
		double marginalH0 = likelihood(thetaH0, 0.0) +
				-0.5 * Math.log (2 * Math.PI) - 0.5 * Math.pow(thetaH0, 2) +
				-0.5 * Math.log(Math.abs(IH0(thetaH0)));
		double logEBF = marginalH1 - marginalH0;
		double log10BF = Math.log10(Math.exp(logEBF));
		return new BayesfactorOut(snpIndex, log10BF, thetaH1[1]);
	}

	
	
	private double likelihood(double mu, double beta) {
		double l = 0.0;
		for (int i = 0; i < probabilityLong.length; i++) {
			l += Math.log(f(alleles[i], mu, beta, probabilityLong[i]));
		}
		return l;
	}
	
	private double[] mleH1 () {
		double[] theta = new double[2];
		theta[0] = 0.0;
		theta[1] = 0.0;
		for (int i = 1; i < K; i++) {
			double[] change = multi(inverse(IH1(theta[0], theta[1])), UH1(theta[0], theta[1]));
			theta[0] += -1.0 * change[0];
			theta[1] += -1.0 * change[1];
		}
		return theta;
	}
	
	private double mleH0 () {
		double theta = 0.0;
		for (int i = 1; i < K; i++) {
			double change = UH0(theta) / IH0(theta);
			theta += -1.0 * change;
		}
		return theta;
	}
	
	private double[] UH1 (double mu, double beta) {
		double[] theta = new double[2];
		theta[0] = -1.0 * mu;
		theta[1] = -1.0 * beta / Math.pow(sigmaBeta, 2);
		for (int i = 0; i < probabilityLong.length; i++) {
			double fI = f(alleles[i], mu, beta, probabilityLong[i]);
			theta[0] += firstDerivFMu (alleles[i], mu, beta, probabilityLong[i]) / fI;
			theta[1] += firstDerivFBeta (alleles[i], mu, beta, probabilityLong[i]) / fI;
		}
		return theta;
	}

	private double UH0 (double mu) {
		double theta = -1.0 * mu;
		for (int i = 0; i < probabilityLong.length; i++) {
			double fI = f(alleles[i], mu, 0.0, probabilityLong[i]);
			theta += firstDerivFMu (alleles[i], mu, 0.0, probabilityLong[i]) / fI;
		}
		return theta;
	}

	private double[][] IH1 (double mu, double beta) {
		double[][] out = new double[2][2];
		out[0][0] = -1.0;
		out[0][1] = 0.0;
		out[1][0] = 0.0;
		out[1][1]  = -1.0 / Math.pow(sigmaBeta, 2);
		for (int i = 0; i < probabilityLong.length; i++) {
			double fI = f(alleles[i], mu, beta, probabilityLong[i]);
			double fDFMu = firstDerivFMu (alleles[i], mu, beta, probabilityLong[i]);
			double fDFBeta = firstDerivFBeta (alleles[i], mu, beta, probabilityLong[i]);
			double sDFMu = secondDerivFMu (alleles[i], mu, beta, probabilityLong[i]);
			double sDFBeta = secondDerivFBeta (alleles[i], mu, beta, probabilityLong[i]);
			double fDFBetaMu = firstDerivFBetaMu (alleles[i], mu, beta, probabilityLong[i]);
			out[0][0] += (sDFMu * fI - Math.pow(fDFMu, 2)) / Math.pow(fI, 2);
			out[1][1] += (sDFBeta * fI - Math.pow(fDFBeta, 2)) / Math.pow(fI, 2);
			out[0][1] += (fDFBetaMu * fI - fDFMu * fDFBeta) / Math.pow(fI, 2);
		}
		out[1][0] = out[0][1];
		return out;
	}

	private double IH0 (double mu) {
		double out = -1.0;
		for (int i = 0; i < probabilityLong.length; i++) {
			double fI = f(alleles[i], mu, 0.0, probabilityLong[i]);
			double fDFMu = firstDerivFMu (alleles[i], mu, 0.0, probabilityLong[i]);
			double sDFMu = secondDerivFMu (alleles[i], mu, 0.0, probabilityLong[i]);
			out += (sDFMu * fI - Math.pow(fDFMu, 2)) / Math.pow(fI, 2);;
		}
		return out;
	}
	
	private double det (double[][] A) {
		return A[0][0] * A[1][1] - A[1][0] * A[0][1];
	}

	private double[][] inverse (double[][] A) {
		double[][] inverse = new double[2][2];
		double d = det(A);
		inverse[0][0] = A[1][1] / d;
		inverse[0][1] = -1.0 * A[0][1] / d;
		inverse[1][0] = -1.0 * A[1][0] / d;
		inverse[1][1] = A[0][0] / d;
		return inverse;
	}
	
	private double[] multi (double[][] A, double[] B) {
		double[] out = new double[2];
		out[0] = A[0][0] * B[0] + A[0][1] * B[1];
		out[1] = A[1][0] * B[0] + A[1][1] * B[1];
		return out;
	}
	
	private double f (boolean allele, double mu, double beta, double probLong) {
		double p1 = Math.exp(mu + beta) / (1.0 + Math.exp(mu + beta));
		double p0 = Math.exp(mu) / (1.0 + Math.exp(mu));
		double out = 0.0;
		if (allele) {
			out = p1 * probLong + p0 * (1.0 - probLong);
		} else {
			out = 1.0 - (p1 * probLong + p0 * (1.0 - probLong));
		}
		return out;
	}
	
	private double p1FDBeta (double mu, double beta) {
		return Math.exp(mu + beta) / Math.pow(1.0 + Math.exp(mu + beta), 2);
	}

	private double p1SDBeta (double mu, double beta) {
		return Math.exp(mu + beta) * (1.0 - Math.exp(mu + beta)) /
				Math.pow(1 + Math.exp(mu + beta), 3);
	}

	private double p0FDBeta  (double mu, double beta) {
		return 0.0;
	}

	private double p1FDMu (double mu, double beta) {
		return Math.exp(mu + beta) / Math.pow(1 + Math.exp(mu + beta), 2);
	}

	private double p1SDMu (double mu, double beta) {
		return Math.exp(mu + beta) * (1 - Math.exp(mu + beta)) /
				Math.pow(1 + Math.exp(mu + beta), 3);
	}

	private double p0FDMu (double mu, double beta) {
		return Math.exp(mu) / Math.pow(1.0 + Math.exp(mu), 2);
	}

	private double p0SDMu (double mu, double beta) {
		return Math.exp(mu) * (1 - Math.exp(mu)) / Math.pow(1 + Math.exp(mu), 3);
	}

	private double p1FDMuFDBeta (double mu, double beta) {
		return Math.exp(mu + beta) * (1 - Math.exp(mu + beta)) / Math.pow(1 + Math.exp(mu + beta), 3);
	}

	private double firstDerivFMu (boolean allele, double mu, double beta, double probLong) {
		double out = 0.0;
		if (allele) {
			out = p1FDMu (mu, beta) * probLong + p0FDMu (mu, beta) * (1.0 - probLong);
		} else {
			out = -1.0 * (p1FDMu (mu, beta) * probLong + p0FDMu (mu, beta) * (1.0 - probLong));
		}
		return out;
	}

	private double firstDerivFBeta (boolean allele, double mu, double beta, double probLong) {
		double out = 0.0;
		if (allele) {
			out = p1FDBeta (mu, beta) * probLong;
		} else {
			out = -1.0 * p1FDBeta (mu, beta) * probLong;
		}
		return out;
	}

	private double secondDerivFMu (boolean allele, double mu, double beta, double probLong) {
		double out = 0.0;
		if (allele) {
			out = p1SDMu (mu, beta) * probLong + p0SDMu (mu, beta) * (1.0 - probLong);
		} else {
			out = -1.0 * (p1SDMu (mu, beta) * probLong + p0SDMu (mu, beta) * (1.0 - probLong));
		}
		return out;
	}

	private double secondDerivFBeta (boolean allele, double mu, double beta, double probLong) {
		double out = 0.0;
		if (allele) {
			out = p1SDBeta (mu, beta) * probLong;
		} else {
			out = -1.0 * p1SDBeta (mu, beta) * probLong;
		}
		return out;
	}

	private double firstDerivFBetaMu (boolean allele, double mu, double beta, double probLong) {
		double out = 0.0;
		if (allele) {
			out = p1FDMuFDBeta (mu, beta) * probLong;
		} else {
			out = -1.0 * p1FDMuFDBeta (mu, beta) * probLong;
		}
		return out;
	}	
}
