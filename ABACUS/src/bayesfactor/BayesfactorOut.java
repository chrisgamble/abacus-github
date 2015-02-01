package bayesfactor;

public class BayesfactorOut {
	private int snp;
	private double bayesFactorIndAdj;
	private double betaHatH1IndAdj;
	
	public BayesfactorOut (int snp, double bayesFactorIndAdj, double betaHatH1IndAdj) {
		this.snp = snp;
		this.bayesFactorIndAdj = bayesFactorIndAdj;
		this.betaHatH1IndAdj = betaHatH1IndAdj; 
	}
	
	public int getSnp () {
		return snp;
	}
	
	public double getBFAdj () {
		return bayesFactorIndAdj;
	}
	
	public double getBetaAdj () {
		return betaHatH1IndAdj;
	}
	
}
