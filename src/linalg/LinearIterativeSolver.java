package linalg;

public interface LinearIterativeSolver extends LinearSolverInterface {
	static final int NO_STATISTICS=0; // do not show statistics
	static final int SHOW_ONLY_TURNS=1; // show iterations
	static final int SHOW_ALL=2; // shall all details
	
	public static boolean COUNT=true; // chose to count iterations
	public static final int[] GLOBAL_COUNTS = new int[1024]; // iterations counter
	
	public void setVerbosity(int verbosity);
	
	/** Set the numeric tolerance */
	public void setPrecision(double precision);
	
	/** Get the numeric tolerance */
	public double getPrecision();
	
	/** Set the starting vector (initial guess) */
	public void setStartingVector(double[] x);

	/** 
	 * Run the iterative solver k times (solve the system Ax=b)
	 * 
	 * @param b  right side vector
	 * @param k  number of iterations to perform 
	 * @return	the solution of Ax=b
	 * */
	public double[] runKtimes(double[] b, int k);
	
	/** 
	 * Check whether Ax=b
	 * 
	 *  @param b
	 *  		right side vector
	 *  @param x
	 *  		solution to check
	 *  @param precision
	 *  		numeric tolerance
	 *  
	 *  @return true if dist(Ax, b) < precision
	 **/
	public boolean checkSolution(double[] x, double[] b, double precision);
}
