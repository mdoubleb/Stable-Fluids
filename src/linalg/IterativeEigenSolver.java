package linalg;

public interface IterativeEigenSolver extends EigenSolver {

	static final int NO_STATISTICS=0; // do not show statistics
	static final int SHOW_ONLY_TURNS=1; // show iterations
	static final int SHOW_ALL=2; // show all details

	/**
	 * Set numeric tolerance
	 * 
	 * @param prec
	 * 			the numeric tolerance
	 */
	public void setPrecision(double prec);

	/**
	 * Choose the information to compute and to show during the eigenvalue computation
	 */
	public void setVerbosity(int verbosity);

}
