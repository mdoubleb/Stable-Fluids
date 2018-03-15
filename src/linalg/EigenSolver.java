package linalg;

/**
 * @author Luca Castelli Aleardi
 * 
 * Interface providing methods for for computing the Eigenvalue Decomposition of a Laplacian matrix
 */
public interface EigenSolver {

	/**
	 * Compute the first 'k' eigenvalues and corresponding eigenvectors
	 * 
	 * @param k  number of eigenvalues to compute
	 */
	public void computeEigenvalueDecomposition(int k);
	
	/**
	 * Return the eigenvalues
	 * 
	 * @return  an array storing the eigenvalues
	 */
	public double[] getEigenvalues();
	
	/**
	 * Return the eigenvectors
	 * 
	 * @return  an array storing the eigenvector
	 */
	public double[][] getEigenvectors();

}
