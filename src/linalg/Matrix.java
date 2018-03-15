package linalg;

/**
 * Defines basic methods that any matrix representation must implement
 */
public interface Matrix {
	
	/**
	 * Computes the product of the matrix A and a given vector X
	 * 
	 * @param x
	 *            the input vector
	 * @return the product Ax
	 */
	public double[] times(double[] x);
	
	/**
	 * Getter returning the height of the matrix
	 * @return the height of the matrix
	 */
	public int getHeight();

	/**
	 * Getter returning the height of the matrix
	 * @return the width of the matrix
	 */
	public int getWidth();
}
