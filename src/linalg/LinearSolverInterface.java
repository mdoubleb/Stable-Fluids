package linalg;

public interface LinearSolverInterface {
	public int getHeight();
	
	/** Return the input matrix (representing the linear system to solve) */
	public Matrix getMatrix();

	/** Return the solution of Ax=b */
	public double[] solve(double[] b);

	/** Return the product Ax */
	public double[] times(double[] x);

}
