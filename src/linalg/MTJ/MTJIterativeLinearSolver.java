package linalg.MTJ;


import linalg.LinearIterativeSolver;
import linalg.Matrix;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.AbstractIterationMonitor;
import no.uib.cipr.matrix.sparse.CG;
import no.uib.cipr.matrix.sparse.DefaultIterationMonitor;
import no.uib.cipr.matrix.sparse.IterationMonitor;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;

/**
 * @author Luca Castelli Aleardi
 * 
 * Wrapper class for the Conjugate Gradient Solver, based on MTJ library (using sparse matrices)
 */
public class MTJIterativeLinearSolver implements LinearIterativeSolver {

	MTJSparseMatrix laplacian; // laplacian matrix of the graph (based on MTJ library)
	double[] start; // starting vector (initial guess)
	double precision; // numeric tolerance
	DenseVector startingVector;
	
	IterativeSolver itSolver; // iterative solver implemented in MTJ library (Conjugate Gradient)

	private int verbosity=0; // decide whether to show convergence statistics
	private DefaultIterationMonitor m;

	/**
	 * Initialize the Conjugate Gradient solver of MTJ library
	 * 
	 * @param A
	 * 			the input matrix
	 */	
	public MTJIterativeLinearSolver(MTJSparseMatrix A, double precision) {
		int n=A.getHeight(); // matrix size
		this.start=new double[n]; // initial guess for the CG iterator
		this.laplacian=A; // create a sparse matrix of size nxn
		this.precision=precision;
		
		// initialize the starting vector (initial guess)
		for (int i=0; i<n; i++) {
			this.start[i]=Math.random()*10; 
		}
		this.startingVector = new DenseVector(start); 

		this.itSolver=new CG(this.startingVector); // set iterative solver
		//itSolver.setIterationMonitor(new SimpleIterationMonitor(500)); // set number of iterations
	
		this.m=new DefaultIterationMonitor();
		m.setMaxIterations(500);
		this.m.setRelativeTolerance(precision);
		this.m.setNormType(Vector.Norm.Two);
		this.m.setAbsoluteTolerance(precision);
		itSolver.setIterationMonitor((IterationMonitor) m); // set number of iterations

	}
	
	/**
	 * Solve linear system Lx=b (where L is the graph laplacian)
	 * 
	 * @param b right hand side vector
	 * 
	 * @return the vector solution x[]
	 */
	public double[] solve(double[] b) {
		// Solve using CG
		//itSolver.setIterationMonitor(new SimpleIterationMonitor( )); // set number of iterations
		DenseVector v;
		v=new DenseVector(b);
		try {
			itSolver.solve(this.laplacian.A, v, v);
		}
		catch (IterativeSolverNotConvergedException e){}
		
		if(this.verbosity>0) {
			int iter=m.iterations();
			double residual=m.residual();
			System.err.println("  Conjugate Gradient solved in "+iter+" turns, dist = "+residual);
		}

		double[] result=MTJSparseMatrix.toArray(v);
		return result;
	}
	
	/** 
	 * Return the product Ax
	 * */
	public double[] times(double[] x) {
		return this.laplacian.times(x);
	}
	
	/** 
	 * Run the iterative solver k times (solve the system Ax=b)
	 * 
	 * @param b
	 * 			right side vector
	 * @param k
	 * 			number of iterations to perform
	 * 
	 * @return	the solution of Ax=b
	 * 
	 * */
	public double[] runKtimes(double[] b,int k) {
		throw new Error("To be implemented");
	}
	
	/**
	 * Set the starting vector (initial guess for the CG method)
	 * 
	 * @param x starting vector
	 */
	public void setStartingVector(double[] start) {
		this.start=start.clone();
	}
	
	public void setPrecision(double precision){
		this.precision = precision;
		//this.m.setRelativeTolerance(precision);
	}
	
	public double getPrecision(){
		return this.precision;
	}
	
	/**
	 * Getter returning the height of the matrix
	 * @return the height of the matrix
	 */
	public int getHeight() {
		return this.laplacian.getHeight();
	}
	
	/**
	 * Return the input matrix (representing the linear system to solve)
	 */
	public Matrix getMatrix() {
		return this.laplacian;
	}
	
	public void setVerbosity(int verbosity){
		this.verbosity=verbosity;
	}
	
	public String toString() {
		return this.laplacian.toString();
	}
	
	/** 
	 * Check whether Ax=b
	 * 
	 *  @param b
	 *  		right side vector
	 *  @param x
	 *  		solution to check
	 *  @param prec
	 *  		numeric tolerance
	 *  
	 *  @return true if dist(Ax, b) < precision
	 **/
	public boolean checkSolution(double[] x, double[] b, double prec) {
		throw new Error("To be completed");
	}

	public static class SimpleIterationMonitor extends AbstractIterationMonitor {
	    private int max;
	 
	    SimpleIterationMonitor(int max) {
	       this.max = max;
	       System.out.println(""+this.getNormType());
	     }
	    protected boolean convergedI(double r, Vector x) throws IterativeSolverNotConvergedException {
	       return convergedI(r);
	     }
	    protected boolean convergedI(double r) throws IterativeSolverNotConvergedException {
	       return iter >= max;
	     }
	}

}
