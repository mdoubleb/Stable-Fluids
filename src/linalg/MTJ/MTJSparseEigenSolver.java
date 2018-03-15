package linalg.MTJ;

import linalg.IterativeEigenSolver;
import no.uib.cipr.matrix.DenseMatrix; // MTJ library
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import sparse.eigenvolvers.java.*; // SEJ library

/**
 * @author Luca Castelli Aleardi
 * 
 * Eigenvalue decomposition of a Laplacian matrix, based on MTJ/SEJ libraries (with sparse matrices)
 * It makes use of conjugate gradient preconditioner and Compressed Row Representation for sparse matrices
 */
public class MTJSparseEigenSolver implements IterativeEigenSolver {
		MTJSparseMatrix laplacian; // sparse matrix (based on MTJ library)
	 
		private boolean withPreconditioning=false; // run with or without preconditioner
		private int blockSize; // number of eigenvectors to compute
		private int maxIterations=2000; // maximum number of iterations
		private double tolerance=1e-5; // numeric tolerance
		private int verbosity=0; // for showing statistics
		DenseMatrix blockVectorX;
		int n;
		
		Lobpcg decomposition; // eigenvalue decomposition (based on SEJ library)
		
		/**
		 * Initialize the matrix
		 */	
		public MTJSparseEigenSolver(MTJSparseMatrix m) {
			n=m.getHeight(); // matrix size
			this.laplacian=m; // sparse matrix of size nxn
			
			this.decomposition=new Lobpcg(); // initialize eigensolver
		}
		
		/**
		 * Choose the information to compute and to show during the eigenvalue computation
		 */
		public void setVerbosity(int verbosity) {
			this.verbosity = verbosity;
		}
		
		/**
		 * Set numeric tolerance
		 * 
		 * @param prec
		 * 			the numeric tolerance
		 */
		public void setPrecision(double prec) {
			this.tolerance=prec;
		}

		public void computeEigenvalueDecomposition(int k) {
			System.out.print("Computing eigenvalue decomposition... ");
			blockVectorX=(DenseMatrix) Matrices.random(n,k);
			
			long startTime=System.nanoTime(), endTime; // for evaluating time performances
			
			if(this.withPreconditioning==false) // run without preconditioning
				this.decomposition.runLobpcg(blockVectorX, this.laplacian.A,tolerance,maxIterations,verbosity); // run without preconditioner
			else { // run with preconditioning
				OperatorPrecCG operT= new OperatorPrecCG(this.laplacian.A);
				int innerIterations=1;
				operT.setCGNumberIterations(innerIterations); // Need to try different values to see best convergence
				this.decomposition.runLobpcg(blockVectorX,this.laplacian.A,null,operT,tolerance,maxIterations,verbosity); // LCA
			}
			
			endTime=System.nanoTime();
	        double duration=(double)(endTime-startTime)/1000000000.;

	        System.out.println("done ("+duration+" seconds)");
	        System.out.print("number of computed eigenvectors: "+blockVectorX.numColumns());
	        System.out.println("\t numeric precision: "+tolerance);
		}
		
		public double[] getEigenvalues() {
			double[] result=this.decomposition.getEigenvalues();
			return result;
		}
		
		public double[][] getEigenvectors() {
			Matrix m=this.decomposition.getEigenvectors();
			int nRows=m.numRows();
			int nColumns=m.numColumns();
			double[][] result=new double[nColumns][nRows]; // the transpose of m
			for(int i=0;i<nColumns;i++) { // compute the transpose
				for(int j=0;j<nRows;j++) {
					result[i][j]=m.get(j, i);
				}
			}
			return result;
		}
		
		private double[][] toArray(Matrix m) {
			int nRows=m.numRows();
			int nColumns=m.numColumns();
			double[][] result=new double[nRows][nColumns];
			for(int i=0;i<nRows;i++) {
				for(int j=0;j<nColumns;j++) {
					result[i][j]=m.get(i, j);
				}
			}
			return result;
		}
		
		public String toString() {
			return this.laplacian.toString();
		}

}
