package linalg.MTJ;

import linalg.Matrix;

// import MTJ library for linear algebra computations
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

import java.util.*;

/**
 * @author Luca Castelli Aleardi
 * 
 * Wrapper class for sparse matrix implementations based on MTJ library
 */
public class MTJSparseMatrix implements Matrix {
	 
	/** matrix (sparse representation, based on MTJ library) */	
		public CompRowMatrix A;
		
		public MTJSparseMatrix(CompRowMatrix A) {
			this.A=A;
			
		}
		
		/**
		 * Computes the product of the matrix A and a given vector X
		 * 
		 * @param x
		 *            the input vector
		 * @return the product Ax
		 */
		public double[] times(double[] x) {
			DenseVector result=new DenseVector(x.length);
			A.mult(new DenseVector(x), result);
			return toArray(result);
		}
		
		/**
		 * Getter returning the height of the matrix
		 * @return the height of the matrix
		 */
		public int getHeight() {
			return this.A.numRows();
		}

		/**
		 * Getter returning the width of the matrix
		 * @return the width of the matrix
		 */
		public int getWidth() {
			return this.A.numColumns();
		}

		/**
		 * Convert the matrix to a double array
		 * 
		 * @return the double array corresponding to the matrix
		 */
		public double[][] toArray(){
			int n=A.numRows();
			int m=A.numColumns();
			double[] AA=A.getData();
			double[][] result=new double[n][m];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++)
					result[i][j]=AA[i+j*n];
			}
			return result;
		}
		
		/**
		 * Convert a dense vector to a one dimensional array
		 * 
		 * @return array corresponding to the array
		 */
		public static double[] toArray(DenseVector v){
			int n=v.size();
			double[] result=new double[n];
			for (int i = 0; i < n; i++)
					result[i]=v.get(i);
			return result;
		}
		
		/**
		 * Extract a column vector from the matrix
		 * 
		 * @param k
		 * 			the index of the column to extract
		 * @return
		 * 			one dimensional array corresponding to the k-th column
		 */
		public double[] extractColumn(int k){
			int n=A.numRows();
			int m=A.numColumns();
			if(k<0 || k>=m)
				throw new Error("Error: wrong column index "+k);
			
			double[] result=new double[n]; // create a one dimensional vector of size n
			for (int i = 0; i < n; i++)
					result[i]=A.get(i, k);
			
			return result;
		}
		
		/**
		 * Extract a row vector from the matrix
		 * 
		 * @param k
		 * 			the index of the column to extract
		 * @return
		 * 			one dimensional array corresponding to the k-th row
		 */
		public double[] extractRow(int k){
			int n=A.numRows();
			int m=A.numColumns();
			if(k<0 || k>=n)
				throw new Error("Error: wrong column index "+k);
			
			double[] result=new double[m]; // create a one dimensional vector of size n
			for (int j = 0; j < m; j++)
					result[j]=A.get(k, j);
			
			return result;
		}
		
		public String toString() {
			return this.A.toString();
		}

}
