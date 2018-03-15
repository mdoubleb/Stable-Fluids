package linalg;

/**
 * This class provides methods for printing vectors and matrices (with a given precision)
 * @author Luca Castelli Aleardi
 *
 */
public class MatrixUtilities {

	private static int numberPrecision=7; // number of digits to print (for a double)
	private static int vectorPrecision=4; // number of digits to print (for a vector)
	private static int matrixPrecision=3; // number of digits to print (for a matrix)
	
	/**
	 * Return an "approximation" (as String) of a given real number
	 */
	public static String approxNumber(double a) {
		String format="%."+numberPrecision+"f";
		String s=String.format(format,a);
		return s;
	}
	
	/**
	 * Print a (column) vector
	 */
	public static void print(double[] a) {
		for(int i=0;i<a.length;i++) {
			String format="%."+vectorPrecision+"f";
			String s=String.format(format,a[i]);
			System.out.println(s+"\t");
		}
		System.out.println();
	}
	
	/**
	 * Print a (column) vector
	 */
	public static void print(int[] a) {
		for(int i=0;i<a.length;i++) {
			System.out.println(a[i]+"\t");
		}
		System.out.println();
	}
	
	/**
	 * Print two column vectors (of same size)
	 */
	public static void print(double[] a, double[] b) {
		for(int i=0;i<a.length;i++) {
			String format="%."+vectorPrecision+"f";
			String sa=String.format(format,a[i]);
			String sb=String.format(format,b[i]);
			System.out.println(sa+"\t\t"+sb);
		}
		System.out.println();
	}

	/**
	 * Print a (dense) matrix
	 */
	public static void print(double[][] m) {
		for(int i=0;i<m.length;i++) {
			for(int j=0;j<m[i].length;j++) {
				String format="%."+matrixPrecision+"f";
				String s=String.format(format, m[i][j]);
				System.out.print(s+" \t");
			}
			System.out.println("");
		}
		System.out.println("");
	}

	/**
	 * Print the first k columns of a matrix
	 */
	public static void print(double[][] m, int k) {
		for(int i=0;i<m.length;i++) {
			for(int j=0;j<k;j++) {
				String format="%."+matrixPrecision+"f";
				String s=String.format(format, m[i][j]);
				System.out.print(s+" \t");
			}
			System.out.println("");
		}
		System.out.println("");
	}

}
