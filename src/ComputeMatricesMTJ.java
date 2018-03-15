import linalg.MTJ.MTJSparseMatrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

/**
 * @author benjaminbarral
 * Solve matrices with Periodic Boundary conditions
 */

public class ComputeMatricesMTJ {
	MTJSparseMatrix Laplacien, Diffusion;
	float dx,dy;
	int N;
	float dt;
	float visc;
	
	public ComputeMatricesMTJ(int n, float dx, float dy, float dt, float visc){
		N=n*n;
		this.dx=dx;
		this.dy=dy;
		this.dt=dt;
		this.visc=visc;
		
		int[][] mA=new int[N][5]; // indices of non-zero entries for each row		
		for(int i = 0; i<N;i++) {	
			int[] m = new int[5];
			m[0]=i;
			int k = i/n;
			int l = i%n;
			if(l!=0) m[1]=i-1;
			else m[1]=i+n-1;
			if(l!=n-1) m[2]= i+1;
			else m[2]=i-n+1;
			
			if(k!=0) m[3]=i-n;
			else m[3]=N-n+l;
			if(k!=n-1)m[4]=i+n;
			else m[4]=l;
			mA[i]=m;
		}
		CompRowMatrix A=new CompRowMatrix(N, N, mA); // create a sparse matrix of size nxn
		CompRowMatrix D=new CompRowMatrix(N, N, mA.clone());
		

		double Cx = 1./(dx*dx);
		double Cy = 1./(dy*dy);
		//Compute A
		for(int i = 0; i<N;i++) {
			int k = i/n;
			int l = i%n;
			A.set(i,i,-2*Cx -2*Cy);
			
			if(l!=0) A.set(i,i-1,Cx);
			else A.set(i,i+n-1,Cx);
			if(l!=n-1) A.set(i,i+1,Cx);
			else A.set(i,i-n+1,Cx);
			
			if(k!=0) A.set(i,i-n,Cy);
			else A.set(i,N-n+l,Cy);
			if(k!=n-1) A.set(i,i+n,Cy);
			else A.set(i,l,Cy);
		}
		
		Cx = -visc*dt/(dx*dx);
		Cy = -visc*dt/(dy*dy);
		//Compute D
		for(int i = 0; i<N;i++) {
			int k = i/n;
			int l = i%n;
			D.set(i,i,1 -2*Cx -2*Cy);
			
			if(l!=0) D.set(i,i-1,Cx);
			else D.set(i,i+n-1,Cx);
			if(l!=n-1) D.set(i,i+1,Cx);
			else D.set(i,i-n+1,Cx);
			
			if(k!=0) D.set(i,i-n,Cy);
			else D.set(i,N-n+l,Cy);
			if(k!=n-1) D.set(i,i+n,Cy);
			else D.set(i,l,Cy);
		}
		Laplacien = new MTJSparseMatrix(A);
		Diffusion = new MTJSparseMatrix(D);
	}
	public static double[] DX(double[] U, int n, double dx) {
		int N = n*n;
		double[] U2 = new double[N];
		for(int i = 0; i<U2.length;i++) {
			int l = i%n;
			if (l!=0) U2[i]-=U[i-1];
			else U2[i]-=U[i+n-1];
			if (l!=n-1) U2[i]+=U[i+1];
			else U2[i]+=U[i-n+1];
			
			U2[i]=U2[i]/(2*dx);
		}
		return U2;
	}
	
	public static double[] DY(double[] U, int n, double dy) {
		int N = n*n;
		double[] U2 = new double[N];
		for(int i = 0; i<U2.length;i++) {
			int l = i%n;
			int k = i/n;
			if (k!=0) U2[i]-=U[i-n];
			else U2[i]-=U[N-n+l];
			if (k!=n-1) U2[i]+=U[i+n];
			else U2[i]+=U[l];
			
			U2[i]=U2[i]/(2*dy);
		}
		return U2;
	}
}
