import linalg.MTJ.MTJSparseMatrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

/**
 * @author benjaminbarral
 * Solve matrices with Smooth Plastic Boundary conditions
 */
public class ComputeMatricesMTJ_BC2 {
	MTJSparseMatrix DiffusionX, Laplacien, DiffusionY;
	float dx,dy;
	int N;
	float dt;
	float visc;
	
	boolean dirichletCL=false;    // equivalent to periodic boundary conditions
	
	public ComputeMatricesMTJ_BC2(int n, float dx, float dy, float dt, float visc){
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
		CompRowMatrix L=new CompRowMatrix(N, N, mA); // create a sparse matrix of size nxn
		CompRowMatrix DX=new CompRowMatrix(N, N, mA);
		CompRowMatrix DY=new CompRowMatrix(N, N, mA);
		

		double Cx = -visc*dt/(dx*dx);
		double Cy = -visc*dt/(dy*dy);
		//Compute DX
		for(int i = 0; i<N;i++) {
			int k = i/n;
			int l = i%n;
			if(k!=0 && k!=n-1) DX.set(i,i,1-2*Cx -2*Cy);
			else DX.set(i,i,1-2*Cx-Cy);
			
			if(l!=0) DX.set(i,i-1,Cx);
			else DX.set(i,i+n-1,0.);
			if(l!=n-1) DX.set(i,i+1,Cx);
			else DX.set(i,i-n+1,0.);
			
			if(k!=0) DX.set(i,i-n,Cy);
			else DX.set(i,N-n+l,0.);
			if(k!=n-1) DX.set(i,i+n,Cy);
			else DX.set(i,l,0.);
		}
		
		//Compute DY
		for(int i = 0; i<N;i++) {
			int k = i/n;
			int l = i%n;
			if(l!=0 && l!=n-1) DY.set(i,i,1 -2*Cx -2*Cy);
			else DY.set(i,i,1 -Cx -2*Cy);
			
			if(l!=0) DY.set(i,i-1,Cx);
			else DY.set(i,i+n-1,0.);
			if(l!=n-1) DY.set(i,i+1,Cx);
			else DY.set(i,i-n+1,0.);

			if(k!=0) DY.set(i,i-n,Cy);
			else DY.set(i,N-n+l,0.);
			if(k!=n-1) DY.set(i,i+n,Cy);
			else DY.set(i,l,0.);
		}
		
		Cx = 1./(dx*dx);
		Cy = 1./(dy*dy);
		//Compute L
		for(int i = 0; i<N;i++) {
			int k = i/n;
			int l = i%n;
			L.set(i,i,-2*Cx -2*Cy);
			
			if(l!=0 && l!=n-1) { 
				L.set(i,i-1,Cx);
				L.set(i,i+1,Cx);
			}
			else {
				if(l==0) L.set(i,i+1,2*Cx);
				else L.set(i, i-n+1, 0.);
				if(l==n-1) L.set(i,i-1,2*Cx);
				else L.set(i, i+n-1, 0.);
			}
			
			if(k!=0 && k!=n-1) {L.set(i,i-n,Cy);
			L.set(i,i+n,Cy);}
			else {
				if(k!=n-1) L.set(i,i+n,Cy);
				else  L.set(i,l,0.);
				if(k!=0) L.set(i,i-n,Cy);
				else L.set(i,N-n+l,0.);
			}

		}

		DiffusionX = new MTJSparseMatrix(DX);
		Laplacien = new MTJSparseMatrix(L);
		DiffusionY = new MTJSparseMatrix(DY);
	}
	public static double[] DX(double[] U, int n, double dx, boolean tilde) {
		int N = n*n;
		double[] U2 = new double[N];
		for(int i = 0; i<U2.length;i++) {
			int l = i%n;
			int k = i/n;
			if (l!=0&& (!tilde || l!=n-1)) U2[i]-=U[i-1];   //AJOUTE
			if (l!=n-1 && (!tilde || l!=0)) U2[i]+=U[i+1];
			
			U2[i]=U2[i]/(2*dx);
		}
		return U2;
	}
	
	public static double[] DY(double[] U, int n, double dy, boolean tilde) {
		int N = n*n;
		double[] U2 = new double[N];
		for(int i = 0; i<U2.length;i++) {
			int l = i%n;
			int k = i/n;
			if (k!=0 && (!tilde || k!=n-1) ) U2[i]-=U[i-n];
			if (k!=n-1 && (!tilde || k!=0)) U2[i]+=U[i+n];
			
			U2[i]=U2[i]/(2*dy);
		}
		return U2;
	}
}
