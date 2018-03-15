
import linalg.MTJ.MTJIterativeLinearSolver;


/**
 * @author benjaminbarral and mikhailbessa
 * Mixed Solver for the velocity field and the scalar field (dust density)
 */
public class ScalarSolver extends AbstractSolver {
	
	int n,N;
	float dx,dy;
	float dt;
	float visc; //viscosity
	float kappa;  //diffusivity constant
	double precisionDiff; //precision for the diffusion step
	double precisionDiss; //precision for the dissipation step
	double[][] VX,VY; //velocities
	double[] Q;     //for projection step
	double[] FX,FY;//force
	double[] source;
	double[][] scalar;
	Grid displayGrid;
	int counter; //number of steps
	int BC; //boundary conditions
	MTJIterativeLinearSolver solverDiffX, solverDiffY, solverDiff, solverDiss, solverDiffScalar;
	int scalarColor;
	
	public ScalarSolver (double[] initialVelocityX, double[] initialVelocityY, int n, float dx, float dy, float dt, float visc, float kappa, Grid grid, double precision,int bC, int scalarColor){
		this.n=n;
		N=n*n;
		this.dx=dx;
		this.dy=dy;
		this.dt=dt;
		this.BC=bC;
		this.visc=visc;
		this.kappa=kappa;
		this.scalarColor=scalarColor;
		this.displayGrid=grid;
		this.precisionDiff = Math.pow(10, -15);
		this.precisionDiss = 0.0001;
		VX=new double[2][N];
		VY=new double[2][N];
		for (int i=0; i<N; i++){
			VX[0][i]=initialVelocityX[i];
			VY[0][i]=initialVelocityY[i];
		}
		scalar = new double[2][N];
		source= new double[N];
		Q=new double[N];
		if(BC==0) {
			ComputeMatricesMTJ m = new ComputeMatricesMTJ(n, dx, dy, dt, visc);
			solverDiff = new MTJIterativeLinearSolver(m.Diffusion, precisionDiff);
			ComputeMatricesMTJ mScalar = new ComputeMatricesMTJ(n, dx, dy, dt, kappa);
			solverDiffScalar = new MTJIterativeLinearSolver(mScalar.Diffusion, precisionDiff);
			solverDiss = new MTJIterativeLinearSolver(m.Laplacien, 0.000001);
		}
		else if(BC==1) {
			ComputeMatricesMTJ_BC2_Scalar m = new ComputeMatricesMTJ_BC2_Scalar(n, dx, dy, dt, visc,kappa);
			solverDiffX = new MTJIterativeLinearSolver(m.DiffusionX, precisionDiff);
			solverDiss = new MTJIterativeLinearSolver(m.Laplacien, 0.000001);
			solverDiffY = new MTJIterativeLinearSolver(m.DiffusionY, precisionDiff);
			solverDiffScalar = new MTJIterativeLinearSolver(m.DiffusionScalar, precisionDiff);
		}
		counter=0;
	}
	
	public void routine(){
		addForce(FX,FY);
		addSource(source);
		Transport();
		swapVelocities();
		swapScalar();
		Diffuse();
		DiffuseScalar();
		Project();
		counter++;
	}

	public void addForce(double[] FX, double[] FY) {
		VX[0]=plus(VX[0],times(FX,dt));
		VY[0]=plus(VY[0],times(FY,dt));
	}
	
	public void addSource(double[] source) {
		scalar[0]=plus(scalar[0],times(source,dt));
		/*for (int i=0; i<N; i++){
			if (scalar[0][i]>255){
				scalar[0][i]=255;
			}
		}*/
	}
	
	public void Transport(){
		Point_2 p;
		int[][][] textures0=displayGrid.textures.clone();
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				p=Trace(i,j);
				linearInterpolate(p.getX(),p.getY(), i,j, true, textures0);
			}
		}
	}
	
	Point_2 Trace(int i, int j){
		double x=(j+0.5)*dx;
		double y=(i+0.5)*dy;
		double ux=VX[0][i*n+j];
		double uy=VY[0][i*n+j];
		x=x-dt*ux;
		y=y-dt*uy;
		return new Point_2(x,y);
	}
	
	void linearInterpolate(double x, double y, int i0, int j0, boolean setTexture, int[][][] textures0){
		//int j=(int)((x-dx/2)/dx);
		//double t1=(x-(j+0.5)*dx)/dx;
		int j=(x>0)?(int)((x-dx/2)/dx) : (int)((x-dx/2)/dx)-1;
		double t1=(x-(j+0.5)*dx)/dx;
		if (j>=n){
			j=j%n;
			if (BC==1) return;    
		}
		if (j<0){
			if (BC==1){
				if (j==-1 && t1>0.5) t1=1.;
				else return;
			}
			int q=j/n;
			j=(j+(-q+1)*n)%n;
		}
		t1=(j==n-1 && BC==1)? 0. : t1;
		//int i=(int)((y-dy/2)/dy);
		//double t=(y-(i+0.5)*dy)/dy;
		int i=(y>0)?(int)((y-dy/2)/dy) : (int)((y-dy/2)/dy)-1;
		double t=(y-(i+0.5)*dy)/dy;
		if (i>=n){
			i=i%n;
			if (BC==1) return;
		}
		if (i<0){
			if (BC==1){
				if (i==-1 && t>0.5) t=1.;
				else return;
			}
			int q=i/n;
			i=(i+(-q+1)*n)%n;
		}
		t=(i==n-1 && BC==1)? 0. : t;
		int i2=(i==n-1)? 0:i+1;
		int j2=(j==n-1)? 0:j+1;
		double s1a=VX[0][i*n+j],s1b=VX[0][i*n+j2];
		double s2a=VX[0][i2*n+j],s2b=VX[0][i2*n+j2];
		double s1=t1*s1b+(1-t1)*s1a;
		double s2=t1*s2b+(1-t1)*s2a;
		double uxBis=t*s2+(1-t)*s1;
		//double uxBis=VX[0][i*n+j];
		VX[1][i0*n+j0]= uxBis;
		s1a=VY[0][i*n+j];
		s1b=VY[0][i*n+j2];
		s2a=VY[0][i2*n+j];
		s2b=VY[0][i2*n+j2];
		s1=t1*s1b+(1-t1)*s1a;
		s2=t1*s2b+(1-t1)*s2a;
		double uyBis=t*s2+(1-t)*s1;
		//double uyBis=UY[i*n+j];
		VY[1][i0*n+j0]= uyBis;
		
		s1a=scalar[0][i*n+j];
		s1b=scalar[0][i*n+j2];
		s2a=scalar[0][i2*n+j];
		s2b=scalar[0][i2*n+j2];
		s1=t1*s1b+(1-t1)*s1a;
		s2=t1*s2b+(1-t1)*s2a;
		double scalarBis=t*s2+(1-t)*s1;
		//double uyBis=UY[i*n+j];
		scalar[1][i0*n+j0]= scalarBis;
		
		if (setTexture){
			for (int k=0; k<3; k++){
				displayGrid.setTexture(i0, j0, k, scalarColor);
				displayGrid.setTransparency(i0, j0, Math.min(255, (int)scalarBis));
			}
		}
	}
	
	
	public void Diffuse(){
		//System.out.println("Diffusing...");
		if(BC==0) {
			if (!equalNull(VX[0])){
				VX[0]=solverDiff.solve(VX[0]);
			}
			if (!equalNull(VY[0])){
				VY[0]=solverDiff.solve(VY[0]);
			}
		}
		else if(BC==1){
			VX[0]=solverDiffX.solve(VX[0]);
			VY[0]=solverDiffY.solve(VY[0]);	
		}
		//System.out.println("Diffusing done.");
	}
	
	public void DiffuseScalar(){
		//System.out.println("Diffusing (scalar) ...");
		if(BC==0) {
			if (!equalNull(scalar[0])){
				scalar[0]=solverDiffScalar.solve(scalar[0]);
			}
		}
		else if(BC==1){
			scalar[0]=solverDiffScalar.solve(scalar[0]);
		}
		//System.out.println("Diffusing (scalar) done.");
	}
	
	public boolean equalNull(double[] tab){
		int i=0;
		while (i<tab.length){
			if (Math.abs(tab[i])>Math.pow(10,-15))
				return false;
			i++;
		}
		return true;
	}
	
	public void Project(){
		//System.out.println("Dissipating...");
		if(BC==0) {
			double[] q = plus(ComputeMatricesMTJ.DX(VX[0], n, dx),ComputeMatricesMTJ.DY(VY[0],n,dy));
			if (!equalNull(q))
				Q=solverDiss.solve(q);
			else
				Q=q;
			VX[0]=minus(VX[0],ComputeMatricesMTJ.DX(Q,n,dx));
			VY[0]=minus(VY[0],ComputeMatricesMTJ.DY(Q,n,dy));
		}
		else if(BC==1) {
			double[] q = plus(ComputeMatricesMTJ_BC2.DX(VX[0], n, dx, false),ComputeMatricesMTJ_BC2.DY(VY[0],n,dy, false));
			if (!equalNull(q))
				Q=solverDiss.solve(q);
			else
				Q=q;
			VX[0]=minus(VX[0],ComputeMatricesMTJ_BC2.DX(Q,n,dx,true));
			VY[0]=minus(VY[0],ComputeMatricesMTJ_BC2.DY(Q,n,dy,true));
		}
		//System.out.println("Dissipating done.");
	}
	
	public void setForce(double[] Fx,double[] Fy){
		this.FX=Fx;
		this.FY=Fy;
	}
	
	
	

	public void printVelocity(){
		System.out.println("Step number "+counter + " :");
		System.out.print("Ux= ( ");
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				System.out.print(VX[0][i*n+j]+ " ");
			}
			if (i==n-1){
				System.out.println(")");
			}
			else
				System.out.println("");
		}
		System.out.print("Uy= ( ");
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				//System.out.print(VY[0][i*n+j]+ " ");
			}
			if (i==n-1){
				//System.out.println(")");
			}
			else {
				//System.out.println("");
			}
		}
	}

	public void swapVelocities (){
		double[] Mint= VX[0];
		VX[0]=VX[1];
		VX[1]=Mint;
		double[] Mint2= VY[0];
		VY[0]=VY[1];
		VY[1]=Mint2;
	}
	
	public void swapScalar() {
		double[] Mint= scalar[0];
		scalar[0]=scalar[1];
		scalar[1]=Mint;
	}
	public double[] plus (double[] array1, double[] array2){
		double[]res=new double[array1.length];
		for (int i=0; i<array1.length; i++){
			res[i]=array1[i]+array2[i];
		}
		return res;
	}
	
	public double[] minus (double[] array1, double[] array2){
		double[]res=new double[array1.length];
		for (int i=0; i<array1.length; i++){
			res[i]=array1[i]-array2[i];
		}
		return res;
	}
	
	public double[] times(double[] array, double x){
		double[] res=new double[array.length];
		for (int i=0; i<array.length; i++){
			res[i]=array[i]*x;
		}
		return res;
	}

	@Override
	public void setSource(double[] source) {
		this.source= source;
		
	}
}
