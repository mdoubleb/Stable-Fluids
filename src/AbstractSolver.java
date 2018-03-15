
/**
 * 
 * @author benjaminbarral
 * Class for the fluid solver : can be fluid motion solver OR mixed substance introduction+velocity solver
 */
public abstract class AbstractSolver {
	public abstract void routine();
	public abstract void setForce(double[] FX, double[] FY);
	public abstract void setSource(double[] source);
	public void changeColor() {
		// TODO Auto-generated method stub
		
	}
}
