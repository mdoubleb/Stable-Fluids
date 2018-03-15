/**
 * A class for representing Vectors in a d-dimensional space
 *
 * @author Luca Castelli Aleardi (INF555, 2012)
 */
public class Vector_d implements Vector_ {
	  public double[] coordinates;

	  public Vector_d() {}
	  
	  public Vector_d(double[] coord) { 
		  this.coordinates=new double[coord.length];
		  for(int i=0;i<coord.length;i++)
			  this.coordinates[i]=coord[i];
	  }

	  public Vector_d(Point_d a, Point_d b) { 
		  this.coordinates=new double[a.dimension()];
		  for(int i=0;i<a.dimension();i++)
			  this.coordinates[i]=b.getCartesian(i)-a.getCartesian(i);
	  }
	  
	  public float[] toFloat() {
		  float[] result=new float[this.dimension()];
		  for(int i=0;i<this.dimension();i++)
			  result[i]=(float)coordinates[i];
		  return result;
	  }
	  
	  public double[] toDouble() {
		  double[] result=new double[this.dimension()];
		  for(int i=0;i<this.dimension();i++)
			  result[i]=coordinates[i];
		  return result;
	  }

	  public String toString() {
		  String result="(";
		  for(int i=0;i<dimension()-1;i++)
			  result=result+this.getCartesian(i)+",";
		  return result+this.getCartesian(dimension()-1)+")";
	  }
	  
	  public int dimension() { return this.coordinates.length;}
	  
	  public double getCartesian(int i) {
		  return this.coordinates[i];
	  }
	  public void setCartesian(int i, double x) {
		  this.coordinates[i]=x;
	  }
	    
	  public Vector_d sum(Vector_ v) {
		  throw new Error("a' completer");
	  }
	  
	  public Vector_d difference(Vector_ v) {
		  throw new Error("a' completer");
	  }
	  
	  public Vector_d opposite() {
		  throw new Error("a' completer");
	  }
	  
	  public double innerProduct(Vector_ v) {
		  throw new Error("a' completer");
	  }

	  public Vector_d divisionByScalar(double s) {
	  	if(s==0.) throw new Error("error: division by zero");
		  throw new Error("a' completer");
	  }
	  
	  public Vector_d multiplyByScalar(double s) {
		  throw new Error("a' completer");
	  }
	  
	  public double squaredLength() {
		  throw new Error("a' completer");
	  }
	  
	  public Vector_d crossProduct(Vector_ b) {
		  throw new Error("a' completer");
	  }
	  
	  public boolean equals(Object o) {
		  Point_d p = (Point_d) o;
		  for(int i=0;i<dimension();i++) {
			  if(this.coordinates[i]!=p.getCartesian(i))
				  return false;
		  }
		  return true;
	  }


}
