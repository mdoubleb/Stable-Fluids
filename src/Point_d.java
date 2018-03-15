/**
 * A class for representing Points in a d-dimensional space
 *
 * @author Luca Castelli Aleardi (INF555, 2012)
 */
public class Point_d implements Point_ {
	  public double[] coordinates;

	  public Point_d() {}

	  public Point_d(int d) {
		  this.coordinates=new double[d];
	  }
	  
	  public Point_d(double[] coord) { 
		  this.coordinates=new double[coord.length];
		  for(int i=0;i<coord.length;i++)
			  this.coordinates[i]=coord[i];
	  }

	  /**
	   * Construct a new point (copying point p)
	   */
	  public Point_d(Point_d p) { 
		  this.coordinates=new double[p.dimension()];
		  for(int i=0;i<p.dimension();i++)
			  this.coordinates[i]=p.getCartesian(i);
	  }

	  public double[] toDouble() {
		  throw new Error("a' completer");
	  }

	  public float[] toFloat() {
		  throw new Error("a' completer");
	  }
	  
	  public Point_ toCartesian() {
		  if(this.coordinates[this.dimension()-1]==0.)
			  throw new Error("Error: division by 0.");
		  double[] result=new double[dimension()-1];
		  for(int i=0;i<dimension()-1;i++)
			  result[i]=getCartesian(i)/getCartesian(this.dimension()-1);
		  
		  return new Point_d(result);
	  }

	  public Point_ toHomogeneous() {
		  throw new Error("method not defined");
	  }

	  public void barycenter(Point_[] points) {
	  	double[] x=new double[dimension()];
	  	for(int i=0;i<points.length;i++) {
	  		for(int j=0;j<dimension();j++)
	  			x[j]=x[j]+points[i].getCartesian(j);
	  	}
	  	for(int j=0;j<dimension();j++)
	  		this.coordinates[j]=x[j]/points.length;
	  }

	  /**
	   * Return a linear combination of a set of points.
	   * Coefficients must sum to 1
	   */
	  public void linearCombination(Point_d[] points, double[] coefficients) {
		  throw new Error("to be completed");
	  }
	    
	  public void translateOf(Vector_ v) {
		  if(v.dimension()!=this.dimension())
			  throw new Error("Point_d error: wrong dimension");
		  for(int i=0;i<dimension();i++)
			  this.coordinates[i]=this.coordinates[i]+v.getCartesian(i);
	  }

	  public boolean equals(Object o) {
		  Point_d p = (Point_d) o;
		  for(int i=0;i<dimension();i++) {
			  if(this.coordinates[i]!=p.getCartesian(i))
				  return false;
		  }
		  return true;
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

	  public void setOrigin() {
		  for(int i=0;i<this.coordinates.length;i++)
			  this.coordinates[i]=0.;
	  }
	    
	  public Vector_ minus(Point_ b){
		  throw new Error("a' completer");
	  }

	    public double squareDistance(Point_d p) {
	    	double d = 0.;
	    	
	    	int dim = Math.min (this.dimension(), p.dimension());
	    	for (int i=0; i<dim; i++)
	    	    d += (this.getCartesian(i)-p.getCartesian(i)) * 
	    	    (this.getCartesian(i)-p.getCartesian(i));
	    	return d;
	    }

	  public int compareTo(Point_ o) {
		  throw new Error("a' completer");
	  }
	  
}
