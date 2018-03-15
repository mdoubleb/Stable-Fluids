/**
 * A class for representing Points in 3D (with floating points coordinates)
 *
 * @author Luca Castelli Aleardi (INF555, 2012)
 */
public class Point_3 extends Point_2{
  public double z;

  public Point_3() {}
  
  public Point_3(double x,double y,double z) { 
  	this.x=x; 
  	this.y=y; 
  	this.z=z;
  }
  
  public Point_3(double[] coordinates) { 
	  	this.x=coordinates[0]; 
	  	this.y=coordinates[1]; 
	  	this.z=coordinates[2];
  }

  /**
   * Construct a new point (copying point p)
   */
  public Point_3(Point_3 p) { 
  	this.x=p.getCartesian(0); 
  	this.y=p.getCartesian(1); 
  	this.z=p.getCartesian(2);
  }
  
  public void barycenter(Point_ [] points) {
  	double x_=0., y_=0., z_=0;
  	for(int i=0;i<points.length;i++) {
  		x_ += points[i].getCartesian(0);
  		y_ += points[i].getCartesian(1);
  		z_ += points[i].getCartesian(2);
  	}
  	this.x = x_/points.length;
  	this.y = y_/points.length;
  	this.z = z_/points.length;
  }

  /**
   * Return a linear combination of a set of points.
   * Coefficients must sum to 1
   */
  public void linearCombination(Point_[] points, double[] coefficients) {
	  	double x_=0., y_=0., z_=0.;
	  	for(int i=0;i<points.length;i++) {
	  		x_=x_+(points[i].getCartesian(0)*coefficients[i]);
	  		y_=y_+(points[i].getCartesian(1)*coefficients[i]);
	  		z_=z_+(points[i].getCartesian(2)*coefficients[i]);
	  	}
	  	this.x = x_;
	  	this.y = y_;
	  	this.z = z_;
	  }

  /**
   * Return a linear combination of a set of points.
   * Coefficients must sum to 1
   */
  public static Point_3 linearCombination(Point_3 [] points, double[] coefficients) {
  	double x_=0., y_=0., z_=0;
  	for(int i=0;i<points.length;i++) {
  		x_=x_+(points[i].getX()*coefficients[i]);
  		y_=y_+(points[i].getY()*coefficients[i]);
  		z_=z_+(points[i].getZ()*coefficients[i]);
  	}
  	return new Point_3(x_,y_,z_);
  }
  
  /**
   * Return the z coordinate of the point
   */
  public double getZ() {return z; }

  /**
   * Set the z coordinate 
   */
  public void setZ(double z) {this.z=z; }
  
  public float[] toFloat() {
	  float[] result=new float[3];
	  result[0]=(float)x;
	  result[1]=(float)y;
	  result[2]=(float)z;
	  return result;
  }
  
  public double[] toDouble() {
	  double[] result=new double[3];
	  result[0]=x;
	  result[1]=y;
	  result[2]=z;
	  return result;
  }

  /**
   * Return a 2D point (with cartesian coordinates).
   * The current point lives in 2D (with homogeneous coordinates)
   */
  public Point_ toCartesian() {
	  if(z==0.) throw new Error("Error: division by 0");
	  return new Point_2(x/z, y/z);
  }

  /**
   * Return the point in homogeneous coordinates (point in 4D)
   */
  public Point_ toHomogeneous() {
	  double[] hCoordinates=new double[4];
	  hCoordinates[0]=x;
	  hCoordinates[1]=y;
	  hCoordinates[2]=z;
	  hCoordinates[3]=1.;
	  return new Point_d(hCoordinates);
  }

  public void translateOf(Vector_ v) {
    this.x=x+v.getCartesian(0);
    this.y=y+v.getCartesian(1);
    this.z=z+v.getCartesian(2);
  }

  public void multiply (double n) {
	    this.x*=n;
	    this.y*=n;
	    this.z*=n;
	  }

  public boolean equals(Object o) { 
	  Point_ p = (Point_) o;
	  return this.x==p.getCartesian(0) && this.y==p.getCartesian(1) && 
	  this.z==p.getCartesian(2); 
  }

  /**
   * Return the Euclidean distance of two points
   */
  public double distanceFrom(Point_3 p) {
    double dX=p.getX()-x;
    double dY=p.getY()-y;
    double dZ=p.getZ()-z;
    return Math.sqrt(dX*dX+dY*dY+dZ*dZ);
  }
  
  public double squareDistance(Point_3 p) {
    double dX=p.getX()-x;
    double dY=p.getY()-y;
    double dZ=p.getZ()-z;
    return dX*dX+dY*dY+dZ*dZ;
  }

  public String toString() {return "("+x+","+y+","+z+")"; }
  public int dimension() { return 3;}
  
  public double getCartesian(int i) {
  	if(i==0) return x;
  	else if(i==1) return y;
  	return z;
  } 
  public void setCartesian(int i,double x) {
  	if(i==0) this.x=x;
  	else if(i==1) this.y=x;
  	else this.z=x;
  }

  public void setOrigin() {
	  this.x=0.;
	  this.y=0.;
	  this.z=0.;
  }

  public Vector_ minus(Point_ b){
	  return new Vector_3(b.getCartesian(0)-x, b.getCartesian(1)-y, b.getCartesian(2)-z);
  }
  
  public int compareTo(Point_ o) {
	  Point_3 p = (Point_3) o;
	  if (this.getX() < p.getX())
		  return -1;
	  if (this.getX() > p.getX())
		  return 1;
	  if (this.getY() < p.getY())
		  return -1;
	  if (this.getY() > p.getY())
		  return 1;
	  if (this.getZ() < p.getZ())
		  return -1;
	  if (this.getZ() > p.getZ())
		  return 1;
	  return 0;
  }

  
}




