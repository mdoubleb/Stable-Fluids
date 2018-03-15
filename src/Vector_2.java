public class Vector_2 implements Vector_{
  public double x,y;

  public Vector_2() {}
  
  public Vector_2(double x,double y) { 
  	this.x=x; 
  	this.y=y;
  }
  
  public Vector_2(double[] coordinates) { 
	  	this.x=coordinates[0]; 
	  	this.y=coordinates[1];
  }

  public Vector_2(Point_2 a, Point_2 b) { 
	  	this.x=b.getX()-a.getX(); 
	  	this.y=b.getY()-a.getY(); 
  }

  public double getX() {return x; }
  public double getY() {return y; }
  
  public void setX(double x) {this.x=x; }
  public void setY(double y) {this.y=y; }
  
  public float[] toFloat() {
	  float[] result=new float[2];
	  result[0]=(float)x;
	  result[1]=(float)y;
	  return result;
  }
  
  public double[] toDouble() {
	  double[] result=new double[2];
	  result[0]=x;
	  result[1]=y;
	  return result;
  }  

  public boolean equals(Vector_ v) { 
  	return (this.x==v.getCartesian(0) && this.y==v.getCartesian(1)); 
  }

  public String toString() {return "["+x+","+y+"]"; }
  public int dimension() { return 2;}
  
  public double getCartesian(int i) {
  	if(i==0) return x;
  	return y;
  } 
  
  public void setCartesian(int i, double x) {
  	if(i==0) this.x=x;
  	else this.y=x;
  }
    
  public Vector_2 sum(Vector_ v) {
  	return new Vector_2(this.x+v.getCartesian(0),
  						this.y+v.getCartesian(1));  	
  }
  
  public Vector_2 difference(Vector_ v) {
  	return new Vector_2(v.getCartesian(0)-x,
  						v.getCartesian(1)-y);  	
  }
  
  public Vector_2 opposite() {
  	return new Vector_2(-x,-y);  	
  }
  
  public double innerProduct(Vector_ v) {
  	return this.x*v.getCartesian(0)+this.y*v.getCartesian(1);  	
  }

  public Vector_2 divisionByScalar(double s) {
  	return new Vector_2(x/s,y/s);  	
  }
  
  public Vector_2 multiplyByScalar(double s) {
  	return new Vector_2(x*s,y*s);  	
  }
  
  public double squaredLength() {
  	return innerProduct(this);  	
  }
  
}




