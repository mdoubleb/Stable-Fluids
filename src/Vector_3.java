public class Vector_3 extends Vector_2 {
  public double z;

  public Vector_3() {}
  
  public Vector_3(double x,double y,double z) { 
  	this.x=x; 
  	this.y=y; 
  	this.z=z; 
  }
  
  public Vector_3(double[] coordinates) { 
	  	this.x=coordinates[0]; 
	  	this.y=coordinates[1]; 
	  	this.z=coordinates[2];
  }

  public Vector_3(Point_3 a, Point_3 b) { 
  	this.x=b.getX()-a.getX(); 
  	this.y=b.getY()-a.getY(); 
  	this.z=b.getZ()-a.getZ(); 
  }
  
  public double getZ() {return z; }
  
  public void setZ(double z) {this.z=z; }
  
  public boolean equals(Vector_ v) { 
  	return(this.x==v.getCartesian(0) 
  		&& this.y==v.getCartesian(1) 
  		&& this.z==v.getCartesian(2)); 
  }
  
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

  public String toString() {return "["+x+","+y+","+z+"]"; }
  public int dimension() { return 3;}
  
  public double getCartesian(int i) {
  	if(i==0) return x;
  	else if(i==1) return y;
  	else return z;
  } 
  
  public void setCartesian(int i, double x) {
  	if(i==0) this.x=x;
  	else if(i==1) this.y=x;
  	else this.z=x;
  }
    
  public Vector_3 sum(Vector_ v) {
  	return new Vector_3(this.x+v.getCartesian(0),
  						this.y+v.getCartesian(1), 
  						this.z+v.getCartesian(2));  	
  }
  
  public Vector_3 difference(Vector_ v) {
  	return new Vector_3(v.getCartesian(0)-x,
  						v.getCartesian(1)-y, 
  						v.getCartesian(2)-z);  	
  }
  
  public Vector_3 opposite() {
  	return new Vector_3(-x,-y,-z);  	
  }
  
  public double innerProduct(Vector_ v) {
  	return this.x*v.getCartesian(0)
  		  +this.y*v.getCartesian(1)
  		  +this.z*v.getCartesian(2);  	
  }

  public Vector_3 divisionByScalar(double s) {
  	if(s==0.) throw new Error("error: division by zero");
  	return new Vector_3(x/s,y/s,z/s);  	
  }
  
  public Vector_3 multiplyByScalar(double s) {
  	return 
  		new Vector_3(x*s,y*s,z*s);  	
  }
  
  public double squaredLength() {
  	return innerProduct(this);  	
  }
  
  public Vector_3 crossProduct(Vector_ b) {
  	return 
  		new Vector_3(y*b.getCartesian(2)-z*b.getCartesian(1),
  					 z*b.getCartesian(0)-x*b.getCartesian(2),
  					 x*b.getCartesian(1)-y*b.getCartesian(0));
  }

  
}




