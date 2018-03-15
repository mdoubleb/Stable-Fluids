import processing.core.*;

import java.util.*;

/**
 * 
 * @author benjaminbarral and mikhailbessa
 * A simple interface to handle interaction with the fluid
 */

public class Draw extends PApplet {

	///Display
	Grid grid;
	final int frameRate = 16;  // can be set larger depending on the processor
	final int n = 110;   //number of pixels = n*n
	final int N=n * n;   //N=n*n
	final int Lx = 600, Ly = 600;   //width and length of the window
	
	///Viscosity and Diffusivity settings
	//Viscosity
	final float conv=(float) ((1440f/0.3)*(1440f/0.3));   //conversion rate for distance unit : actually depends on your screen
	final float visc_water = 8.84f * pow(10, -7) *conv;
	final float visc_air = 1.56f * pow(10, -5)* conv;
	final float visc_interm=120f * pow(10,-6)*conv;
	final float visc_glycerin=1182f * pow(10,-6)* conv;;//1182 * pow(10,-6)* conv;
	final float visc_mercury=0.1f * pow(10,-6)* conv;
	float visc;
	int viscosityInteger=0;  //parameter for setting the viscosity
	//Diffusivity
	final float kappa_CO2_air=16f*pow(10,-6)*conv;
	final float kappa_CO2_water=0.016f*pow(10,-6)*conv;
	float kappa_virtual=20f;
	float kappa=kappa_CO2_water;
	int kappaInteger=0;

	///Solver and simulation
	AbstractSolver solver;  //the simulation solver
	final double precision = Math.pow(10, -5);   //precision for the linear solver
	final float dx=((float) (Lx)) / n, dy=((float) (Ly)) / n;;   //space discretization parameter
	final float dt = 1f / frameRate;;   //time discretization parameter
	boolean isSimulating;
	int stepCounter;
	
	int[][][] initialTextures, textures;  //rendering textures
	int initialColors = 0;   //parameter for setting the color environment
	
	int simulationMode = 1;   //fluid motion VS substance introduction
	int BC = 1;   //boundary conditions : "smooth plastic" and "periodic"
	
	// Substance introduction settings
	double[] source;
	double[] scalar;
	final int cloud_color=255;  //simulating a cloud in the air
	final int smoke_color=150;   //simulating smoke in the air
	int scalarColor;    //color of the substance 
	int scalarColorInteger=0;   //parameter for setting the color of the substance
	
	/// Force settings
	double[] FXbis, FYbis;   //Force arrays
	final float g = 9.81f;   //gravity; used for scaling the force
	boolean isClicking, endClick, afterClick;   //force clicking booleans
	Point_2 pFstart, pFend;   //mouse positions for force
	boolean infinitesimal = true;   //force computed throughout the mouse path 
	//Scrollbar Force handling
	boolean useScrollbar=false;
	HScrollbar hs1;  

	public void setupMethod() {
		isSimulating = true;
		frameRate(frameRate);
		if (useScrollbar){
			size(Lx, Ly+16); // size of the window
			hs1 = new HScrollbar(0, Ly+8, Lx, 16, 16);
		}
		else size(Lx, Ly);
		initialTextures = new int[n][n][3];
		setInitialTextures();
		textures = initialTextures.clone();
		setViscosity();
		setKappa();
		setScalarColor();
		FXbis = new double[N];
		FYbis = new double[N];
		scalar = new double[N];
		source = new double[N];

		this.grid = new Grid(this, n, dx, dy, textures,simulationMode);

		double[] initialVelocityX = new double[n * n];
		double[] initialVelocityY = new double[n * n]; 
		if (simulationMode == 1)
			this.solver = new Solver(initialVelocityX, initialVelocityY, n, dx, dy, dt, visc, this.grid, precision, BC);
		else
			this.solver = new ScalarSolver(initialVelocityX, initialVelocityY, n, dx, dy, dt, visc, kappa, this.grid,
					precision, BC, scalarColor);
		
		endClick = false;
		afterClick = false;
		
		stepCounter = 0;

		printInformation();
	}

	public void setup() {
		setupMethod();
	}

	public void draw() {
		//background(119, 181, 254); 
		if (simulationMode==0) background(119, 181, 254);    //sky Color for simulationMode=0
		else background(0);
		grid.drawGrid();
		if (isSimulating) {
			solver.setSource(source);
			getForce();
			solver.setForce(FXbis, FYbis);
			solver.routine();
			grid.drawGrid();
			stepCounter++;
			if (useScrollbar){
				hs1.update();
				hs1.display();
			}
		}
	}
	
	void setInitialTextures() {   //3 possible color fluid environments
		initialColors = initialColors % 3;
		switch (initialColors) {
		case (0):
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					initialTextures[i][j][0] = (int) (5 * i * (50. / (n)))%255;// +2*j);
					initialTextures[i][j][1] = (int) 130;// +3*j);
					initialTextures[i][j][2] = (int) (200);
				}
			}
			break;
		case (1):
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					initialTextures[i][j][0] = (int) (21 + (1. + 0.97*i) / n * (234.))%255;// +2*j);
					initialTextures[i][j][1] = (int) (56 + (1. + 0.97*i) / n * 200.)%255;// +3*j);
					initialTextures[i][j][2] = (int) (139 + (1. + 0.97*i) / n * 116.)%255;
				}
			}
			break;
		case (2):
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					initialTextures[i][j][0] = (int) (100+ 2*i)%255;
					initialTextures[i][j][1] = (int) (100+2*j)%255;
					initialTextures[i][j][2] = (int) (i+j)%255;
				}
			}
			break;
		}
	}
	
	void setViscosity() {  //3 possible viscosities
		viscosityInteger = viscosityInteger % 4;
		switch (viscosityInteger) {
		case (0):
			visc=visc_water;
			break;
		case (1):
			visc=visc_air;
			break;
		case (2):
			visc=visc_interm;
			break;
		case (3):
			visc=visc_glycerin;
			break;
		}
		printInformation();
	}
	
	void setKappa(){
		kappaInteger=kappaInteger%3;
		switch(kappaInteger){
		case (0):
			kappa=kappa_CO2_air;
			break;
		case (1):
			kappa=kappa_CO2_water;
			break;
		case (2):
			kappa=kappa_virtual;
			break;
		}
			
	}
	
	void setScalarColor(){   //2 possible substance colors
		scalarColorInteger=scalarColorInteger%2;
		switch(scalarColorInteger){
		case(0):
			scalarColor=cloud_color;break;
		case(1):
			scalarColor=smoke_color;break;
		}
	}
	
	public void reSetup(boolean changeBoundaryConditions, boolean changeInfinitesimalMode,
			boolean changeSimulationMode) {  //
		if (changeBoundaryConditions) {
			BC = 1 - BC;
		}
		if (changeSimulationMode) {
			simulationMode = 1 - simulationMode;
		}
		if (changeInfinitesimalMode) {
			infinitesimal = !infinitesimal;
		}

		setupMethod();
	}

	public void keyPressed() {
		switch (key) {
		case ('z'):
		case ('Z'):
			isSimulating = !isSimulating;
			break; // stop the simulation
		case ('a'):
		case ('A'):
			reSetup(false, false, false);
			break; // reset with previous parameters
		case ('b'):
		case ('B'):
			reSetup(true, false, false);
			break; // reset with other boundary conditions
		case ('i'):
		case ('I'):
			reSetup(false, true, false);
			break; // reset with other type of mouse force (infinitesimal VS total path)
		case ('m'):
		case ('M'):
			reSetup(false, false, true);
			break; // reset with other simulation mode : fluid motion VS susbtance introduction
		case ('c'):
		case ('C'):
			initialColors++;
			reSetup(false, false, false);
			break; // reset with other fluid environment colors
		case ('v'):
		case ('V'):
			viscosityInteger++;
			reSetup(false, false, false);
			break; // reset with other viscosities
		case ('s'):
		case ('S'):
			scalarColorInteger++;
			reSetup(false, false, false);
			break; // reset with other substance colors
		case ('q'):
		case ('Q'):
			useScrollbar=!useScrollbar;
			reSetup(false, false, false);
			break; // reset with using scrollbar or not
		case ('k'):
		case ('K'):
			kappaInteger++;
			reSetup(false, false, false);
			break; // reset with other viscosities
		}
	}
	
	private void printInformation() {
		System.out.println("[Commands]");
		System.out.println("[z : stop/start], [a : reset with current parameters], [c: change initial colors], [v: change viscosity], [i : change mouse mode]");
		System.out.println("[q: handle intensity of the mouse force with scrollbar], [m : change simulation mode], [b : change boundary conditions], [s: change substance introduced]");
		System.out.println();
		System.out.println("[Current parameters]");
		String parameters="";
		if (simulationMode == 0)
			parameters+="Simulation mode: Introduction of substance;  ";
		else
			parameters+="Simulation mode: Fluid motion;  ";
		if (BC == 0)
			parameters+="Boundary conditions: periodic;  ";
		else
			parameters+="Boundary conditions: smooth plastic;  ";
		if (infinitesimal)
			parameters+="Mouse mode: tracked;  ";
		else
			parameters+="Mouse mode: stroke;  ";
		System.out.println(parameters);
		String viscosityString="Viscosity: ";
		switch (viscosityInteger) {
		case (0):
			viscosityString+="viscosity of water;  ";
			break;
		case (1):
			viscosityString+="viscosity of air;  ";
			break;
		case (2):
			viscosityString+="intermediate viscosity ";
			break;
		case (3):
			viscosityString+="viscosity of glycerin";
			break;
		}
		if (simulationMode==0){
			viscosityString+=";  Diffusion of : ";
			switch(scalarColorInteger){
			case(0):
				viscosityString+="cloud in ";break;
			case(1):
				viscosityString+="smoke in ";break;
			}
			switch(kappaInteger){
			case(0):
				viscosityString+="air. ";break;
			case(1):
				viscosityString+="water. ";break;
			case(2):
				viscosityString+="intermediate environment. ";break;
			}
		}
		System.out.println(viscosityString);
		System.out.println("");
	}


	public void mousePressed() { // we start the click
		Point_2 p = new Point_2(mouseX, mouseY);
		if (mouseButton == LEFT) pFstart = p; 
		else if (mouseButton == RIGHT && simulationMode == 0) {
			pFstart = p;
			setSource(p);
		}
	}

	public void mouseDragged() {  //clicking
		Point_2 p = new Point_2(mouseX, mouseY);
		if (mouseButton == LEFT && infinitesimal) {
			pFend = p;
			isClicking = true;
		} 
		else if (mouseButton == RIGHT && simulationMode == 0 ) {
			pFstart = p;
			setSource(p);
		}
	}

	public void mouseReleased() {  //end of click
		if (mouseButton == LEFT) {
			Point_2 p = new Point_2(mouseX, mouseY);
			if (!infinitesimal) {
				pFend = p;
				endClick = true;
			} 
			else {
				pFend = p;
				isClicking = false;
			}
		}
		else if (mouseButton == RIGHT && simulationMode == 0 ) {
			setSource(null);
		}
	}

	void getForce() {
		if (!infinitesimal) {
			if (endClick) {
				getForceClicked(false);
				endClick = false;
				afterClick = true;
			} else if (afterClick) { // removes the force after having set it once
				afterClick = false;
				getForceClicked(true);
			}
		} else {
			if (isClicking) {
				getForceClicked(false);
				afterClick = true;
				pFstart = pFend;
			} else if (afterClick) { // removes the force after having set it once
				afterClick = false;
				getForceClicked(true);
			}
		}
	}
	
	private void setSource(Point_2 p) {   
		//introduction of substance with a certain local density, within 
		//a certain area around point clicked
		if (p != null) {
			int i0 = (int) (p.y / dy - 0.5);
			int j0 = (int) (p.x / dx - 0.5);
			if (i0 < 5 || i0 >= n - 5 || j0 < 5 || j0 >= n - 5)
				return;
			for (int i = -5; i < 5; i++) {
				for (int j = -5; j < 5; j++) {
					if(p.distanceFrom(grid.points[i0+i][j0+j])<4*dx) source[(i0 + i) * n + j + j0] = 150;
				}
			}
		} else {
			for (int i = 0; i < N; i++) {
				source[i] = 0;
			}
		}

	}

	double squaredDistancePointVector(Vector_2 u, Point_2 originPoint, Point_2 point) {
		// distance from point to line passing by originPoint directed by u
		Vector_2 v = (Vector_2) point.minus(originPoint);
		double s = v.innerProduct(u) / u.squaredLength();
		return (s >= 0 && s <= 1) ? ((u.multiplyByScalar(s)).difference(v)).squaredLength() : 2 * Lx * Lx; // we give a
	}
	
	void getForceClicked(boolean reset) { //a "weaping" force due to mouse clicking moves
		if (reset) {   //after click : reset force to zero
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					FXbis[i * n + j] = 0;
					FYbis[i * n + j] = 0;
				}
			}
		} else { // we apply the "weaping" force within a certain width, and for points 
					// located between the two mouse positions
			Vector_2 strengthVector = (Vector_2) pFend.minus(pFstart);
			Point_2 p;
			Vector_2 strength;
			double strengthInitialIntensity;
			if(useScrollbar) strengthInitialIntensity=hs1.getForce();
			else strengthInitialIntensity= (infinitesimal) ? 15.5 * g : 1.55 * g;
			double strengthIntensity;
			double widthWeaped =(simulationMode==1)? ((double) (Lx) / 30) * ((double) (Ly) / 30) : ((double) (Lx) / 20) * ((double) (Ly) / 20);
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					p = grid.points[i][j];
					strengthIntensity = (squaredDistancePointVector(strengthVector, pFstart, p) < widthWeaped)
							? strengthInitialIntensity
							: 0;
					strength = strengthVector.multiplyByScalar(-strengthIntensity);
					FXbis[i * n + j] = strength.getX();
					FYbis[i * n + j] = strength.getY();
				}
			}
		}
	}
	
	///Class for the scrollbar handling the intensity of the mouse force
	class HScrollbar {
		  int swidth, sheight;    // width and height of bar
		  float xpos, ypos;       // x and y position of bar
		  float spos, newspos;    // x position of slider
		  float sposMin, sposMax; // max and min values of slider
		  int loose;              // how loose/heavy
		  boolean over;           // is the mouse over the slider?
		  boolean locked;
		  float ratio;

		  HScrollbar (float xp, float yp, int sw, int sh, int l) {
		    swidth = sw;
		    sheight = sh;
		    int widthtoheight = sw - sh;
		    ratio = (float)sw / (float)widthtoheight;
		    xpos = xp;
		    ypos = yp-sheight/2;
		    spos = xpos + swidth/2 - sheight/2;
		    newspos = spos;
		    sposMin = xpos;
		    sposMax = xpos + swidth - sheight;
		    loose = l;
		  }

		  void update() {
		    if (overEvent()) {
		      over = true;
		    } else {
		      over = false;
		    }
		    if (mousePressed && over) {
		      locked = true;
		    }
		    if (!mousePressed) {
		      locked = false;
		    }
		    if (locked) {
		      newspos = constrain(mouseX-sheight/2, sposMin, sposMax);
		    }
		    if (abs(newspos - spos) > 1) {
		      spos = spos + (newspos-spos)/loose;
		    }
		  }

		  float constrain(float val, float minv, float maxv) {
		    return min(max(val, minv), maxv);
		  }

		  boolean overEvent() {
		    if (mouseX > xpos && mouseX < xpos+swidth &&
		       mouseY > ypos && mouseY < ypos+sheight) {
		      return true;
		    } else {
		      return false;
		    }
		  }

		  void display() {
		    noStroke();
		    fill(204);
		    rect(xpos, ypos, swidth, sheight);
		    if (over || locked) {
		      fill(0, 0, 0);
		    } else {
		      fill(102, 102, 102);
		    }
		    rect(spos, ypos, sheight, sheight);
		  }
		  
		  float getForce(){
			  float forceMin=0.1f*g;
			  float forceMax=15*g;
			  float res=forceMin+(forceMax-forceMin)*spos/swidth;
			  return res;
		  }
		}

}
