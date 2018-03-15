/**
 * A pixel grid for the fluid simulation
 * 
 * @author Luca Castelli Aleardi (INF555, 2012)
 *
 */
public class Grid{
	
	Draw frame; // drawing frame
	Point_2[][] points;  //centers of pixels = "particles"
	int[][][] textures;  //the textures of the pixels
	int[][] transparencies;  //transparency coordinate for substance introduction
	int n; 
	float dx;
	float dy;
	float dt;
	
	public Grid(Draw frame, int n, float dx, float dy, int[][][] textures, int sM) {
		this.textures=textures;
		this.frame=frame;
		this.n=n;
		this.dx = dx;
		this.dy=dy;
		this.points = new Point_2[n][n];
		for(int i = 0; i<n;i++) {
			for(int j = 0; j<n;j++) {
				points[i][j]= new Point_2((j+0.5)*dx, (i+0.5)*dy);
			}
		}
		transparencies=new int[n][n];
		if (sM==1){
			for(int i = 0; i<n;i++) {
				for(int j = 0; j<n;j++) {
					transparencies[i][j]= 255;
				}
			}
		}
	}
	
	public void drawSquare(int[] texture, Point_2 p, int transparency){
		//this.frame.stroke(texture[0],texture[1],texture[2],transparency);
		this.frame.noStroke();
		this.frame.fill(texture[0],texture[1],texture[2], transparency);
		this.frame.rect((float)p.x - dx/2, (float)p.y -dy/2, dx, dy);
	}
	
	public void drawGrid() {
		//gaussian();
		for(int i = 0; i<n;i++) {
			for(int j = 0; j<n;j++) {
				drawSquare(textures[i][j],points[i][j], transparencies[i][j]);
			}
		}
	}
	public void gaussian() {
		int[][][] texturesBis = new int[n][n][3];
		for(int i = 1; i<n-1;i++) {
			for(int j = 1; j<n-1;j++) {
				for(int k = 0;k<2;k++) {
					for(int l = -3;l<4;l++) {
						for(int m = -3;m<4;m++) {
							texturesBis[i][j][k]+=1./9*textures[i+l][j+m][k];
						}
					}
				}
			}
		}
		this.textures=texturesBis;
	}
	
	public void drawPoint(Point_2 p) {
		this.frame.ellipse((float)p.getX(), (float)p.getY(), 5, 5);
	}
	
	public void setTexture(int i, int j, int k, int t){
		textures[i][j][k]=t;
	}
	
	public void setTransparency(int i, int j, int t){
		transparencies[i][j]=t;
	}
}
