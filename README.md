# Stable-Fluids : a Java implementation of a physics-based rendering algorithm to simulate fluid motion

An interactive UI that simulates fluid motion, based on fluid dynamics equations. Implementation of Jos Stam's [paper](http://movement.stanford.edu/courses/cs448-01-spring/papers/stam.pdf). The user can either set liquids in motion, create and move gas particles, with different visual environments and boundary conditions.

## Usage
### Dependencies
The program uses the following libraries : 
-[Processing](https://processing.org)
-[MTJ (Matrix Toolkit Java)](https://github.com/fommil/matrix-toolkits-java)
You have to add the following JARs to the build path of your project : 
- core.jar (Processing)
- sparse-eig-01.1.jar and mtj-0.9.14.jar (MTJ)
For convenience we included the JARs in the [lib](lib/) folder of the repo as well as the source code for the [MTJ library](src/linalg/).

### Running the code
Run the ```Draw.java``` class.

### Interacting with the program
You can interact with the UI by using keys and mouse actions : 
- z : stop/restart
- left click & drag : add forces to set the fluid in motion
- a : reset with same parameters
- g : handle intensity of the mouse force with scrollbar
- b : change boundary conditions (periodic or smooth plastic)
- v : change viscosity 
- i : change mouse mode (tracked or straight strokes)
- m : change the simulation mode : MODE 1 : Fluid motion (liquid-like) or MODE 2: Introduction of particles (gas)
- c : change the visual style (different colors)
When using MODE 2 :
- right click & drag : add sources
- s : change substance

## References
Jos Stam. Stable fluids. In Proceedings of the 26th annual conference on Computer graphics and interactive techniques, pages 121–128. ACM Press/Addison-Wesley Publishing Co., 1999. [Link to paper](http://movement.stanford.edu/courses/cs448-01-spring/papers/stam.pdf)

Authors : Mikhaïl Bessa and Benjamin Barral
