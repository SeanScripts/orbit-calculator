import java.awt.Canvas;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;

import javax.swing.JFrame;

public class OrbitDisplay extends Canvas implements MouseListener, KeyListener, MouseMotionListener {
	public static void main(String[] args) {
		JFrame f = new JFrame("Orbit Visualization");
		f.setDefaultCloseOperation(3);
		f.setSize(1000, 1000);
		OrbitDisplay display = new OrbitDisplay();
		f.add(display);
		f.setVisible(true);
		// Basic draw loop, no adjustment to try to match the frame rate.
		boolean running = true;
		while (running) {
			display.step();
			display.repaint();
			try {
				Thread.sleep(16);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	int step = 0;
	ArrayList<Body> bodies;
	ArrayList<Orbit> orbits;
	
	// Controls
	boolean forward = false;
	boolean backward = false;
	boolean left = false;
	boolean right = false; 
	boolean up = false;
	boolean down = false;
	
	// Camera
	Point3D<Double> cp;
	//double cx = 10.0;
	//double cy = 0.0;
	//double cz = 0.0;
	//double cspeed = 1.0;
	Point3D<Double> cs;
	Point3D<Double> cv;
	//double cvx = 0.0;
	//double cvy = 0.0;
	//double cvz = 0.0;
	double friction = 0.5;
	double yaw = 0.0;
	double pitch = 0.0;
	double rr = 0.1;
	
	double scale = 1;
	
	public OrbitDisplay() {
		// Initialize here
		addKeyListener(this);
		addMouseListener(this);
		addMouseMotionListener(this);
		
		cp = new Point3D<Double>(10.0, 0.0, 0.0);
		cs = Point3D.X.scale(0.5);
		cv = Point3D.ORIGIN;
		
		bodies = new ArrayList<Body>();
		orbits = new ArrayList<Orbit>();
		
		Body test = new Body(1);
		test.setP(Point3D.ORIGIN);
		bodies.add(test);
	}
	
	public void step() {
		step++;
		// Time step here
		if (forward) {
			cv = cv.add(cs.rotateXY(pitch).rotateXZ(yaw));
			//cvx += -cv*Math.cos(pitch)*Math.sin(yaw);
			//cvy += cv*Math.sin(pitch);
			//cvz += -cv*Math.cos(pitch)*Math.cos(yaw);
		}
		if (backward) {
			cv = cv.add(cs.rotateXY(pitch).rotateXZ(yaw).negate());
			//cvx += cv*Math.cos(pitch)*Math.sin(yaw);
			//cvy += -cv*Math.sin(pitch);
			//cvz += cv*Math.cos(pitch)*Math.cos(yaw);
		}
		if (left) {
			cv = cv.add(cs.rotateXZ(yaw+Math.PI/2.0));
			//cvx += -cv*Math.cos(yaw);
			//cvz += cv*Math.sin(yaw);
		}
		if (right) {
			cv = cv.add(cs.rotateXZ(yaw-Math.PI/2.0));
			//cvx += cv*Math.cos(yaw);
			//cvz += -cv*Math.sin(yaw);
		}
		if (up) {
			cv = cv.add(cs.rotateXY(Math.PI/2.0));
			//cvy += cv;
		}
		if (down) {
			cv = cv.add(cs.rotateXY(-Math.PI/2.0));
			//cvy += -cv;
		}
		
		cp = cp.add(cv);
		//cx += cvx;
		//cy += cvy;
		//cz += cvz;
		
		cv = cv.scale(friction);
		//cvx *= friction;
		//cvy *= friction;
		//cvz *= friction;
	}
	
	// Graphics
	
	public void update(Graphics g) {
		paint(g);
	}
	
	public void paint(Graphics g) {
		Image mem = createImage(1000, 1000);
		Graphics g2 = (Graphics)(mem.getGraphics());
		
		g2.setColor(Color.BLACK);
		g2.fillRect(0, 0, 1000, 1000);
		
		// Draw graphics here
		
		g2.setColor(Color.WHITE);
		g2.drawString(cp.toString(), 50, 50);
		//g2.drawString("x: "+(int)(cx*100)/100, 50, 50);
		//g2.drawString("y: "+(int)(cy*100)/100, 50, 100);
		//g2.drawString("z: "+(int)(cz*100)/100, 50, 150);
		
		g2.drawString(""+step, 50, 200);
		
		Point3D<Double> refx = Point3D.Z.negate().rotateXY(pitch).rotateXZ(yaw);
		Point3D<Double> refy = Point3D.Y.rotateXY(pitch).rotateXZ(yaw);
		Point3D<Double> refz = Point3D.X.rotateXY(pitch).rotateXZ(yaw);
		
		for (Body body : bodies) {
			//TODO
			int tx = (int)(500 + scale*body.getP().sub(cp).dot(refx));
			int ty = (int)(500 - scale*body.getP().sub(cp).dot(refy));
			double tz = scale*body.getP().sub(cp).dot(refz);
			g2.drawString(""+tz, 100, 100);
			int size = (int)(100/tz);
			g2.fillOval(tx-size, ty-size, 2*size, 2*size);
		}
		//g2.fillOval(500, 500, 10, 10);
		
		g.drawImage(mem, 0, 0, null);
	}
	
	// Interaction
	
	public void keyPressed(KeyEvent e) {
		if (e.getKeyCode() == KeyEvent.VK_W) {
			forward = true;
		}
		if (e.getKeyCode() == KeyEvent.VK_A) {
			left = true;
		}
		if (e.getKeyCode() == KeyEvent.VK_S) {
			backward = true;
		}
		if (e.getKeyCode() == KeyEvent.VK_D) {
			right = true;
		}
		if (e.getKeyCode() == KeyEvent.VK_SPACE) {
			up = true;
		}
		if (e.getKeyCode() == KeyEvent.VK_SHIFT) {
			down = true;
		}
	}
	
	public void keyReleased(KeyEvent e) {
		if (e.getKeyCode() == KeyEvent.VK_W) {
			forward = false;
		}
		if (e.getKeyCode() == KeyEvent.VK_A) {
			left = false;
		}
		if (e.getKeyCode() == KeyEvent.VK_S) {
			backward = false;
		}
		if (e.getKeyCode() == KeyEvent.VK_D) {
			right = false;
		}
		if (e.getKeyCode() == KeyEvent.VK_SPACE) {
			up = false;
		}
		if (e.getKeyCode() == KeyEvent.VK_SHIFT) {
			down = false;
		}
	}
	
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	public void mouseDragged(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	public void mouseMoved(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
}
