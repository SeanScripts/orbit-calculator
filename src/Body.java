
public class Body {
	// Mass (kg)
	double mass;
	// Position (m)
	double x;
	double y;
	double z;
	// Velocity (m/s)
	double vx;
	double vy;
	double vz;
	// Pole direction (unit vector)
	double px;
	double py;
	double pz;
	// Rotation angle
	double r;
	
	public Body(double mass) {
		this.mass = mass;
		// Other things need to be initialized.
		x = y = z = 0;
		vx = vy = vz = 0;
	}
	
	public void update(double t) {
		
	}
	
	// Getters and setters
	
	public double getMass() {
		return mass;
	}

	public void setMass(double mass) {
		this.mass = mass;
	}

	public double getX() {
		return x;
	}

	public void setX(double x) {
		this.x = x;
	}

	public double getY() {
		return y;
	}

	public void setY(double y) {
		this.y = y;
	}

	public double getZ() {
		return z;
	}

	public void setZ(double z) {
		this.z = z;
	}

	public double getVx() {
		return vx;
	}

	public void setVx(double vx) {
		this.vx = vx;
	}

	public double getVy() {
		return vy;
	}

	public void setVy(double vy) {
		this.vy = vy;
	}

	public double getVz() {
		return vz;
	}

	public void setVz(double vz) {
		this.vz = vz;
	}
	
	public double getPx() {
		return px;
	}

	public void setPx(double px) {
		this.px = px;
	}

	public double getPy() {
		return py;
	}

	public void setPy(double py) {
		this.py = py;
	}

	public double getPz() {
		return pz;
	}

	public void setPz(double pz) {
		this.pz = pz;
	}
	
	public double getR() {
		return r;
	}
	
	public void setR(double r) {
		this.r = r;
	}
}
