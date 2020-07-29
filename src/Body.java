
public class Body {
	// Mass (kg)
	double mass;
	// Position (m)
	Point3D<Double> p;
	// Velocity (m/s)
	Point3D<Double> v;
	// Axis direction (unit vector)
	Point3D<Double> a;
	// Rotation angle
	double r;
	
	public Body(double mass) {
		this.mass = mass;
		// Other things need to be initialized.
		p = Point3D.ORIGIN;
		v = Point3D.ORIGIN;
		a = Point3D.Z;
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

	public Point3D<Double> getP() {
		return p;
	}

	public void setP(Point3D<Double> p) {
		this.p = p;
	}

	public Point3D<Double> getV() {
		return v;
	}

	public void setV(Point3D<Double> v) {
		this.v = v;
	}

	public Point3D<Double> getA() {
		return a;
	}

	public void setA(Point3D<Double> a) {
		this.a = a;
	}
	
	public double getR() {
		return r;
	}
	
	public void setR(double r) {
		this.r = r;
	}
}
