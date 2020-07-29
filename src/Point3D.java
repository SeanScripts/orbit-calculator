import java.util.ArrayList;

/**
 * A point with three numeric dimensions of type N.
 */
public class Point3D<N extends Number> implements Comparable<Point3D<N>>
{
	private N x;
	private N y;
	private N z;
	
	public static final Point3D<Double> ORIGIN = new Point3D<Double>(0.0, 0.0, 0.0);
	public static final Point3D<Double> X = new Point3D<Double>(1.0, 0.0, 0.0);
	public static final Point3D<Double> Y = new Point3D<Double>(0.0, 1.0, 0.0);
	public static final Point3D<Double> Z = new Point3D<Double>(0.0, 0.0, 1.0);
	
	public static final double maxError = 1E-6;
	
	/** Creates a new point at (0, 0, 0) */
	public Point3D() {}
	/**
	 * Creates a new point with the specified coordinates
	 * @param a - the x value
	 * @param b - the y value
	 * @param c - the z value
	 */
	public Point3D(N a, N b, N c)
	{
		x = a; y = b; z = c;
	}
	
	public Point3D<N> clone() {
		return new Point3D<N>(this.x, this.y, this.z);
	}
	
	public static Point3D<Double> randomUnitVector() {
		return new Point3D<Double>(-1 + 2*Math.random(), -1 + 2*Math.random(), -1 + 2*Math.random()).normalize();
	}
	
	/**
	 * Casts this point to an integral 3D point by truncation.
	 * @return the truncated point
	 */
	public Point3D<Integer> truncate() {
		return new Point3D<Integer>(x.intValue(), y.intValue(), z.intValue());
	}
	
	/**
	 * Casts this point to an integral 3D point such that -0.2 would become -1, and 0.2 would become 0, for example.
	 * @return the cast point
	 */
	public Point3D<Integer> castUnique() {
		return new Point3D<Integer>(
			(x.doubleValue() < 0 && x.doubleValue() - x.intValue() != 0) ? x.intValue()-1 : x.intValue(), 
			(y.doubleValue() < 0 && y.doubleValue() - y.intValue() != 0) ? y.intValue()-1 : y.intValue(), 
			(z.doubleValue() < 0 && z.doubleValue() - z.intValue() != 0) ? z.intValue()-1 : z.intValue()
		);
	}
	
	/**
	 * Considers this point a 3D vector and normalizes it such that it has magnitude 1 in whatever direction it was facing.
	 * @return the normalized vector
	 */
	public Point3D<Double> normalize() {
		double d = dist();
		if (d == 0)
			return new Point3D<Double>(0.0, 0.0, 0.0);
		return new Point3D<Double>(x.doubleValue()/d, y.doubleValue()/d, z.doubleValue()/d);
	}
	
	/**
	 * Negates this point's coordinates.
	 * @return the negated point
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public Point3D<N> negate() {
		if (x instanceof Float || x instanceof Double) {
			return new Point3D(
				-this.x.doubleValue(), 
				-this.y.doubleValue(), 
				-this.z.doubleValue()
			);
		}
		return new Point3D(
			-this.x.intValue(), 
			-this.y.intValue(), 
			-this.z.intValue()
		);
	}
	
	/**
	 * Scale this point by a factor
	 * @param amt - the factor to scale it by
	 * @return the scaled point
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public Point3D<N> scale(N amt) {
		if (x instanceof Float || x instanceof Double || amt instanceof Float || amt instanceof Double) {
			return new Point3D(
				this.x.doubleValue() * amt.doubleValue(), 
				this.y.doubleValue() * amt.doubleValue(), 
				this.z.doubleValue() * amt.doubleValue()
			);
		}
		return new Point3D(
			this.x.intValue() * amt.intValue(), 
			this.y.intValue() * amt.intValue(), 
			this.z.intValue() * amt.intValue()
		);
	}
	
	/**
	 * Scale this point by a factor
	 * @param amt - the factor to scale it by
	 * @return the scaled point
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public Point3D<N> scale(String dir, N amt) {
		if (x instanceof Float || x instanceof Double || amt instanceof Float || amt instanceof Double) {
			return new Point3D(
				(dir.equals("X")) ? this.x.doubleValue() * amt.doubleValue() : this.x.doubleValue(), 
				(dir.equals("Y")) ? this.y.doubleValue() * amt.doubleValue() : this.y.doubleValue(), 
				(dir.equals("Z")) ? this.z.doubleValue() * amt.doubleValue() : this.z.doubleValue()
			);
		}
		return new Point3D(
			(dir.equals("X")) ? this.x.intValue() * amt.intValue() : this.x.intValue(), 
			(dir.equals("Y")) ? this.y.intValue() * amt.intValue() : this.y.intValue(), 
			(dir.equals("Z")) ? this.z.intValue() * amt.intValue() : this.z.intValue()
		);
	}
	
	/**
	 * Scale this point by a factor
	 * @param amt - the factor to scale it by
	 * @return the scaled point
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public Point3D<N> scale(String dir, N amtPlus, N amtMinus) {
		if (x instanceof Float || x instanceof Double || amtPlus instanceof Float || amtPlus instanceof Double || amtMinus instanceof Float || amtMinus instanceof Double) {
			return new Point3D(
				(dir.equals("X")) ? ( (this.x.doubleValue() >= 0) ? this.x.doubleValue() * amtPlus.doubleValue() : this.x.doubleValue() * amtMinus.doubleValue() ) : this.x.doubleValue(), 
				(dir.equals("Y")) ? ( (this.y.doubleValue() >= 0) ? this.y.doubleValue() * amtPlus.doubleValue() : this.y.doubleValue() * amtMinus.doubleValue() ) : this.y.doubleValue(), 
				(dir.equals("Z")) ? ( (this.z.doubleValue() >= 0) ? this.z.doubleValue() * amtPlus.doubleValue() : this.z.doubleValue() * amtMinus.doubleValue() ) : this.z.doubleValue()
			);
		}
		return new Point3D(
			(dir.equals("X")) ? ( (this.x.intValue() >= 0) ? this.x.intValue() * amtPlus.intValue() : this.x.intValue() * amtMinus.intValue() ) : this.x.intValue(), 
			(dir.equals("Y")) ? ( (this.y.intValue() >= 0) ? this.y.intValue() * amtPlus.intValue() : this.y.intValue() * amtMinus.intValue() ) : this.y.intValue(), 
			(dir.equals("Z")) ? ( (this.z.intValue() >= 0) ? this.z.intValue() * amtPlus.intValue() : this.z.intValue() * amtMinus.intValue() ) : this.z.intValue()
		);
	}
	
	/**
	 * Finds the distance in three dimensions to another point.
	 * @param p - the other point
	 * @return the distance
	 */
	public double dist(Point3D<N> p) {
		if (p != null) {
			return Math.sqrt( 
				  Math.pow(p.x.doubleValue() - x.doubleValue(), 2) 
				+ Math.pow(p.y.doubleValue() - y.doubleValue(), 2) 
				+ Math.pow(p.z.doubleValue() - z.doubleValue(), 2)
			);
		}
		return -1.0f;
	}
	
	/**
	 * Finds the distance of this point from the origin.
	 * @return the distance
	 */
	public double dist() {
		return Math.sqrt( x.doubleValue()*x.doubleValue() + y.doubleValue()*y.doubleValue() + z.doubleValue()*z.doubleValue() );
	}
	
	/**
	 * Sets the x, y, and  values to be equivalent to another point
	 * @param p - the other point
	 */
	public void set(Point3D<N> p) {
		if (p != null) {
			x = p.x; y = p.y; z = p.z;
		}
	}
	
	/**
	 * Adds the coordinates of this point to another point.
	 * @param p - the other point
	 * @return the resulting coordinates as a new point
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public Point3D<N> add(Point3D<N> p) {
		if (p.x instanceof Float && x instanceof Float) {
			return new Point3D(
				this.x.floatValue() + p.x.floatValue(), 
				this.y.floatValue() + p.y.floatValue(), 
				this.z.floatValue() + p.z.floatValue()
			);
		}
		else if (p.x instanceof Float || p.x instanceof Double || x instanceof Float || x instanceof Double) {
			return new Point3D(
				this.x.doubleValue() + p.x.doubleValue(), 
				this.y.doubleValue() + p.y.doubleValue(), 
				this.z.doubleValue() + p.z.doubleValue()
			);
		}
		return new Point3D(
			this.x.intValue() + p.x.intValue(), 
			this.y.intValue() + p.y.intValue(), 
			this.z.intValue() + p.z.intValue()
		);
	}
	
	/**
	 * Subtracts the coordinates of another point from this point.
	 * @param p - the other point
	 * @return the resulting coordinates as a new point
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public Point3D<N> sub(Point3D<N> p) {
		if (p.x instanceof Float || p.x instanceof Double || x instanceof Float || x instanceof Double) {
			return new Point3D(
				this.x.doubleValue() - p.x.doubleValue(), 
				this.y.doubleValue() - p.y.doubleValue(), 
				this.z.doubleValue() - p.z.doubleValue()
			);
		}
		return new Point3D(
			this.x.intValue() - p.x.intValue(), 
			this.y.intValue() - p.y.intValue(), 
			this.z.intValue() - p.z.intValue()
		);
	}
	
	/**
	 * Scales the coordinates of this point by another point
	 * @param p - the other point
	 * @return the resulting coordinates as a new point
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public Point3D<N> scaleBy(Point3D<N> p) {
		if (p.x instanceof Float && x instanceof Float) {
			return new Point3D(
				this.x.floatValue() * p.x.floatValue(), 
				this.y.floatValue() * p.y.floatValue(), 
				this.z.floatValue() * p.z.floatValue()
			);
		}
		else if (p.x instanceof Float || p.x instanceof Double || x instanceof Float || x instanceof Double) {
			return new Point3D(
				this.x.doubleValue() * p.x.doubleValue(), 
				this.y.doubleValue() * p.y.doubleValue(), 
				this.z.doubleValue() * p.z.doubleValue()
			);
		}
		return new Point3D(
			this.x.intValue() * p.x.intValue(), 
			this.y.intValue() * p.y.intValue(), 
			this.z.intValue() * p.z.intValue()
		);
	}
	
	/**
	 * "( x, y, z )"
	 */
	public String toString() {	
		if (x instanceof Float || x instanceof Double) {
			double tx = round(3, x.doubleValue());
			double ty = round(3, y.doubleValue());
			double tz = round(3, z.doubleValue());
			
			return "( "+tx+", "+ty+", "+tz+" )";
		}
		
		return "( "+x+", "+y+", "+z+" )";
	}
	
	/**
	 * "( x, y, z )"
	 */
	public String toString(int dec) {	
		if (x instanceof Float || x instanceof Double) {
			double tx = round(dec, x.doubleValue());
			double ty = round(dec, y.doubleValue());
			double tz = round(dec, z.doubleValue());
			
			return tx+" "+ty+" "+tz;
		}
		
		return x+" "+y+" "+z;
	}
	
	/**
	 * Compares this point to another point. If it returns zero, they are the same. 
	 */
	public int compareTo(Point3D<N> o) {
		if ( z.doubleValue() == o.z.doubleValue() ) {
			if ( y.doubleValue() == o.y.doubleValue() ) {
				if ( x.doubleValue() == o.x.doubleValue() ) {
					return 0;
				}
				return Double.compare(x.doubleValue(), o.x.doubleValue());
			}
			return Double.compare(y.doubleValue(), o.y.doubleValue());
		}
		return Double.compare(z.doubleValue(), o.z.doubleValue());
	}
	
	public N getX() {
		return x;
	}
	
	public N getY() {
		return y;
	}
	
	public N getZ() {
		return z;
	}
	
	public void setX(N value) {
		x = value;
	}
	
	public void setY(N value) {
		y = value;
	}
	
	public void setZ(N value) {
		z = value;
	}
	
	public static Point3D<Double> average(ArrayList<Point3D<Double>> list) {
		Point3D<Double> res = new Point3D<Double>(0.0, 0.0, 0.0);
		if (list.size() == 0)
			return res;
		for (Point3D<Double> point : list) {
			res = res.add(point);
		}
		return res.scale(1.0/(double)list.size());
	}
	
	public static Point3D<Double> midpoint(ArrayList<Point3D<Double>> list) {
		double minx = 10000000;
		double maxx = -10000000;
		double miny = 10000000;
		double maxy = -10000000;
		double minz = 10000000;
		double maxz = -10000000;
		for (Point3D<Double> pt : list) {
			if (pt.getX() < minx)
				minx = pt.getX();
			if (pt.getX() > maxx)
				maxx = pt.getX();
			if (pt.getY() < miny)
				miny = pt.getY();
			if (pt.getY() > maxy)
				maxy = pt.getY();
			if (pt.getZ() < minz)
				minz = pt.getZ();
			if (pt.getZ() > maxz)
				maxz = pt.getZ();
		}
		return new Point3D<Double>((minx+maxx)/2.0, (miny+maxy)/2.0, (minz+maxz)/2.0);
	}
	
	/**
	 * Finds the dot product of the two vectors
	 * @param o - the other vector
	 * @return this dot o
	 */
	public double dot(Point3D<Double> o) {
		return this.x.doubleValue()*o.getX().doubleValue() + this.y.doubleValue()*o.getY().doubleValue() + this.z.doubleValue()*o.getZ().doubleValue();
	}
	
	/**
	 * Finds the cross product of the two vectors
	 * @param o - the other vector
	 * @return this cross o
	 */
	public Point3D<Double> cross(Point3D<Double> o) {
		return new Point3D<Double>(
			this.y.doubleValue()*o.getZ().doubleValue() - this.z.doubleValue()*o.getY().doubleValue(),
			this.z.doubleValue()*o.getX().doubleValue() - this.x.doubleValue()*o.getZ().doubleValue(),
			this.x.doubleValue()*o.getY().doubleValue() - this.y.doubleValue()*o.getX().doubleValue()
		);
	}
	
	public static double[] generateRandomRotation(double maxRotationRate) {
		double[] res = {-maxRotationRate + 2*maxRotationRate*Math.random(), -maxRotationRate + 2*maxRotationRate*Math.random(), -maxRotationRate + 2*maxRotationRate*Math.random()};
		return res;
	}
	
	public static Point3D<Double> applyRandomRotation(double time, double[] rotationParameters, Point3D<Double> point) {
		double[] rotationParametersAtTime = {rotationParameters[0]*time, rotationParameters[1]*time, rotationParameters[2]*time, 
			rotationParameters[3]*time, rotationParameters[4]*time, rotationParameters[5]*time};
		return point.rotate(rotationParametersAtTime[0], rotationParametersAtTime[1], rotationParametersAtTime[2]);
	}
	
	/**
	 * Applies a rotation to a point about another point
	 * @param about - the point about which to rotate
	 * @param xy - the angle in the xy plane
	 * @param xz - the angle in the xz plane
	 * @param xw - the angle in the xw plane
	 * @param yz - the angle in the yz plane
	 * @param yw - the angle in the yw plane
	 * @param zw - the angle in the zw plane
	 * @return the rotated point
	 */
	public Point3D<Double> rotate(double xy, double xz, double yz) {
		//order to rotate is xy, xz, xw, yz, yw, zx
		return this.rotateXY(xy).rotateXZ(xz).rotateYZ(yz);
	}
	
	public Point3D<Double> rotateXY(double angle) {
		return new Point3D<Double>(
			Math.cos(angle)*this.x.doubleValue() - Math.sin(angle)*this.y.doubleValue(),
			Math.sin(angle)*this.x.doubleValue() + Math.cos(angle)*this.y.doubleValue(),
			this.z.doubleValue()
		);
	}
	
	public Point3D<Double> rotateXZ(double angle) {
		return new Point3D<Double>(
			Math.cos(angle)*this.x.doubleValue() - Math.sin(angle)*this.z.doubleValue(), 
			this.y.doubleValue(), 
			Math.sin(angle)*this.x.doubleValue() + Math.cos(angle)*this.z.doubleValue()
		);
	}
	
	public Point3D<Double> rotateYZ(double angle) {
		return new Point3D<Double>(
			this.x.doubleValue(), 
			Math.cos(angle)*this.y.doubleValue() - Math.sin(angle)*this.z.doubleValue(),
			Math.sin(angle)*this.y.doubleValue() + Math.cos(angle)*this.z.doubleValue()
		);
	}
	
	public Point3D<Double> reflect(String axis) {
		if (axis.equals("X")) {
			return new Point3D<Double>(-this.x.doubleValue(), this.y.doubleValue(), this.z.doubleValue());
		}
		else if (axis.equals("Y")) {
			return new Point3D<Double>(this.x.doubleValue(), -this.y.doubleValue(), this.z.doubleValue());
		}
		else if (axis.equals("Z")) {
			return new Point3D<Double>(this.x.doubleValue(), this.y.doubleValue(), -this.z.doubleValue());
		}
		else {
			return new Point3D<Double>(this.x.doubleValue(), this.y.doubleValue(), this.z.doubleValue());
		}
	}
	
	public static Point3D<Double> convertIntegerToDouble(Point3D<Integer> i){
		return new Point3D<Double>(i.getX().doubleValue(), i.getY().doubleValue(), i.getZ().doubleValue());
	}
	
	public static double round(int decPlaces, double num) {
		long temp = Math.round(num*(double)Math.pow(10,decPlaces));
		return temp/(double)Math.pow(10,decPlaces);
	}
	
	public static Point3D<Double> angleIndexNormal(double angle, double index, int gear) {
		double inda = index/(double)gear;
		Point3D<Double> res = new Point3D<Double>(0.0, 0.0, 1.0);
		return res.rotateXZ(-Math.PI*angle/180.0).rotateXY(2*Math.PI*inda);
	}
	
	public static Point3D<Double> angleIndexU(double angle, double index, int gear) {
		double inda = index/(double)gear;
		Point3D<Double> res = new Point3D<Double>(1.0, 0.0, 0.0);
		return res.rotateXZ(-Math.PI*angle/180.0).rotateXY(2*Math.PI*inda);
	}
	
	public static Point3D<Double> indexV(double index, int gear) {
		double inda = index/(double)gear;
		Point3D<Double> res = new Point3D<Double>(0.0, 1.0, 0.0);
		return res.rotateXY(2*Math.PI*inda);
	}
	
	//radians
	public static Point3D<Double> getU(double angle, double index) {
		Point3D<Double> res = new Point3D<Double>(1.0, 0.0, 0.0);
		return res.rotateXZ(-angle).rotateXY(index);
	}
	
	//radians
	public static Point3D<Double> getV(double index) {
		Point3D<Double> res = new Point3D<Double>(0.0, 1.0, 0.0);
		return res.rotateXY(index);
	}
	
	//radians
	public static Point3D<Double> getNormal(double angle, double index) {
		Point3D<Double> res = new Point3D<Double>(0.0, 0.0, 1.0);
		return res.rotateXZ(-angle).rotateXY(index);
	}
	
	public static Point3D<Double> getDirection(double angle, double index, double yaw, double pitch) {
		return new Point3D<Double>(Math.cos(yaw)*Math.cos(pitch), Math.sin(yaw)*Math.cos(pitch), Math.sin(pitch)).rotateXZ(Math.PI/2-angle).rotateXY(index).negate();
	}
	
	public static double pointLineDist(Point3D<Double> p, Point3D<Double> a, Point3D<Double> n) {
		Point3D<Double> tmp = a.sub(p);
		return (tmp.sub(n.scale(tmp.dot(n)))).dist();
	}

	public boolean equals(Point3D<N> p) {
		return (this.dist(p) < maxError);
	}
	
	public Point3D<Double> scaleZ(double amtPlus, double amtMinus, double minz, double maxz) {
		double tz = z.doubleValue();
		//Requirements:
		// if z >= minz, tz = z
		// otherwise scale the difference of z and minz by amtMinus and set that to tz
		// so that tz = minz - (minz - z)*amtMinus
		// if z <= maxz, tz = z
		// otherwise scale the difference of z and maxz by amtPlus and set that to tz
		// so that tz = maxz + (z - maxz)*amtPlus
		if (z.doubleValue() < minz)
			tz = minz - (minz - z.doubleValue())*amtMinus;
		else if (z.doubleValue() > maxz)
			tz = maxz + (z.doubleValue() - maxz)*amtPlus;
		return new Point3D<Double>(x.doubleValue(), y.doubleValue(), tz);
	}
	
	public Point3D<Double> translateZ(double amt) {
		double tz = 0;
		if (z.doubleValue() > 0)
			tz = z.doubleValue() + amt;
		else if (z.doubleValue() < 0)
			tz = z.doubleValue() - amt;
		return new Point3D<Double>(x.doubleValue(), y.doubleValue(), tz);
	}
	
	//Helper methods that avoid creating new points until the end.
	public static Point3D<Double> getReflectionDir(Point3D<Double> dir, Point3D<Double> normal, double c) {
		//dir.add(intersectionTriangle.n.scale(2*c)).normalize();
		double px = dir.getX() + normal.getX()*2*c;
		double py = dir.getY() + normal.getY()*2*c;
		double pz = dir.getZ() + normal.getZ()*2*c;
		double norm = Math.sqrt(px*px + py*py + pz*pz);
		px /= norm;
		py /= norm;
		pz /= norm;
		return new Point3D<Double>(px, py, pz);
	}
	
	public static Point3D<Double> getRefractionDir(Point3D<Double> dir, Point3D<Double> normal, double r, double refractionFactor) {
		//(dir.scale(r).add(intersectionTriangle.n.scale(refractionFactor))).normalize();
		double px = dir.getX()*r + normal.getX()*refractionFactor;
		double py = dir.getY()*r + normal.getY()*refractionFactor;
		double pz = dir.getZ()*r + normal.getZ()*refractionFactor;
		double norm = Math.sqrt(px*px + py*py + pz*pz);
		px /= norm;
		py /= norm;
		pz /= norm;
		return new Point3D<Double>(px, py, pz);
	}
	
	public Point3D<Double> rotate(String dir, double amt) {
		if (dir.equals("X")) {
			return this.rotateYZ(amt);
		}
		else if (dir.equals("Y")) {
			return this.rotateXZ(amt);
		}
		else if (dir.equals("Z")) {
			return this.rotateXY(amt);
		}
		else {
			return new Point3D<Double>(this.getX().doubleValue(), this.getY().doubleValue(), this.getZ().doubleValue());
		}
	}
}