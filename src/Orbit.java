public class Orbit {
	static double G = 6.6743015E-11; //Gravitational parameter, in m^3 kg^-1 s^-2
	
	Body primary;
	double at; // axial tilt
	double tp; // tilt phase
	double rr; // rotation rate
	double rp; // rotation phase
	// precession/nutation is not included (yet, anyway).
	
	double a; // semi-major axis
	double e; // eccentricity
	//double m; // mass of smaller object
	//double M; // mass of larger object
	double i; // inclination
	double W; // longitude of ascending node
	double w; // argument of periapsis (longitude of periapsis is W + w)
	double p; // phase angle
	
	Body secondary;
	
	// Vis-Viva Equation
	// v^2 = GM(2/r - 1/a)
	// r = 1/(v^2 / (2GM) + 1/(2a))
	
	// r_apo/r_peri = (1 + e)/(1 - e)
	
	// r_peri = a(1 - e)
	// r_apo = a(1 + e)
	
	// Orbit equation (polar form)
	// r^-3 = -theta'^2/(G M) * 1/(1 + e cos(theta))
	// r = -cbrt(GM/(theta')^2 * (1 + e cos(theta)))
	// cos(theta) = x
	// theta = arccos(x)
	// theta' = -x'/sqrt(1-x^2)
	
	// How to get theta for a particular r?
	// r = (h^2 / GM) / (1 + e cos theta)
	// where h = r v cos phi, where phi is the "flight path angle"...
	// and h is a conserved quantity?
	// or alternatively h = sqrt(G M l), l = b^2 / a
	// recall e = sqrt(1 - b^2 / a^2)
	// so a^2 (1 - e^2) = b^2
	// therefore: h = sqrt(G M a (1 - e^2))
	// r = a (1 - e^2)/(1 + e cos(theta))
	// ...Oh.
	// theta = arccos((a(1 - e^2)/r - 1)/e)
	
	// How does time work into this?
	// Can calculate dx/dt using velocity, and that can get dtheta/dt
	// But that doesn't really give you an equation that you can just plug time into and get position, does it?
	
	/**
	 * Creates a two-body orbit with the given parameters.
	 * @param primary - main body, and the one described by these parameters
	 * @param secondary - the body the primary is orbiting around
	 * @param t - axial tilt
	 * @param tp - tilt phase
	 * @param rr - rotation rate
	 * @param rp - rotation phase
	 * @param a - semi-major axis (m)
	 * @param e - eccentricity (0 - 1)
	 * @param i - inclination (radians)
	 * @param W - longitude of ascending node (radians)
	 * @param w - argument of periapsis (radians)
	 * @param p - phase (radians)
	 */
	public Orbit(Body primary, Body secondary, double at, double tp, double rr, double rp, double a, double e, double i, double W, double w, double p) {
		this.primary = primary;
		this.secondary = secondary;
		this.at = at;
		this.tp = tp;
		this.rr = rr;
		this.rp = rp;
		this.a = a;
		this.e = e;
		this.i = i;
		this.W = W;
		this.w = w;
		this.p = p;
	}

	public double getVelocity(double r) {
		return Math.sqrt(G * M() * (2/r - 1/a));
	}
	
	public double getPeriapsis() {
		return a*(1-e);
	}
	
	public double getApoapsis() {
		return a*(1+e);
	}
	
	public double getPeriapsisVelocity() {
		return getVelocity(getPeriapsis());
	}
	
	public double getApoapsisVelocity() {
		return getVelocity(getApoapsis());
	}
	
	public double getPeriod() {
		// Hopefully m is small relative to M.
		return 2*Math.PI*Math.sqrt(a*a*a/(G*(M()+m())));
	}
	
	public Point3D<Double> getPositionAtTime(double t) {
		Point3D<Double> pos = Point3D.ORIGIN;
		//double[] pos = {0.0, 0.0, 0.0};
		double T = getPeriod();
		double mean_anomaly = 2*Math.PI*t/T + p; //equal to the phase angle at time 0, and otherwise increases by 2pi each orbit at a constant rate
		// Calculate eccentric anomaly E, with mean anomaly M = E - e sin E
		// Newton's method
		double ecc_anomaly = mean_anomaly;
		// Might do a fixed number of iterations instead because this should converge quickly.
		double iterations = 0;
		while(true) {
			iterations++;
			double delta = (ecc_anomaly - e * Math.sin(ecc_anomaly) - mean_anomaly)/(1 - e * Math.cos(ecc_anomaly));
			ecc_anomaly -= delta;
			if (Math.abs(delta) < 1E-6) break;
		}
		System.out.println("Converged in "+iterations+" iterations.");
		// P and Q give a coordinate system in the plane of the orbit, with P pointing toward periapsis
		double P = a * (Math.cos(ecc_anomaly) - e);
		double Q = a * Math.sin(ecc_anomaly) * Math.sqrt(1 - e*e);
		Point3D<Double> PQ0 = new Point3D<Double>(P, Q, 0.0);
		// Apply rotations to get this into a 3D coordinate system
		pos = PQ0.rotateXY(w);
		//pos[0] = Math.cos(w)*P - Math.sin(w)*Q;
		//pos[1] = Math.sin(w)*P + Math.cos(w)*Q;
		pos = pos.rotateYZ(i);
		//pos[2] = Math.sin(i)*pos[1];
		//pos[1] = Math.cos(i)*pos[1]; //TODO: This seems wrong.
		pos = pos.rotateXY(W);
		//double tmp = pos[0];
		//pos[0] = Math.cos(W)*tmp - Math.sin(W)*pos[1];
		//pos[1] = Math.sin(W)*tmp + Math.cos(W)*pos[1]; // TODO: Again, this seems wrong.
		return pos;
	}
	
	public Point3D<Double> getAxisAtTime(double t) {
		// Start with the normal axial tilt
		Point3D<Double> axis = new Point3D<Double>(Math.sin(at), 0.0, Math.cos(at));
		//double[] axis = {Math.sin(at), 0.0, Math.cos(at)};
		// Rotate in xy to account for the tilt phase
		axis = axis.rotateXY(tp);
		//axis[1] = Math.sin(tp)*axis[0];
		//axis[0] = Math.cos(tp)*axis[0];
		// Then apply rotations to account for the orientation of the orbit
		axis = axis.rotateXY(w);
		//double tmp = axis[0];
		//axis[0] = Math.cos(w)*tmp - Math.sin(w)*axis[1];
		//axis[1] = Math.sin(w)*tmp + Math.cos(w)*axis[1];
		axis = axis.rotateYZ(i);
		//tmp = axis[1];
		//axis[1] = Math.cos(i)*tmp - Math.sin(i)*axis[2];
		//axis[2] = Math.sin(i)*tmp + Math.cos(i)*axis[2];
		axis = axis.rotateXY(W);
		//tmp = axis[0];
		//axis[0] = Math.cos(W)*tmp - Math.sin(W)*axis[1];
		//axis[1] = Math.sin(W)*tmp + Math.cos(W)*axis[1];
		return axis;
	}
	
	public double getRotationAtTime(double t) {
		return rp + rr*t;
	}
	
	public double getSiderealDayLength() {
		return 2*Math.PI/rr;
	}
	
	public double getSiderealDays() {
		double T = getPeriod();
		return T/getSiderealDayLength();
	}
	
	public double getSolarDays() {
		return getSiderealDays() - 1;
	}
	
	public double getSolarDayLength() {
		// Not as computationally efficient as it could be, but more readable?
		double siderealDayLength = getSiderealDayLength();
		double siderealDays = getSiderealDays();
		double solarDays = getSolarDays();
		return siderealDayLength * siderealDays/solarDays;
	}
	
	//Earth sidereal and solar day
	//86164 s 
	//86400 s
	
	// Getters and setters
	
	public double getSemimajorAxis() {
		return a;
	}

	public void setSemimajorAxis(double a) {
		this.a = a;
	}

	public double getEccentricity() {
		return e;
	}

	public void setEccentricity(double e) {
		this.e = e;
	}

	public double M() {
		return secondary.getMass();
	}

	public double m() {
		return primary.getMass();
	}

	public double getInclination() {
		return i;
	}

	public void setInclination(double i) {
		this.i = i;
	}

	public double getLongitudeOfAscendingNode() {
		return W;
	}

	public void setLongitudeOfAscendingNode(double W) {
		this.W = W;
	}

	public double getArgumentOfPeriapsis() {
		return w;
	}

	public void setArgumentOfPeriapsis(double w) {
		this.w = w;
	}

	public double getPhase() {
		return p;
	}

	public void setPhase(double p) {
		this.p = p;
	}
}
