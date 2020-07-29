
public class OrbitCalculator {
	public static void main(String[] args) {
		Planet Earth = new Planet(5.9722E24);
		Star Sun = new Star(1.98847E30);
		Orbit test = new Orbit(Earth, Sun, Math.PI*23.5/180.0, 0.0, 7.292115E-5, 0.0, 1.4959787E11, 0.01671022, 0.0, 0.0, 0.0, 0.0);
		double period = test.getPeriod();
		System.out.println(period);
		Point3D<Double> pos = test.getPositionAtTime(period*0.25);
		System.out.println(pos);
		//System.out.println("x: "+pos[0]);
		//System.out.println("y: "+pos[1]);
		//System.out.println("z: "+pos[2]);
		double dist = pos.dist(); //Math.sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
		System.out.println(dist);
		double solarDayLength = test.getSolarDayLength();
		double solarDays = test.getSolarDays();
		System.out.println("Solar day length: "+solarDayLength);
		System.out.println("Solar days: "+solarDays);
		System.out.println(convertTime(solarDayLength));
		double siderealDayLength = test.getSiderealDayLength();
		System.out.println(convertTime(siderealDayLength));
		//Seems to be working okay. Some slight error because the values for masses and distances are not exact.
	}
	
	public static String convertTime(double t) {
		int d = (int)(t/86400);
		int h = (int)((t % 86400)/3600);
		int m = (int)((t % 3600)/60);
		int s = (int)(t % 60);
		//int ms = (int)((1000*t)%1000);
		return d+" d "+h+" h "+m+" m "+s+" s";
	}
}
