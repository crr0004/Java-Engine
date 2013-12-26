/**
 * 
 */
package me.tempus.math;

/**
 * @author Chris
 *
 */
public class RootFinding {

	/**
	 * 
	 * @param c The coefficients of each term in the equation
	 * @return
	 * a = c[1]
	 * b = c[2]
	 * c = c[3]
	 */
	public static double[] solveCubic2(float[] c){
		final double[] returnValues = new double[3];		
		final int size = c.length;
		final float a3 = c[1] /3;
		
		if(c[0] != 1){
			final double c0 = c[0];
			for(int i = 0; i < 4; i++){
				c[i] /= c0;
			}
		}
		
		final float p = (float) (Math.pow((c[1]), 2) /-3.0) + c[2];
		final float q = (float) (2*Math.pow(c[1], 3) /27.0) - ((c[1]) * c[2])/3 + c[3];
		
		//final float d = (-4 * (p*p*p)) - (27 * (p * p));
		
		final float pDash = p/3;
		final float qDash = q/2;
		
		final float d = -108 * ((pDash*pDash*pDash) + (qDash*qDash));
		final float dDash = d/108;
		
		final float r =  pow((-qDash + sqrt(-dDash)), 1.0/3);
		final float s =  pow((-qDash - sqrt(-dDash)), 1.0/3);
		
		if(dDash < 0){
			returnValues[0] = (r + s) - a3;
		}else if(dDash == 0){
			returnValues[0] = (2 * r) -a3;
			returnValues[1] = (-r) -a3;
		}else if(dDash > 0){
			final double m = Math.sqrt(-p/3);
			final double theta = Math.acos(-qDash / 2 * (m*m*m)) /3;
			final double mDash = 2 * Math.sqrt(-pDash);
			final double radianConst = (2 * Math.PI) /3;
			returnValues[0] = ((float) (mDash * Math.cos(theta))) -a3;
			returnValues[1] = ((float) (mDash * Math.cos(theta + radianConst))) -a3;
			returnValues[2] = ((float) (mDash * Math.cos(theta - radianConst))) -a3;
		}
		
		return returnValues;
	}
	
	public static double[] solveCubic(double[] c){
		final double[] ret = new double[3];
		
		if(c[0] != 1){
			final double c0 = c[0];
			for(int i = 0; i < 4; i++){
				c[i] /= c0;
			}
		}
		
		final double a3 = c[1]/3;
		
		final double p = (c[2]/3) - ((c[1] * c[1]) /9);
		final double q = ((c[1]*c[1]*c[1])/27) - ((c[1] * c[2]) /6) + (c[3] /2);
		final double disc = (p*p*p) + (q*q);
		
		if(disc >= 0){
			final double sqrtDisc = Math.sqrt(disc);
			final double r = Math.pow(-q + sqrtDisc, 1/3.0);
			if(disc == 0){
				ret[0] = (2*r) - a3;
				ret[1] = (-r) - a3;
				return ret;
			}else{
				final double s = Math.pow(-q - sqrtDisc, 1/3.0);
				ret[0] = (r+s) - a3;
				return ret;
			}
		}else{
			final double ang = Math.acos(-q / Math.sqrt((-p*p*p)));
			final double r = 2 * Math.sqrt(-p);
			final double twoPi = 2 * Math.PI;
			for(int i = -1; i < 2; i++){
				ret[i +1] = r * Math.cos((ang - twoPi*i) /3);
				ret[i +1] -= a3;
				
			}
			return ret;
		}
	}
	
	public static float pow(double d, double e){
		if(d < 0){
			return (float) -Math.pow(Math.abs(d), e);
			
		}
		return (float) Math.pow(d, e);
	}
	
	public static double sqrt(double a){
		if(a < 0){
			return -Math.sqrt(Math.abs(a));
		}else{
			return Math.sqrt(a);
		}
	}
	
}
