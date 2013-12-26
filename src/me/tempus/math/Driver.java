package me.tempus.math;

/**
 * A class used for testing
 * @author Chris
 *
 */
public class Driver {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		/*
		double[][] result = MatrixMath.solveEigenVectors(new double[][]{
				{2, 1, -1},
				{-3, -1, 2},
				{-2, 1, 2}
		}, new double[]{5.14005,-2.14005,1});
		*/
		//System.out.println(result);
		
		double[][] m = new double[][]{
				{2 - 3.21432, 1, -1},
				{-3, -1 - 3.21432, 2},
				{-2, 1, 2 - 3.21432}
		};
		MatrixMath.solveLUDecomp(MatrixMath.LUDecomp(m, 3), m, new double[]{0, 0, 0});
	}
}
