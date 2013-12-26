/**
 * 
 */
package me.tempus.math;



/**
 * @author Chris
 *
 */
public class PCA {

	/**
	 * The average of a list of vectors
	 * @param vectors
	 * @return A vector containing the average of all the vectors
	 */
	public static Vector3f vectorAverage(Vector3f[] vectors){
		float xSum = 0;
		float ySum = 0;
		float zSum = 0;
		final int size = vectors.length;
		for(Vector3f vector: vectors){
			xSum += vector.x;
			ySum += vector.y;
			zSum += vector.z;
		}		
		
		return new Vector3f(xSum / size, ySum /size, zSum /size);
	}
	
	/**
	 * Setups a covariance matrix for a set of vectors
	 * @param vectors The vectors to up the Covariance Matrix with
	 * @param vectorAverage The average of the vectors
	 * @return Float Covariance Matrix
	 */
	public static float[] setUpCovarianceMatrix(Vector3f[] vectors, Vector3f vectorAverage){
		final float[] matrix = new float[9];
		
		/**
		 * Setup the matrix values
		 */
		final int size = vectors.length;
		for(Vector3f vector : vectors){
			matrix[0] += square(vector.x - vectorAverage.x); // 1x1
			matrix[1] += (vector.x - vectorAverage.x) * (vector.y - vectorAverage.y); // 1x2
			matrix[2] += (vector.x - vectorAverage.x) * (vector.z - vectorAverage.z); // 1x3
			matrix[4] += square(vector.y - vectorAverage.y); //2x2
			matrix[5] += (vector.y - vectorAverage.y) * (vector.z - vectorAverage.z); // 2x3
			matrix[8] += square(vector.z - vectorAverage.z); // 3x3
		}
		matrix[3] = matrix[1];
		matrix[6] = matrix[2];
		matrix[7] = matrix[5];
		
		/**
		 * Average the matrix
		 */
		for(int i = 0; i < matrix.length; i++){
			matrix[i] /= size;
		}
		return matrix;
	}
	
	public static float square(float value){
		return value * value;
	}
	
}
