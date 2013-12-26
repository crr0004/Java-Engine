/**
 * 
 */
package me.tempus.math;

/**
 * @author Chris
 *
 */
public class MatrixMath {

	private static final double EPSILON = 1e-10;
	
	/**
	 * Solves the eigenVectors of a matrix given its eigenvalues
	 * Notes: http://www.sosmath.com/matrix/eigen2/eigen2.html
	 * @param matrix Row major
	 * @param eigenValues The eigenvalues of the matrix
	 * @return At most 3 vectors double[i] is the ith eigenvector, such that [i][j] is the eigenvectors jth element
	 */
	public static double[][] solveEigenVectors(double[][] matrix, double[] eigenValues){
		
		final double[][] results = new double[matrix.length][matrix[0].length];
		
		for(int i = 0; i < results.length; i++){
			double[][] m = new double[matrix.length][matrix[0].length];
			
			//System.arraycopy(matrix, 0, m, 0, matrix.length*matrix[0].length);
			for(int j = 0; j < matrix.length; j++){
				for(int k = 0; k < matrix[0].length; k++){
					m[j][k] = matrix[j][k];
				}
			}
			
			
			final double eigenValue = eigenValues[i];
			m[0][0] -= eigenValue;
			m[1][1] -= eigenValue;
			m[2][2] -= eigenValue;
			
			results[i] = solveLUDecomp(LUDecomp(m, m.length), m, new double[]{0,0,0});
			
			//System.out.println(result);
		}
		
		
		return results;		
	}
	
	/**
	 * Solves for the coefficients of the characteristic polynomial
	 * Uses: http://en.wikipedia.org/wiki/Characteristic_polynomial#Characteristic_equation
	 * @param matrix Matrix to find the polynomial for
	 * @return The coefficients of the characteristic polynomial of a 3x3 matrix.
	 * The order of the array is respective to x^3, x^2, x^1 and the constant
	 */
	public static float[] cp3(float[] matrix){
		final float[] values = new float[4];
		values[0] = -1;
		final float trace = trace3(matrix);
		values[1] = trace;
		values[2] = -((trace * trace) - (trace3(matrixSquare(matrix, 3, 3))))/2;
		values[3] = det3(matrix);
		return values;
	}

	/**
	 * Determinate for a 3x3 matrix
	 * @param m Matrix to solve for
	 * @return Scalar value of the determinate
	 */
	public static float det3(float[] m){
		return (m[0]*m[4]*m[8]) - (m[0]*m[5]*m[7]) - (m[1]*m[3]*m[8]) + (m[1]*m[5]*m[6]) + (m[2]*m[3]*m[7]) - (m[2]*m[4]*m[6]);
	}
	/**
	 * Determinate for a 3x3 matrix
	 * @param m Matrix to solve for
	 * @return Scalar value of the determinate
	 */
	public static double det3(double[] m){
		return (m[0]*m[4]*m[8]) - (m[0]*m[5]*m[7]) - (m[1]*m[3]*m[8]) + (m[1]*m[5]*m[6]) + (m[2]*m[3]*m[7]) - (m[2]*m[4]*m[6]);
	}
	
	public static double det3(double[][] m){
		return (m[0][0]*m[1][1]*m[2][2]) + (m[0][1]*m[1][2]*m[2][0]) + (m[0][2]*m[1][0]*m[2][1]) - (m[0][2]*m[1][1]*m[2][0]) - (m[0][1]*m[1][0]*m[2][2]) - (m[0][0]*m[1][2]*m[2][1]);
	}
	
	/**
	 * Sqaure a matrix
	 * @param matrix
	 * @param n Row size
	 * @param m Column Size
	 * @return The square of the matrix
	 */
	public static float[] matrixSquare(float[] matrix, int n, int m){
		final float[] result = new float[matrix.length];
		
		for(int i = 0; i < n; i++){
			for(int j = 0; j < m; j++){
				result[getPosistion(i, j, m)] = getMatrixSum(matrix, matrix, m, m, i, j);
			}			
		}
		
		return result;
	}
	
	/**
	 * Raise a matrix to the power of 3
	 * @param matrix Matrix
	 * @param pow The power to raise the matrix to
	 * @param n Row size of the matrix
	 * @param m Column size of the matrix
	 * @return The product matrix
	 */
	public static float[] matrixPow3(float[] matrix, int n, int m){
		float[] result = new float[matrix.length];

		for(int i = 0; i < n; i++){
			for(int j = 0; j < m; j++){
				result[getPosistion(i, j, m)] = getMatrixSum(matrix, matrix, m, m, i, j);
			}			
		}
		final float[] resultTwo = new float[matrix.length];
		for(int i = 0; i < n; i++){
			for(int j = 0; j < m; j++){
				resultTwo[getPosistion(i, j, m)] = getMatrixSum(result, matrix, m, m, i, j);
			}			
		}
		return resultTwo;
	}
		/**
		 * Trace of a nxn size matrix
		 * @param matrix A square matrix
		 * @param n Size of the matrix
		 * @return The trace of the matrix
		 */
		public static float trace(float[] matrix, int n){
			float trace = 0;
			for(int i = 0; i < n; i++){
				trace += matrix[getPosistion(i, i, n)];
			}
			return trace;
		}

		/**
		 * @param matrix
		 * @return Trace for a 3x3 matrix
		 */
		public static float trace3(float[] matrix){
			return matrix[0] + matrix[4] + matrix[8];
		}

		/**
		 * @param matrix
		 * @return Trace for a 2x2 matrix
		 */
		public static float trace2(float[] matrix){
			return matrix[0] + matrix[3];
		}

		/**
		 * Returns a new 3x3 identity matrix
		 * @return identity matrix
		 */
		public static float[] identityMatrix3(){
			return new float[]{
				1,0,0,
				0,1,0,
				0,0,1
			};
			
		}
		
		/**
		 * Multiples a matrix by a scalar value
		 * @param matrix
		 * @param scalar
		 * @param m Row size of the matrix
		 * @param n Column size of the matrix
		 * @return The resulting matrix
		 */
		public static float[] scalarMultiple(float[] matrix, float scalar, int m, int n){
			float[] result = new float[matrix.length];
			
			for(int i = 0; i < m; i++){
				for(int j = 0; j < n; j++){
					result[getPosistion(i, j, n)] = matrix[getPosistion(i, j, n)] * scalar;
				}
			}
			
			return result;
		}
		
		/**
		 * Subtracts one square matrix from another (m1 - m2)
		 * @param m1 Matrix
		 * @param m2 Matrix
		 * @param m Row size of resulting matrix
		 * @param n Column size of result matrix
		 * @return The result matrix
		 */
		public static float[] subtractMatrix(float[] m1, float[] m2, int m, int n){
			float[] result = new float[m1.length];
			
			for(int i = 0; i < m; i++){
				for(int j = 0; j < n; j++){
					result[getPosistion(i, j, n)] = m1[getPosistion(i, j, n)] - m2[getPosistion(i, j, n)];
				}
			}
			
			return result;
		}
		
		/**
		 * TODO Test
		 * @param matrix1 Matrix 1
		 * @param matrix2 Matrix 2
		 * @param n Matrix1 row size
		 * @param m Matrix1 column size
		 * @param o Matrix2 row size
		 * @param p Matrix2 column size
		 * @return The product matrix
		 */
		public static float[] mulitpleMatrix(float[] matrix1, float[] matrix2, int n, int m, int o, int p){
			assert m == o;
			final float[] f = new float[m*p];

			for(int i = 0; i < m; i++){
				for(int j = 0; j < p; j++){
					f[getPosistion(i, j, p-1)] = getMatrixSum(matrix1, matrix2, m, p, i, j);
				}
			}

			return f;
		}

		/**
		 * Used in multiplying matrices. Basically the cross product
		 * @param matrix1
		 * @param matrix2
		 * @param m Matrix1 column size
		 * @param p Matrix2 column size
		 * @param i Current row position of the product matrix
		 * @param j Current column position of the product matrix
		 * @return
		 */
		public static float getMatrixSum(float[] matrix1, float[] matrix2, int m, int p, int i, int j){
			float returnValue = 0;
			for(int k = 0; k < m; k++){
				returnValue += matrix1[getPosistion(i, k, m)] * matrix2[getPosistion(k, j, p)];
			}
			return returnValue;
		}

		/**
		 * 
		 * @param m nxn matrix
		 * @param n size of the matrix
		 * @return L and U such that LU=m
		 */
		public static double[][] LUDecomp(double[][] m, int n){
			
			if(det3(m) == 0){
				//throw new RuntimeException("LU Decompostion can't handle singular matricies");
			}
			
			//http://www.unc.edu/~marzuola/Math547_S13/Math547_S13_Projects/H_Xing_Section001_LeastSquaresCholesky.pdf
			//int n=m.length;
			double[][]L= new double[n][n];
			for(int k = 0; k < n-1; k++){ 
				
				L[k][k] = 1; // add diagonal entry 1
				
				for(int i=k+1; i < n; i++){
					
					L[i][k]= m[i][k] / m[k][k]; // update column of L 
					
					for(int j = k; j <n; j++)
						m[i][j] -= L[i][k] * m[k][j]; //row operation
				}
			}
			
			L[n-1][n-1] = 1;// don't forget the last diagonal entry 
			return L; 
			
		}

		/**
		 * Solves Ly = b, then y = Ux. Solves Ax = b
		 * @param m LU Decomposition matrix
		 * @param b Vector answer
		 * @return
		 */
		@Deprecated
		public static double[] solveLUDecomp(double[][] m, double[] b){
			double[] x = new double[3]; // Vector that is being solved for
			int n = m[0].length;
			
			//Use forward substitution to solve for y
			double[] y = new double[x.length];
			for(int i = 0; i < n; i++){
				double a = b[i];
				for(int j = 0; j < i-1; j++){
					a -= m[i][j]*y[j];
				}
				y[i] = a / m[i][i];
			}
			
			//Use back substitution to solve for x
			for(int i = n-1; i >= 0; i--){
				double a = y[i];
				for(int j = i+1; j < n; j++){
					a -= m[i][j] * x[j];
				}
				x[i] = a / m[i][i];
			}
			
			return x;
		}
		
		/**
		 * Solves Ly = b, then y = Ux. Solves Ax = b
		 * @param m LU Decomposition matrix
		 * @param b Vector answer
		 * @return
		 */
		public static double[] solveLUDecomp(double[][] L, double[][] U, double[] b){
			double[] x = new double[3]; // Vector that is being solved for
			int n = L[0].length;
			
			//Use forward substitution to solve for y
			double[] y = new double[x.length];
			for(int i = 0; i < n; i++){
				y[i] = b[i];
				for(int j = 0; j < i; j++){
					y[i] -= L[i][j]*y[j];
				}
				y[i] /= L[i][i];
			}
			
			
			
			//Use back substitution to solve for x
			for(int i = n-1; i >= 0; i--){
				
				if(U[2][2] <= EPSILON && i == n-1){
					x[i] = 1;
					i = n-2;
				}
				x[i] = y[i];
				for(int j = i+1; j < n; j++){
					x[i] -= U[i][j] * x[j];
				}
				x[i] /= U[i][i] ;
			}
			
			
			
			return x;
		}
		
		/**
		 * Reduced row echelon form
		 * @param matrix The augmented matrix to solve
		 * @param matrixSize The size the matrix's row
		 * @return Returns the solved matrix
		 */
		public static float[] rref2(float[] matrix, int matrixSize){
			int i = 0;
			int j = 0;

			for(;j < matrixSize; j++){
				float lastValue = 0;
				int largestRow = i;
				for(int k = i; k < matrixSize; k++){
					final int pos = getPosistion(k, j, matrixSize+1);
					if(Math.abs(matrix[pos]) > lastValue){
						lastValue = matrix[pos];
						largestRow = k;
					}
				}
				/**
				 * If found largest value in a row, swamp the row
				 */
				if(largestRow != i){
					matrix = swapRows(matrix, matrixSize +1, i, largestRow);
				}

				/**
				 * Multiple row i by 1/mij
				 */
				matrix = multipleRow(matrix, matrixSize+1, i, 1 / matrix[getPosistion(i, j, matrixSize+1)]);

				/**
				 * For each row r, where 1≤r≤n and r≠i, add -Mrj times row i it row r
				 */
				for(int r = 0; r < matrixSize; r++){
					if(r != i){
						final float mrj = -matrix[getPosistion(r, j, matrixSize+1)];
						for(int l = 0; l < matrixSize+1; l++){
							matrix[getPosistion(r, l, matrixSize+1)] += mrj * matrix[getPosistion(i, l, matrixSize+1)];
						}
					}
				}
				i++;
			}

			return matrix;
		}
		
		private static float[] swapRows(float[] matrix, int matrixSize, int row1, int row2) {
			
			for(int i = 0; i < matrixSize; i++){
				final double m1 = matrix[((matrixSize * row1) + i)];
				final double m2 = matrix[((matrixSize * row2) + i)];

				matrix[(matrixSize * row1) + i] = (float) m2;
				matrix[(matrixSize * row2) + i] = (float) m1;			
			}

			return matrix;
		}

		/**
		 * Sourced from: http://introcs.cs.princeton.edu/java/95linear/GaussianElimination.java.html
		 * @param A
		 * @param b
		 * @return
		 */
		@Deprecated
		public static double[][] rref(double[][] A, double[]b){
			int N  = b.length;

	        for (int p = 0; p < N; p++) {

	            // find pivot row and swap
	            int max = p;
	            for (int i = p + 1; i < N; i++) {
	                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
	                    max = i;
	                }
	            }
	            double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
	            double   t    = b[p]; b[p] = b[max]; b[max] = t;

	            // singular or nearly singular
	            if (Math.abs(A[p][p]) <= EPSILON) {
	               // throw new RuntimeException("Matrix is singular or nearly singular");
	            }

	            // pivot within A and b
	            for (int i = p + 1; i < N; i++) {
	                double alpha = A[i][p] / A[p][p];
	                b[i] -= alpha * b[p];
	                for (int j = p; j < N; j++) {
	                    A[i][j] -= alpha * A[p][j];
	                }
	            }
	        }
	        
	     // back substitution
	        double[] x = new double[N];
	        for (int i = N - 1; i >= 0; i--) {
	            double sum = 0.0;
	            for (int j = i + 1; j < N; j++) {
	                sum += A[i][j] * x[j];
	            }
	            x[i] = (b[i] - sum) / A[i][i];
	            System.out.println(x[i]);
	        }
	        
	        return A;
		}
		
		/**
		 * Reduced row echelon form
		 * Uses http://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode
		 * @param matrix The augmented matrix to solve
		 * @param matrixSize The size the matrix's row
		 * @return Returns the solved matrix
		 */
		public static double[] rref(double[] matrix, int matrixSize){
			int m = matrixSize;
			int n = matrixSize+1;
			
			if(det3(matrix) == 0){
				throw new RuntimeException("Det(A) is 0");
			}
			
			for(int k = 0; k < m; k++){
				//Find Pivot
				int pivot = k;
				for(int i = k; i < m; i++){
					if(Math.abs(matrix[getPosistion(i, k, n)]) > Math.abs(matrix[getPosistion(pivot, k, n)])){
						pivot = i;
					}
				}
				
				//swapRows(k, pivot);
				for(int i = 0; i < n; i++){
					final double m1 = matrix[((n * k) + i)];
					final double m2 = matrix[((n * pivot) + i)];

					matrix[(n * k) + i] = m2;
					matrix[(n * pivot) + i] = m1;			
				}
				
				for(int i = k+1; i < m; i++){
					for(int j = k; j < n; j++){
						matrix[getPosistion(i, j, n)] -= matrix[getPosistion(k, j, n)] * (matrix[getPosistion(i, k, n)] / matrix[getPosistion(k, k, n)]);
					}
					matrix[getPosistion(i, k, n)] = 0;
				}
			}

			return matrix;
		}

		/**
		 * Solves the 3x4 augmented matrix
		 * @param matrixToSolve Augmented matrix
		 * @return Vector containing the solution
		 */
		@Deprecated
		public static double[] solve3x3AugmentedMatrix(double[] matrixToSolve){
			rref(matrixToSolve, 3);
			
			//Vector3f result = new Vector3f(0, 0, 0);
			double[] result = new double[3];
			final int matrixSize = 4;
			//Use back solving
			if(matrixToSolve[getPosistion(2, 2, 4)] == 0){ //Prevents zeroing out
				result[2] = 1;
			}else{
				result[2] = matrixToSolve[getPosistion(2, 3, matrixSize)] / matrixToSolve[getPosistion(2, 2, matrixSize)];
			}
			//result[2] = 1;
			result[1] = (matrixToSolve[getPosistion(1, 3, matrixSize)] - (matrixToSolve[getPosistion(1, 2, matrixSize)] * result[2])) / matrixToSolve[getPosistion(1, 1, matrixSize)];
			result[0] = (matrixToSolve[getPosistion(0, 3, matrixSize)] - (matrixToSolve[getPosistion(0, 2, matrixSize)] * result[2]) - (matrixToSolve[getPosistion(0, 1, matrixSize)] * result[1])) / matrixToSolve[getPosistion(0, 0, matrixSize)];
			
			/*
			final int N = 3; //Vector size
			double[] result = new double[N];
	        for (int i = N - 1; i >= 0; i--) {
	            double sum = 0.0;
	            for (int j = i + 1; j < N; j++) {
	                sum += matrixToSolve[getPosistion(i, j, N+1)] * result[j];
	            }
	            result[i] = (matrixToSolve[getPosistion(i, 3, N+1)] - sum) / matrixToSolve[getPosistion(i, i, N+1)];
	        }
			*/
			return result;
		}
		
		/**
		 * Gets the last column of a matrix
		 * Useful for getting the result of an augmented matrix via rref
		 * @param matrix
		 * @param rowSize
		 * @param columnSize
		 * @return Last column of the matrix
		 */
		public static float[] getLastColumn(float[] matrix, int rowSize, int columnSize){
			//TODO Test
			final float[] newMatrix = new float[rowSize];

			for(int i = 0; i < rowSize; i++){
				newMatrix[i] = ((i+1) * columnSize) -1;
			}

			return newMatrix;
		}

		/**
		 * Values are zero indexed
		 * @param matrix
		 * @param matrixSize Column size of the matrix
		 * @param row1 Row one to be swaped with
		 * @param row2 Row two to swap to row one
		 * @return The matrix with the rows swamped
		 */
		public static double[] swapRows(double[] matrix, int matrixSize, int row1, int row2){

			for(int i = 0; i < matrixSize; i++){
				final double m1 = matrix[((matrixSize * row1) + i)];
				final double m2 = matrix[((matrixSize * row2) + i)];

				matrix[(matrixSize * row1) + i] = m2;
				matrix[(matrixSize * row2) + i] = m1;			
			}

			return matrix;
		}

		public static float[] multipleRow(float[] matrix, int columnSize, int row, double scalar){

			if (Double.isInfinite(scalar) || Double.isNaN(scalar)) {
				scalar = 0;
			}
			
			for(int i = 0; i < columnSize; i++){
				matrix[getPosistion(row, i, columnSize)] *= scalar;
			}

			return matrix;
		}

		public static float[] addRow(float[] matrix, int matrixSize, int row, float scalar){
			for(int i = 0; i < matrixSize; i++){
				matrix[(matrixSize * row) + i] += scalar;
			}

			return matrix;
		}

		/**
		 * All values expect matrixSize must be zero indexed
		 * @param i
		 * @param j
		 * @param matrixSize The column size of the matrix, <b> must <u>not</u> be zero indexed </b>
		 * @return The zero indexed posistion of the value in the matrix
		 */
		public static int getPosistion(int i, int j, int matrixSize){
			final int k = ((matrixSize * i) + j);
			return k;
		}
	}
