package me.tempus.math;

public class Vector3f {

	public float x;
	public float y;
	public float z;
	
	
		
	/**
	 * @param x
	 * @param y
	 * @param z
	 */
	public Vector3f(float x, float y, float z) {
		super();
		this.x = x;
		this.y = y;
		this.z = z;
	}
	/**
	 * @return the x
	 */
	public float getX() {
		return x;
	}
	/**
	 * @param x the x to set
	 */
	public void setX(float x) {
		this.x = x;
	}
	/**
	 * @return the y
	 */
	public float getY() {
		return y;
	}
	/**
	 * @param y the y to set
	 */
	public void setY(float y) {
		this.y = y;
	}
	/**
	 * @return the z
	 */
	public float getZ() {
		return z;
	}
	/**
	 * @param z the z to set
	 */
	public void setZ(float z) {
		this.z = z;
	}

	@Override
	public String toString() {
		return String.format("X: %s, Y: %s, Z: %s", x, y, z);
	}

}
