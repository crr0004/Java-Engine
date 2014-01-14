package me.tempus.gameobjects;

import java.io.IOException;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;

import me.tempus.shader.PVM;

import org.lwjgl.BufferUtils;
import org.lwjgl.opengl.GL11;
import org.lwjgl.opengl.GL15;
import org.lwjgl.opengl.GL20;
import org.lwjgl.opengl.GL30;
import org.lwjgl.util.vector.Vector3f;
import org.newdawn.slick.opengl.Texture;
import org.newdawn.slick.opengl.TextureLoader;
import org.newdawn.slick.util.ResourceLoader;

public class Plane {

	private int vboiId;
	private int vaoId;
	
	private Vector3f pos;
	private Vector3f scale;
	private Texture texture;
	
	public Plane(Vector3f pos, Vector3f scale){
		this.pos = pos;
		this.scale = scale;
	}
	
	public void draw(){
		
		PVM.pushMatrix();
		PVM.rotate(45, 0, 0);
		PVM.translate(pos.x, pos.y, pos.z);
		PVM.scale(scale.x, scale.y, scale.z);
		PVM.bufferMatricies();
		PVM.popMatrix();
		
		GL11.glBindTexture(GL11.GL_TEXTURE_2D, texture.getTextureID());
		
		GL30.glBindVertexArray(vaoId);
		GL20.glEnableVertexAttribArray(0); // Again the reference to the shader bound reference
		GL20.glEnableVertexAttribArray(1);
		
		GL15.glBindBuffer(GL15.GL_ELEMENT_ARRAY_BUFFER, vboiId);
		
		GL11.glDrawElements(GL11.GL_TRIANGLES, 2 * 3, GL11.GL_UNSIGNED_INT, 0); // Amount is how many triangles times by the vertex size
		GL20.glDisableVertexAttribArray(1);
		GL20.glDisableVertexAttribArray(0);
		GL30.glBindVertexArray(0);
		
		
	}
	
	public void create(String texturePath){
		
		float[] vertices = new float[]{
				-1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, -1.0f, -1.0f, 0.0f, -1.0f
		};
		
		float[] textureData = new float[]{
				0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f
		};
		
		int[] indicies = new int[]{
			0, 1, 2,
			0, 2, 3
		};
		
		FloatBuffer verticesBuffer = BufferUtils.createFloatBuffer(vertices.length);
		verticesBuffer.put(vertices);
		verticesBuffer.flip();
		
		FloatBuffer textureBuffer = BufferUtils.createFloatBuffer(textureData.length);
		textureBuffer.put(textureData);
		textureBuffer.flip();		
		
		IntBuffer indicesBuffer = BufferUtils.createIntBuffer(indicies.length);
		indicesBuffer.put(indicies);
		indicesBuffer.flip();		
		
		vaoId = GL30.glGenVertexArrays();
		GL30.glBindVertexArray(vaoId);
		
		
		//Bind and buffer the vertex data
		final int vboId = GL15.glGenBuffers();
		GL15.glBindBuffer(GL15.GL_ARRAY_BUFFER, vboId);
		GL15.glBufferData(GL15.GL_ARRAY_BUFFER, verticesBuffer, GL15.GL_STATIC_DRAW);
		GL20.glVertexAttribPointer(0, 3, GL11.GL_FLOAT, false, 0, 0); // First 0 is a reference to the shader binding location
		GL15.glBindBuffer(GL15.GL_ARRAY_BUFFER, 0);
		
		final int vbotId = GL15.glGenBuffers();
		GL15.glBindBuffer(GL15.GL_ARRAY_BUFFER, vbotId);
		GL15.glBufferData(GL15.GL_ARRAY_BUFFER, textureBuffer, GL15.GL_STATIC_DRAW);
		GL20.glVertexAttribPointer(1, 2, GL11.GL_FLOAT, false, 0, 0);
		GL15.glBindBuffer(GL15.GL_ARRAY_BUFFER, 0);	
		
		
		GL30.glBindVertexArray(0);
		
		//Bind and buffer the indices Data
		vboiId = GL15.glGenBuffers();
		GL15.glBindBuffer(GL15.GL_ELEMENT_ARRAY_BUFFER, vboiId);
		GL15.glBufferData(GL15.GL_ELEMENT_ARRAY_BUFFER, indicesBuffer, GL15.GL_STATIC_DRAW);
		GL15.glBindBuffer(GL15.GL_ELEMENT_ARRAY_BUFFER, 0);
		
		try {
			texture = TextureLoader.getTexture("JPG", ResourceLoader.getResourceAsStream(texturePath));
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	public void setPos(Vector3f planePos) {
		// TODO Auto-generated method stub
		this.pos = planePos;
	}
	
	
}
