package me.tempus.shader;

import static org.lwjgl.opengl.GL11.GL_FALSE;

import org.lwjgl.opengl.GL11;
import org.lwjgl.opengl.GL20;

public class VIT_Shader {
	
	protected int shaderID;
	protected int vertID;
	protected int fragID;
	
	private int vVertex = 0;
	
	private int vTexture = 1;
	
	/**
	 * <h1>Shader Constuctor </h1>
	 * 
	 * @param vertPath The path to the vertex shader, (Full path, File.separatorChar for \, and extension)
	 * @param fragPath The path to the fragment shader, (Full path, File.separatorChar for \, and extension)
	 */
	
	public VIT_Shader(){

	}
	
	/**
	 * <h1>Pass the uniform variable data to the shader</h1>
	 * Use GL20.glUniform<br>
	 * <b>Must be overriden to have code</b>
	 */
	
	public void passUniforms(float[] ambientColor, float[] diffuseColor, float[] specularColor, float shine){

	}
	
	
	/**
	 * <h1>Set the uniform locations to variables</h1>
	 * Use GL20.glGetUniformLocation(program, name)<br>
	 * <b>Must be overriden to have code</b>
	 */
	
	public void getUniformLocations(){

	}
	
	/**
	 * <h1> Bind attribute locations (in variables)</h1>
	 * Use GLint glGetAttribLocation(int program, char *name);<br>
	 * Called from within loadShaders<br>
	 * <b>Must be overriden to have code</b>
	 */
	
	public void bindAttributeLocations(){
		GL20.glBindAttribLocation(shaderID, vVertex, "vVertex");
		GL20.glBindAttribLocation(shaderID, vTexture, "vTexture");
	}
	
	/**
	 * <h1>Load the Shaders</h1>
	 * Loads, complies, sets up attribute locations and links shader<br>
	 * <b>Must be called after the setup of OpenGL</b>
	 */
	
	public void load(String vertPath, String fragPath){
		
		shaderID = GL20.glCreateProgram();
		vertID = GL20.glCreateShader(GL20.GL_VERTEX_SHADER);
		fragID = GL20.glCreateShader(GL20.GL_FRAGMENT_SHADER);

		StringBuilder[] shaderSources = ShaderUtilities.loadShaders(vertPath, fragPath);

		GL20.glShaderSource(vertID, shaderSources[0]);

		GL20.glCompileShader(vertID);

		if((GL20.glGetShaderi(vertID, GL20.GL_COMPILE_STATUS)) == GL11.GL_FALSE){
			System.err.println("Compling of vertex shader failed");

		}


		GL20.glShaderSource(fragID, shaderSources[1]);
		GL20.glCompileShader(fragID);

		if((GL20.glGetShaderi(fragID, GL20.GL_COMPILE_STATUS)) == GL11.GL_FALSE){
			System.err.println("Compling of fragment shader failed");
		}

		bindAttributeLocations();
		
		GL20.glAttachShader(shaderID, vertID);
		GL20.glAttachShader(shaderID, fragID);

				
		
		GL20.glLinkProgram(shaderID);
		GL20.glValidateProgram(shaderID);
		
		int projectionMatrixLocation = GL20.glGetUniformLocation(shaderID, "projection_matrix");
		int viewMatrixLocation = GL20.glGetUniformLocation(shaderID, "view_matrix");
		int modelMatrixLocation = GL20.glGetUniformLocation(shaderID, "model_matrix");

		PVM.setLocations(projectionMatrixLocation, viewMatrixLocation, modelMatrixLocation);
		getUniformLocations();
		
		if((GL20.glGetProgrami(shaderID, GL20.GL_LINK_STATUS) == GL_FALSE)){
				
			System.out.println("Link failed");
			System.out.println(GL20.glGetShaderInfoLog(fragID, 10000));
			System.out.println(GL20.glGetShaderInfoLog(vertID, 10000));
			System.out.println(GL20.glGetShaderInfoLog(shaderID, 10000));
		}
	}
	
	public int getShaderID(){
		return shaderID;
	}
	
	public int getFragID(){
		return fragID;
	}
	
	public int getVertID(){
		return vertID;
	}
	

	
}
