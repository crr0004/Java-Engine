package me.tempus.skeleton;

import me.tempus.camera.Camera;
import me.tempus.shader.PVM;
import me.tempus.shader.VI_Shader;

import org.lwjgl.LWJGLException;
import org.lwjgl.input.Keyboard;
import org.lwjgl.opengl.Display;
import org.lwjgl.opengl.DisplayMode;
import org.lwjgl.opengl.GL11;
import org.lwjgl.opengl.GL20;

public class Driver {


	private boolean done = false;
	
	public VI_Shader shader;
	
	final Camera camera = new Camera();
	
	private void createScreen(int width, int height){
		try {
	
			//Display.setFullscreen(fullscreen);
	
			Display.setDisplayMode(new DisplayMode(width, height));
			Display.setTitle("Camera");
				Display.create();
			
			GL11.glViewport(0, 0, width, height);
	
			GL11.glEnable(GL11.GL_TEXTURE_2D);                          // Enable Texture Mapping
			
			GL11.glShadeModel(GL11.GL_SMOOTH);                          // Enables Smooth Color Shading
			GL11.glClearColor(0.5f, 0.5f, 0.5f, 0.0f);                // This Will Clear The Background Color To Black
			GL11.glClearDepth(1.0);                                   // Enables Clearing Of The Depth Buffer
			GL11.glEnable(GL11.GL_DEPTH_TEST);                          // Enables Depth Testing
			GL11.glDepthFunc(GL11.GL_LEQUAL);                           // The Type Of Depth Test To Do
			PVM.setUpProjection(45f, width, height, 0.1f, 100f);
			System.out.println("OpenGL version: " + GL11.glGetString(GL11.GL_VERSION));
		} catch (LWJGLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
		
	private void initOPG(){
		Keyboard.enableRepeatEvents(true);
	}
	
	public void loop(){
		
		createScreen(1024, 1024);
		initOPG();
		setUpShaders();
		
		while(!done){
			update();
			draw();
			
		}
		return;
	}
	
	private void setUpShaders(){
		shader = new VI_Shader("TexturedShader.vert", "TexturedShader.frag");
	}
		
	private void update(){
			
		camera.pollInput();
			
	}
		
	private void draw(){
			
		GL20.glUseProgram(shader.getShaderID());
		GL11.glClear(GL11.GL_COLOR_BUFFER_BIT | GL11.GL_DEPTH_BUFFER_BIT);
		
		PVM.loadIdentity();
		camera.transform();
		
		
		GL20.glUseProgram(0);
		Display.sync(60);
		Display.update();
			
	}
		
		
		
	public static void main(String[] args){
			(new Driver()).loop();
	}
}
