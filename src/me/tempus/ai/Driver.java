package me.tempus.ai;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import me.tempus.camera.Camera;
import me.tempus.collada.ColladaParser;
import me.tempus.gameobjects.Box;
import me.tempus.gameobjects.Geometry;
import me.tempus.shader.PVM;
import me.tempus.shader.VCI_Shader;

import org.lwjgl.LWJGLException;
import org.lwjgl.input.Keyboard;
import org.lwjgl.opengl.Display;
import org.lwjgl.opengl.DisplayMode;
import org.lwjgl.opengl.GL11;
import org.lwjgl.opengl.GL20;
import org.lwjgl.util.vector.Vector3f;

public class Driver {


	private boolean done = false;
	
	private VCI_Shader shader;
	
	private final Camera camera = new Camera();
	private List<Geometry> gameObjects = new ArrayList<Geometry>();
	
	private Box agent;
	private float agentMoveSpeed = 0.1f;
	private Vector3f desitination;
	private Vector3f direction;
	
	private void createScreen(int width, int height){
		try {
	
			//Display.setFullscreen(fullscreen);
	
			Display.setDisplayMode(new DisplayMode(width, height));
			Display.setTitle("Driver");
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
		createObjects();
		
		while(!done){
			update();
			draw();
			
		}
		return;
	}
	
	private void createObjects() {
		final HashMap<String, me.tempus.collada.Geometry> geometryMap = (new ColladaParser()).parse("src/me/tempus/ai/AI_Test_Scene.DAE");
		gameObjects = new ArrayList<Geometry>();
		int i = 0;
		Random r = new Random();
		for(me.tempus.collada.Geometry g : geometryMap.values()){
			final Box geometry = new Box();
			geometry.setName(g.getId());
			geometry.setIndices(g.getIndices());
			geometry.setVertices(g.getPosistions());
			geometry.setTriangleSize(g.getTriangleCount());
			geometry.setModelMatrix(g.getModelMatrix());
			
			geometry.setPos(g.getTranslation());
			geometry.setRotation(new float[]{g.getxRotation(), g.getyRotation(), g.getzRotation()});
			geometry.setScale(g.getScale());
			
			geometry.setColour(new Vector3f((float)r.nextDouble(), (float)r.nextDouble(), (float)r.nextDouble()));
			geometry.setInitialColour(geometry.getColour());
			
			geometry.create();
			
			if(geometry.getTriangleSize() == 0){
				System.err.println("Can't have 0 triangles");
				System.exit(1);
			}
			geometry.setIDs(shader.getVBOID(geometry.getVertices(), geometry.getIndices()));
			gameObjects.add(geometry);
			i++;
		}
		agent = (Box) gameObjects.get(1);
		desitination = gameObjects.get(2).getPos();
		direction = subtract(agent.getPos(), desitination);
		direction.normalise();
	}

	private void setUpShaders(){
		
		shader = new VCI_Shader("shaders/VCI.vert", "shaders/VCI.frag");
		
	}
		
	private void update(){
			
		camera.pollInput();
		if(Math.sqrt(subtract(agent.getPos(), desitination).lengthSquared()) > 1){
			agent.setPos(addVector(agent.getPos(), scaleVector(direction, agentMoveSpeed)));
		}
		
			
	}
	
	private Vector3f addVector(Vector3f a, Vector3f b){
		return new Vector3f(a.x + b.x, a.y + b.y, a.z + b.z);
	}
	
	private Vector3f subtract(Vector3f a, Vector3f b){
		return new Vector3f(b.x - a.x, b.y - a.y, b.z - a.z);
	}
	
	private Vector3f scaleVector(Vector3f a, float c){
		return new Vector3f(a.x * c, a.y * c, a.z * c);
	}
	
	private void draw(){
			
		GL20.glUseProgram(shader.getShaderID());
		GL11.glClear(GL11.GL_COLOR_BUFFER_BIT | GL11.GL_DEPTH_BUFFER_BIT);
		
		PVM.loadIdentity();
		camera.transform();
		
		for(Geometry g : gameObjects){
			g.doTransformation();
			shader.setColour(g.getColour());
			shader.render(g.getVboID(), g.getIndicesID(), g.getTriangleSize());
		}
		
		GL20.glUseProgram(0);
		Display.sync(60);
		Display.update();
			
	}
		
		
		
	public static void main(String[] args){
			(new Driver()).loop();
	}
}
