#version 330

uniform sampler2D sample;

out vec4 vFragColor;
in vec2 vTextureO;



void main(void){ 

   vFragColor = texture(sample, vTextureO);
}