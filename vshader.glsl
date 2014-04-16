attribute vec3 aVertexPosition;
attribute vec2 aTextureCoordinates;

//uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

varying vec2 vTextureCoordinates;
varying vec3 position;

void main() {
  vTextureCoordinates = aTextureCoordinates;
  position = aVertexPosition;
  gl_Position = projectionMatrix * /*modelViewMatrix */ vec4(aVertexPosition, 1.0);
}