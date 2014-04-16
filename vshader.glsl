attribute vec3 aVertexPosition;
attribute vec2 aTextureCoordinates;

//uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

varying vec2 vTextureCoordinates;
varying vec2 position;

void main() {
  vTextureCoordinates = aTextureCoordinates;
  gl_Position = projectionMatrix * /*modelViewMatrix */ vec4(aVertexPosition, 1.0);
  position = vec2(aVertexPosition.x, aVertexPosition.y);
}