#version 130

uniform vec2 uTopLeft;
uniform vec2 uBottomRight;
uniform vec2 uTLPos;
uniform vec2 uBRPos;

smooth out vec2 fragPos;

void main()
{
    vec4 pos = gl_ModelViewProjectionMatrix * gl_Vertex;
    gl_Position = pos;
    gl_TexCoord[0] = gl_TextureMatrix[0] * gl_MultiTexCoord0;
    gl_FrontColor = gl_Color;
    vec2 tl = gl_Vertex.xy - uTLPos;
    vec2 br = gl_Vertex.xy - uBRPos;
    fragPos.x = abs(tl.x) < abs(br.x) ? uTopLeft.x : uBottomRight.x;
    fragPos.y = abs(tl.y) < abs(br.y) ? uTopLeft.y : uBottomRight.y;
}
