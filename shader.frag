#version 130

#define FLT_MAX 10000000000.//3.402823466e+38

uniform sampler2D texture;

vec4 blur(vec2 texCoord, vec2 d)
{
    vec4 colour = vec4(0.f, 0.f, 0.f, 1.f);
    mat3 kernel = mat3(1,2,1,2,4,2,1,2,1)/16;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            colour += kernel[i][j] * texture2D(texture, texCoord + vec2(d.x * (i-1), d.y * (j-1)));
        }
    }
    colour.w = 1.f;
    return colour;
}

vec4 hsv(float h, float s, float v)
{
    float c = s * v;
    h *= 0.954929659;
    float x = c * (1. - abs(mod(h, 2.) - 1.));
    vec4 col = vec4(0, 0, 0, 1);
    if (h <= 1.) col.rgb = vec3(c, x, 0);
    else if (h <= 2.) col.rgb = vec3(x, c, 0);
    else if (h <= 3.) col.rgb = vec3(0, c, x);
    else if (h <= 4.) col.rgb = vec3(0, x, c);
    else if (h <= 5.) col.rgb = vec3(x, 0, c);
    else if (h <= 6.) col.rgb = vec3(c, 0, x);
    float m = v - c;
    col.rgb += vec3(m, m, m);
    return col;
}

vec2 cmul(vec2 a, vec2 b)
{
    return vec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

vec2 cpow(vec2 z, int exponent)
{
    vec2 ret = vec2(1, 0);
    while (exponent > 0)
    {
        if ((exponent & 1) == 1) ret = cmul(ret, z);
        z = cmul(z, z);
        exponent >>= 1;
    }
    return ret;
}

vec2 conj(vec2 z)
{
    return vec2(z.x, -z.y);
}

float cnorm(vec2 z)
{
    return z.x*z.x+z.y*z.y;
}

vec2 cdiv(vec2 a, vec2 b)
{
    float den = 1 / (b.x*b.x + b.y*b.y);
    return vec2(a.x * b.x + a.y * b.y, a.y * b.x - a.x * b.y) * den;
}

float atan2(float y, float x)
{
    if (x > 0) return atan(y/x);
    if (x < 0)
    {
        if (y >= 0) return atan(y/x)+3.14159265;
        else return atan(y/x)-3.14159265;
    }
    if (y < 0) return -1.570796327;
    return 1.570796327;
}

float carg(vec2 z)
{
    return atan2(z.y, z.x);
}

vec3 newton(vec2 chaos, vec2 z, int exponent, int maxIt, float maxDist)
{
    vec2 one = vec2(1, 0);
    vec2 p = cpow(z, exponent);
    vec2 pm = cpow(z, exponent-1);
    int it = 0;
    float pmn = cnorm(p-one);
    while (pmn > maxDist && ++it < maxIt
           && p.x <= FLT_MAX
           && p.y <= FLT_MAX)
    {
        vec2 den = exponent * pm;
        z = z - cdiv(cmul(chaos, p-one), den);
        p = cpow(z, exponent);
        pm = cpow(z, exponent-1);
        pmn = cnorm(p-one);
    }
    return vec3(z.x, z.y, it);
}


vec4 newtonCol(vec2 chaos, vec2 z, int exponent, int maxIt, float maxDist)
{
    vec3 pack = newton(chaos, z, exponent, maxIt, maxDist);
    z = pack.xy;
    int it = int(pack.z);
    float arg = carg(z);
    if (arg > FLT_MAX) return vec4(1, 1, 1, 1);
    return hsv(carg(z)+3.14159265, 1., 1. - 0.05 * (it / maxIt));
}

vec4 hsv(vec3 _hsv)
{
    return hsv(_hsv.x, _hsv.y, _hsv.z);
}

smooth in vec2 fragPos;

uniform vec2 uChaos;
uniform int uExponent;
uniform int uMaxIt;
uniform float uMaxDist;

uniform float lerp;

void main()
{
    vec4 expected = texture2D(texture, gl_TexCoord[0].xy);
    vec4 calculated = newtonCol(uChaos, fragPos, uExponent, uMaxIt, uMaxDist);
    gl_FragColor = mix(expected, calculated, clamp(lerp, 0, 1));
}
