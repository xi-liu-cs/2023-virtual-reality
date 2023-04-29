import * as cg from "../render/core/cg.js";
import {g2} from "../util/g2.js";
import * as global from "../global.js";
import {Gltf2Node} from "../render/nodes/gltf2.js";
import { controllerMatrix, buttonState, joyStickState } from "../render/core/controllerInput.js";
import * as keyboardInput from "../util/input_keyboard.js";
let leftTriggerPrev = false;
let rightTriggerPrev = false;

let prev_time = 0;

let mouseX = 0, mouseY = 0;

document.addEventListener('mousemove', e => {
    let rect = canvas.getBoundingClientRect();
    mouseX = e.clientX - rect.left; /* canvas.offsetLeft */
    mouseY = e.clientY - rect.top; /* canvas.offsetTop */
});
 
export const init = async model => {
    let screen = model.add('cube').scale(7, 7, -7);
    let isHUD = false;
    model.control('h', 'toggle HUD', () => isHUD = ! isHUD);
    model.setRoom(false);
    model.setTable(false);
    model.identity().scale(8.0);
    let posZ = -100.;
    let posX = 0.;
    let posY = 0.;

    model.animate(() => {
        // Add Controllers
        let ml = controllerMatrix.left;
        let mr = controllerMatrix.right;
        let leftTrigger = buttonState.left[0].pressed;
        let rightTrigger = buttonState.right[0].pressed;

        posZ = -model.time;
        if (leftTrigger||rightTrigger){
            if (prev_time==0)
                prev_time = model.time;
            if (leftTrigger&&rightTrigger){
                posY +=1.;
                // posY += 10*model.time-prev_time;
            }
            else if (leftTrigger){
                posX -=1.;
                //posX -= 10*model.time-prev_time;
            }
            else {
                posX +=1.;
                //posX += 10*model.time-prev_time;
            }
        }
        else
            prev_time=0;


        let m = views[0]._viewMatrix, c = .5*Math.cos(model.time), s = .5*Math.sin(model.time);
        if (isHUD)
            model.hud();
        else
            model.setMatrix([m[0],m[4],m[8],0,m[1],m[5],m[9],0,m[2],m[6],m[10],0,0,1.6,-1,1]);
        model.scale(1.,1.,.0001);

        screen.flag('uRayTrace');
        /* model.setUniform('4fv','uL', [.5,.5,.5,1., -.5,-.5,-.5,.2, .7,-.7,0,.2, -.7,.7,0,.2]);
        model.setUniform('4fv','uS', [c,s,0,0, s,0,c,0, 0,c,s,0, -c,-s,0,0]);
        model.setUniform('4fv','uC', [0,0,0,2., 0,1,1,2, 1,0,1,2, 0,1,0,2]); */
        model.setUniform('4fv','uP', [posX,posY,posZ+0.1,0, 0,0,0,0, 0,0,0,0, 0,0,0,0]);

        /* clay.defineMesh('terrain', clay.createGrid(100, 100));
        let terrain = model.add('terrain').color(.1,.3,1).opacity(.7);

        terrain.flag('uTerrainTexture');
        terrain.identity().move(0, -0.7, 0).turnX(-.5 * Math.PI).scale(2);
        terrain.setVertices((u, v) => {
         return [(2 * u - 1) * model.time, (2 * v - 1) * model.time, .1 * cg.noise(30 * u - .3 * c, 30 * v * c, .3 * c)];
        }); */

        let canvas = document.getElementById('anidrawCanvas');
        model.setUniform('2fv', 'resolution', [canvas.width, canvas.height]);
        model.setUniform('2fv', 'mouse', [mouseX, mouseY]);

        model.customShader(`
#define ITR 60
#define FAR 400.

uniform int uRayTrace;
/* uniform vec4 uC[4], uL[4], uS[4]; */
uniform vec4 uP[4];
vec4 light[4], sphere[4];
const int nL = 1;
uniform int uObjTexture;
uniform int uSkyTexture;
uniform int uTerrainTexture;
uniform int uTerrain2Texture;
uniform vec3 uBgColor; /* background color */
uniform vec3 uLd[nL]; /* light direction */
uniform vec3 uLc[nL]; /* light color */
uniform float uLi; /* light intensity */
uniform vec4 uC[4], uL[4], uS[4];
uniform vec2 resolution,
mouse;
vec3 u_sky_color = vec3(.2, .1, .4);
vec3 u_floor_color = vec3(.7, .6, .5);
float fl = 3.;

float turbulence(vec3 p)
{
   float t = 0., f = 1.;
   for(int i = 0; i < 10; ++i)
   {
      t += abs(noise(f * p)) / f;
      f *= 2.;
   }
   return t;
}

float pattern(vec3 v)
{
   const int n = 10;
   float res = 0., f = 1.;
   for(int i = 1; i < n; ++i)
   {
      res += noise(f * v) / float(i);
      f *= float(i);
      f += float(i);
   }
   return res;
}

float ray_sphere(vec3 V, vec3 W, vec4 S)
{
   V -= S.xyz;
   float b = dot(V, W);
   float d = b * b - dot(V, V) + S.w * S.w;
   return d < 0. ? -1. : -b - sqrt(d);
}

vec3 shade_sphere(vec3 p, vec4 s, vec4 c)
{
   vec3 N = normalize(p - s.xyz);
   vec3 color = .1 * c.rgb;
   for(int l = 0; l < 4; l++)
   {
      vec3 lDir = light[l].xyz;
      float lBrightness = light[l].w;
      float t = -1.;
      for (int i = 0; i < 4; ++i)
         t = max(t, ray_sphere(p, lDir, sphere[i]));
      if(t < 0.)
      {
         vec3 R = 2. * N * dot(N, lDir) - lDir;
         color += lBrightness * (c.rgb * .9 * max(0., dot(N, lDir)) + c.a * vec3(pow(max(0., R.z), 10.)));
      }
   }
   return color;
}

float water_height(vec3 vertex)
{
   float n = noise(vertex * 8. + vec3(uTime, uTime, 0.));
   n += noise(vertex * 16. + vec3(uTime, -uTime, 0.));
   return n / 100.;
}

#define triangle_normal(p1, p2, p3)(cross(p2 - p1, p3 - p1))

vec3 water_normal(vec3 vertex)
{/* s1, s2, s2 are 3 points of a surface defined by y = f(x, z) */
   vec3 s1 = vec3(vertex.x - .01, vertex.y, vertex.z);
   vec3 s2 = vec3(vertex.x + .01, vertex.y, vertex.z);
   vec3 s3 = vec3(vertex.x, vertex.y, vertex.z - .01); /* vec3 s1 = vec3(vertex), s2 = vec3(vertex), s3 = vec3(vertex); */
   vec3 h1 = vec3(vertex.x - .01, water_height(s1), vertex.z);
   vec3 h2 = vec3(vertex.x + .01, water_height(s2), vertex.z);
   vec3 h3 = vec3(vertex.x, water_height(s3), vertex.z - .001);
   vec3 n = triangle_normal(h1, h2, h3);
   if(n.y < 0.) n = -n;
   return normalize(n);
}

vec3 water_surface(vec3 W, vec3 vertex)
{
   vec3 c = uBgColor * .5;
   vec3 N = water_normal(vertex);
   vec3 diffuse = vec3(.1);
   vec3 specular = vec3(.5);
   for (int l = 0 ; l < nL ; l++)
   {
      vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
      c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l])) * uLi
            + specular * pow(max(0., dot(R, W)), 20.)) * uLi;
   }
   vec3 rc;
   float rtMin = 10000.;
   if (rtMin < 10000.)
      c += .5 * rc;
   c *= uLi;
   return c;
}

mat2 mm2(float a){float c = cos(a), s = sin(a);return mat2(c,-s,s,c);}
mat2 m2 = mat2(0.934, 0.358, -0.358, 0.934);
float tri(float x){return abs(fract(x)-0.5);}

float heightmap(vec2 p)
{
   p *= .05;
   float z = 2.;
   float rz = 0.;
   for(int i = 1; i < 4; ++i)
   {
      rz += tri(p.x + tri(p.y * 1.5)) / z;
      z = z * -.85;
      p = p * 1.32;
      p *= m2;
   }
   rz += sin(p.y + cos(p.x*.9) + uTime) + cos(p.y + sin(p.x*.9) + uTime);
   rz *= .7;
   return rz*5.;
}

vec3 bary(vec2 a, vec2 b, vec2 c, vec2 p)
{
   vec2 v0 = b - a, v1 = c - a, v2 = p - a;
   float inv_denom = 1.0 / (v0.x * v1.y - v1.x * v0.y)+1e-9;
   float v = (v2.x * v1.y - v1.x * v2.y) * inv_denom;
   float w = (v0.x * v2.y - v2.x * v0.y) * inv_denom;
   float u = 1.0 - v - w;
   return abs(vec3(u,v,w));
}

float map(vec3 p)
{
   vec3 q = fract(p)-0.5;
   vec3 iq = floor(p);
   vec2 p1 = vec2(iq.x-.5, iq.z+.5);
   vec2 p2 = vec2(iq.x+.5, iq.z-.5);
   
   float d1 = heightmap(p1);
   float d2 = heightmap(p2);
   
   float sw = sign(q.x+q.z); 
   vec2 px = vec2(iq.x+.5*sw, iq.z+.5*sw);
   float dx = heightmap(px);
   vec3 bar = bary(vec2(.5*sw,.5*sw),vec2(-.5,.5),vec2(.5,-.5), q.xz);
   return (bar.x*dx + bar.y*d1 + bar.z*d2 + p.y + 3.)*.9;
}

float march(vec3 ro, vec3 rd)
{
   float precis = 0.001;
   float h=precis*2.0;
   float d = 0.;
   for( int i=0; i<ITR; i++ )
   {
      if( abs(h)<precis || d>FAR ) break;
      d += h;
      float res = map(ro+rd*d)*1.1;
      h = res;
   }
   return d;
}

#define DRAG_MULT 0.048
#define ITERATIONS_RAYMARCH 13
#define ITERATIONS_NORMAL 48
/* #define Mouse (iMouse.xy / iResolution.xy) */
#define Mouse (mouse.xy / resolution.xy)
#define Resolution (resolution.xy)
#define Time uTime

vec2 wavedx(vec2 position, vec2 direction, float speed, float frequency, float timeshift)
{
    float x = dot(direction, position) * frequency + timeshift * speed;
    float wave = exp(sin(x) - 1.0);
    float dx = wave * cos(x);
    return vec2(wave, -dx);
}

float getwaves(vec2 position, int iterations)
{
	float iter = 0.0;
    float phase = 6.0;
    float speed = 2.0;
    float weight = 1.0;
    float w = 0.0;
    float ws = 0.0;
    for(int i = 0; i < iterations; i++)
    {
        vec2 p = vec2(sin(iter), cos(iter));
        vec2 res = wavedx(position, p, speed, phase, Time);
        position += p * res.y * weight * DRAG_MULT;
        w += res.x * weight;
        iter += 12.0;
        ws += weight;
        weight = mix(weight, 0.0, 0.2);
        phase *= 1.18;
        speed *= 1.07;
    }
    return w / ws;
}

float raymarchwater(vec3 camera, vec3 start, vec3 end, float depth)
{
    vec3 pos = start;
    float h = 0.0;
    float hupper = depth;
    float hlower = 0.0;
    vec2 zer = vec2(0.0);
    vec3 dir = normalize(end - start);
    float eps = 0.01;
    for(int i=0;i<318;i++)
    {
        h = getwaves(pos.xz * 0.1, ITERATIONS_RAYMARCH) * depth - depth;
        float dist_pos = distance(pos, camera);
        if(h + eps*dist_pos > pos.y)
        {
            return dist_pos;
        }
        pos += dir * (pos.y - h);
        //eps *= 1.01;
    }
    return -1.0;
}

vec3 wave_normal(vec2 pos, float e, float depth)
{
    vec2 ex = vec2(e, 0);
    float H = getwaves(pos.xy * 0.1, ITERATIONS_NORMAL) * depth;
    vec3 a = vec3(pos.x, H, pos.y);
    return (cross(normalize(a-vec3(pos.x - e, getwaves((pos.xy - ex.xy)*0.1, ITERATIONS_NORMAL) * depth, pos.y)), 
                           normalize(a-vec3(pos.x, getwaves((pos.xy + ex.yx )* 0.1, ITERATIONS_NORMAL) * depth, pos.y + e))));
}
mat3 rotmat(vec3 axis, float angle)
{
	float s = sin(angle);
	float c = cos(angle);
	float oc = 1.0 - c;
	return mat3(oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s, 
	oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s, 
	oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}

vec3 getRay(vec2 uv)
{
    uv = (uv * 2.0 - 1.0) * vec2(Resolution.x / Resolution.y, 1.0);
	vec3 proj = normalize(vec3(uv.x, uv.y, 1.0) + vec3(uv.x, uv.y, -1.0) * pow(length(uv), 2.0) * 0.05);	
    if(Resolution.x < 400.0) return proj;
	vec3 ray = rotmat(vec3(0.0, -1.0, 0.0), 3.0 * (Mouse.x * 2.0 - 1.0)) * rotmat(vec3(1.0, 0.0, 0.0), 1.5 * (Mouse.y * 2.0 - 1.0)) * proj;
    return ray;
}

float intersectPlane(vec3 origin, vec3 direction, vec3 point, vec3 normal)
{ 
    return clamp(dot(point - origin, normal) / dot(direction, normal), -1.0, 9991999.0); 
}

vec3 extra_cheap_atmosphere(vec3 raydir, vec3 sundir)
{
	sundir.y = max(sundir.y, -0.07);
	float special_trick = 1.0 / (raydir.y * 1.0 + 0.1);
	float special_trick2 = 1.0 / (sundir.y * 11.0 + 1.0);
	float raysundt = pow(abs(dot(sundir, raydir)), 2.0);
	float sundt = pow(max(0.0, dot(sundir, raydir)), 8.0);
	float mymie = sundt * special_trick * 0.2;
	vec3 suncolor = mix(vec3(1.0), max(vec3(0.0), vec3(1.0) - vec3(5.5, 13.0, 22.4) / 22.4), special_trick2);
	vec3 bluesky= vec3(5.5, 13.0, 22.4) / 22.4 * suncolor;
	vec3 bluesky2 = max(vec3(0.0), bluesky - vec3(5.5, 13.0, 22.4) * 0.002 * (special_trick + -6.0 * sundir.y * sundir.y));
	bluesky2 *= special_trick * (0.24 + raysundt * 0.24);
	return bluesky2 * (1.0 + 1.0 * pow(1.0 - raydir.y, 3.0)) + mymie * suncolor;
}

vec3 getatm(vec3 ray)
{
 	return extra_cheap_atmosphere(ray, normalize(vec3(1.0))) * 0.5;
}

float sun(vec3 ray)
{
 	vec3 sd = normalize(vec3(1.0));   
    return pow(max(0.0, dot(ray, sd)), 528.0) * 110.0;
}
vec3 aces_tonemap(vec3 color)
{	
	mat3 m1 = mat3(
        0.59719, 0.07600, 0.02840,
        0.35458, 0.90834, 0.13383,
        0.04823, 0.01566, 0.83777
	);
	mat3 m2 = mat3(
        1.60475, -0.10208, -0.00327,
        -0.53108,  1.10813, -0.07276,
        -0.07367, -0.00605,  1.07602
	);
	vec3 v = m1 * color;    
	vec3 a = v * (v + 0.0245786) - 0.000090537;
	vec3 b = v * (0.983729 * v + 0.4329510) + 0.238081;
	return pow(clamp(m2 * (a / b), 0.0, 1.0), vec3(1.0 / 2.2));	
}

#define time uTime
#define USE_BMAP true
float PrSphDf(vec3 p, float r);
float PrCapsDf(vec3 p, float r, float h);
float SmoothBump(float lo, float hi, float w, float x);
vec2 Rot2D(vec2 q, float a);
float Hashfv3(vec3 p);
vec3 Hashv3f(float p);
vec3 VaryNf(vec3 p, vec3 n, float f);

vec4 vcId;
vec3 ltPos[2], boatPos[2];
float boatAng[2], dstFar, tCur, htWat, dstBMap;
int idObj;
bool uWat, hitWat;
const int idBoat = 1,
    idBLamp = 2,
    idFLamp = 3;
const float pi = 3.14159;

float ObjDf(vec3 p)
{
    vec3 q;
    float dMin, d;
    dMin = dstFar;
    for(int k = 0; k < 2; k++)
    {
        q = p - boatPos[k];
        q.xz = Rot2D(q.xz, boatAng[k]);
        d = max(PrCapsDf(q, 0.11, 0.25),
            -PrCapsDf(q + vec3(0., -0.02, 0.), 0.1, 0.24));
        if(d < dMin)
        {
            dMin = d;
            idObj = idBoat;
        }
        q.y -= 0.1;
        q.z -= 0.3;
        d = PrSphDf(q, 0.01);
        if(d < dMin)
        {
            dMin = d;
            idObj = idFLamp;
        }
        q.z -= -0.6;
        d = PrSphDf(q, 0.01);
        if(d < dMin)
        {
            dMin = d;
            idObj = idBLamp;
        }
    }
    return dMin;
}

float ObjRay(vec3 ro, vec3 rd)
{
    float dHit, d;
    dHit = 0.;
    for(int j = 0; j < 100; j++)
    {
        d = ObjDf(ro + dHit * rd);
        dHit += d;
        if(d < 0.001 || dHit > dstFar) break;
    }
    return dHit;
}

vec3 ObjNf(vec3 p)
{
    vec4 v;
    vec3 e = vec3(0.001, -0.001, 0.);
    v = vec4(ObjDf(p + e.xxx), ObjDf(p + e.xyy),
        ObjDf(p + e.yxy), ObjDf(p + e.yyx));
    return normalize(vec3(v.x - v.y - v.z - v.w) + 2. * v.yzw);
}

float VPoly(vec3 p)
{
    vec3 ip, fp, g, w, a;
    ip = floor(p);
    fp = fract(p);
    a = vec3(2.);
    for(float gz = -1.; gz <= 1.; gz++)
    {
        for(float gy = -1.; gy <= 1.; gy++)
        {
            for(float gx = -1.; gx <= 1.; gx++)
            {
                g = vec3(gx, gy, gz);
                w = g + 0.7 * Hashfv3(ip + g) - fp;
                a.x = dot(w, w);
                if(a.x < a.y)
                {
                    vcId = vec4(ip + g, a.y - a.x);
                    a = a.zxy;
                }
                else a.z = min(a.z, a.x);
            }
        }
    }
    return a.y;
}

vec3 TrackPath(float t)
{
    return vec3(4.7 * sin(t * 0.15) + 2.7 * cos(t * 0.19), 0., t);
}

float CaveDf(vec3 p)
{
    vec3 hv;
    float s, d;
    s = p.y - htWat;
    p.xy -= TrackPath(p.z).xy;
    p += 0.1 * (1. - cos(2. * pi * (p + 0.2 * (1. - cos(2. * pi * p.zxy)))));
    hv = cos(0.6 * p - 0.5 * sin(1.4 * p.zxy + 0.4 * cos(2.7 * p.yzx))) + 0.5 * sin(time) * noise(p); /* hv = cos(0.6 * p - 0.5 * sin(1.4 * p.zxy + 0.4 * cos(2.7 * p.yzx))); */
    if(USE_BMAP && dstBMap < 10.) hv *= 1. + 0.01 *
        (1. - smoothstep(0., 10., dstBMap)) *
        smoothstep(0.05, 0.4, VPoly(10. * p)) / length(hv);
    d = 0.9 * (length(hv) - 1.1);
    if(!uWat) d = min(d, s);
    return d;
}

float CaveRay(vec3 ro, vec3 rd)
{
    float d, dHit;
    dHit = 0.;
    for(int j = 0; j < 200; j++)
    {
        dstBMap = dHit;
        d = CaveDf(ro + dHit * rd);
        dHit += d;
        if(d < 0.001 || dHit > dstFar) break;
    }
    return dHit;
}

vec3 CaveNf(vec3 p)
{
    vec4 v;
    const vec3 e = vec3(0.001, -0.001, 0.);
    v = vec4(CaveDf(p + e.xxx), CaveDf(p + e.xyy),
        CaveDf(p + e.yxy), CaveDf(p + e.yyx));
    return normalize(vec3(v.x - v.y - v.z - v.w) + 2. * v.yzw);
}

float CaveSShadow(vec3 ro, vec3 rd)
{
    float sh, d, h;
    sh = 1.;
    d = 0.1;
    for(int j = 0; j < 16; j++)
    {
        h = CaveDf(ro + rd * d);
        sh = min(sh, smoothstep(0., 0.05 * d, h));
        d += max(0.2, 0.1 * d);
        if(sh < 0.05) break;
    }
    return 0.4 + 0.6 * sh;
}

vec3 CaveCol(vec3 ro, vec3 rd, vec3 ltDir, float atten)
{
    vec3 col, vn, q, vno;
    float glit;
    VPoly(10. * ro);
    q = ro;
    if(!USE_BMAP) q = 0.004 * floor(250. * q);
    vn = VaryNf(10. * q, CaveNf(q), 1.); /* vn = sin(time) * VaryNf(10. * q, CaveNf(q), 1.); */
    col = (vec3(0.2, 0.1, 0.1) + vec3(0.3, 0.2, 0.1) * Hashv3f(Hashfv3(vcId.xyz))) *
        (1.2 - 0.4 * Hashfv3(100. * ro)) *
        (0.4 + 0.6 * smoothstep(0.05, 1., sqrt(vcId.w))) *
        (0.2 + 0.8 * max(dot(vn, ltDir), 0.) +
            2. * pow(max(dot(normalize(ltDir - rd), vn), 0.), 256.));
    if(!hitWat)
    {
        vno = CaveNf(ro);
        glit = 20. * pow(max(0., dot(ltDir, reflect(rd, vno))), 4.) *
            pow(1. - 0.6 * abs(dot(normalize(ltDir - rd),
                VaryNf(100. * ro, vno, 5.))), 8.);
        col += vec3(1., 1., 0.5) * glit;
    }
    col *= atten * CaveSShadow(ro, ltDir);
    return col;
}

vec3 ObjCol(vec3 ro, vec3 rd, vec3 vn, vec3 ltDir, float atten)
{
    vec4 col4;
    if(idObj == idBoat) col4 = vec4(0.3, 0.3, 0.6, 0.2);
    else if(idObj == idFLamp) col4 = vec4(0., 1., 0., -1.);
    else if(idObj == idBLamp) col4 = vec4(1., 0., 0., -1.);
    if(col4.a >= 0.)
        col4.rgb = col4.rgb * (0.2 + 0.8 * CaveSShadow(ro, ltDir)) *
        (0.1 + 0.9 * atten * max(dot(ltDir, vn), 0.)) +
        col4.a * atten * pow(max(dot(normalize(ltDir - rd), vn), 0.), 64.);
    return col4.rgb;
}

vec3 ShowScene(vec3 ro, vec3 rd)
{
    vec3 col, colR, bgCol, ltVec, vn, roo, rdo, row, vnw;
    float dstCave, dstObj, atten, frFac;
    roo = ro;
    rdo = rd;
    bgCol = (abs(rd.y) < 0.5) ? vec3(0., 0.05, 0.08) : vec3(0.01);
    uWat = false;
    hitWat = false;
    dstCave = CaveRay(ro, rd);
    dstObj = ObjRay(ro, rd);
    if(dstCave < min(dstObj, dstFar) && ro.y + rd.y * dstCave < htWat + 0.001)
    {
        hitWat = true;
        ro += rd * dstCave;
        row = ro;
        vnw = VaryNf(1.5 * ro, vec3(0., 1., 0.), 0.1);
        rd = reflect(rd, vnw);
        ro += 0.01 * rd;
        dstCave = CaveRay(ro, rd);
        dstObj = ObjRay(ro, rd);
    }
    if(min(dstCave, dstObj) < dstFar)
    {
        ltVec = roo + ltPos[0] - ro;
        atten = 1. / (0.1 + dot(ltVec, ltVec));
        if(hitWat) atten *= 3.;
        ltVec = normalize(ltVec);
        ro += min(dstCave, dstObj) * rd;
        if(dstCave < dstObj) col = mix(CaveCol(ro, rd, ltVec, atten), bgCol,
            smoothstep(0.45, 0.99, dstCave / dstFar));
        else col = ObjCol(ro, rd, ObjNf(ro), ltVec, atten);
    }
    else col = bgCol;
    if(hitWat) /* hit water */
    {
        frFac = rdo.y * rdo.y;
        frFac *= frFac;
        if(frFac > 0.005)
        {
            vec3 color_offset = vec3(clamp(-cos(time), -0.2, 0.2), 0., clamp(-sin(time), -0.2, 0.2));
            vec2 uv = gl_FragCoord.xy / resolution.xy;
            float waterdepth = 2.1;
            vec3 wfloor = vec3(0.0, -waterdepth, 0.0);
            vec3 wceil = vec3(0.0, 0.0, 0.0);
            vec3 orig = vec3(0.0, 2.0, 0.0);
            vec3 ray = getRay(uv);
            float hihit = intersectPlane(orig, ray, wceil, vec3(0.0, 1.0, 0.0));
            if(ray.y >= -0.01)
            {
                vec3 C = getatm(ray) * 2.0 + sun(ray);
                //tonemapping
                C = aces_tonemap(C) + color_offset;
                bgCol = C;
            }
            else
            {
                float lohit = intersectPlane(orig, ray, wfloor, vec3(0.0, 1.0, 0.0));
                vec3 hipos = orig + ray * hihit;
                vec3 lopos = orig + ray * lohit;
                float dist = raymarchwater(orig, hipos, lopos, waterdepth);
                vec3 pos = orig + ray * dist;

                vec3 N = wave_normal(pos.xz, 0.001, waterdepth);
                vec2 velocity = N.xz * (1.0 - N.y);
                N = mix(vec3(0.0, 1.0, 0.0), N, 1.0 / (dist * dist * 0.01 + 1.0));
                vec3 R = reflect(ray, N);
                float fresnel = (0.04 + (1.0-0.04)*(pow(1.0 - max(0.0, dot(-N, ray)), 5.0)));
                
                vec3 C = fresnel * getatm(R) * 2.0 + fresnel * sun(R);
                //tonemapping
                C = aces_tonemap(C) + color_offset;
                bgCol = C;
            }
            rd = refract(rdo, vnw, 1. / 1.333);
            ro = row + 0.01 * rd;
            uWat = true;
            dstCave = CaveRay(ro, rd);
            if(min(dstCave, dstObj) < dstFar)
            {
                ltVec = roo + ltPos[1] - ro;
                atten = 1. / (0.1 + dot(ltVec, ltVec));
                ltVec = normalize(ltVec);
                ro += rd * dstCave;
                hitWat = false;
                colR = mix(CaveCol(ro, rd, ltVec, atten), bgCol,
                    smoothstep(0.45, 0.99, dstCave / dstFar));
            }
            else colR = bgCol;
            /* col = mix(col, colR * vec3(0.4, 1., 0.6) * exp(0.02 * ro.y), frFac); */
            col = mix(col, colR * vec3(0.4, 1., 0.6) * exp(0.02 * ro.y), frFac);
        }
    }
    return pow(clamp(col, 0., 1.), vec3(0.8));
}

float PrSphDf(vec3 p, float r)
{
    return length(p) - r;
}

float PrCapsDf(vec3 p, float r, float h)
{
    return length(p - vec3(0., 0., clamp(p.z, -h, h))) - r;
}

float SmoothBump(float lo, float hi, float w, float x)
{
    return (1. - smoothstep(hi - w, hi + w, x)) * smoothstep(lo - w, lo + w, x);
}

vec2 Rot2D(vec2 q, float a)
{
    return q * cos(a) + q.yx * sin(a) * vec2(-1., 1.);
}

const vec4 cHashA4 = vec4(0., 1., 57., 58.);
const vec3 cHashA3 = vec3(1., 57., 113.);
const float cHashM = 43758.54;

float Hashfv2(vec2 p)
{
    return fract(sin(dot(p, cHashA3.xy)) * cHashM);
}

float Hashfv3(vec3 p)
{
    return fract(sin(dot(p, cHashA3)) * cHashM);
}

vec3 Hashv3f(float p)
{
    return fract(sin(vec3(p, p + 1., p + 2.)) *
        vec3(cHashM, cHashM * 0.43, cHashM * 0.37));
}

vec4 Hashv4f(float p)
{
    return fract(sin(p + cHashA4) * cHashM);
}

float Noisefv2(vec2 p)
{
    vec4 t;
    vec2 ip, fp;
    ip = floor(p);
    fp = fract(p);
    fp = fp * fp * (3. - 2. * fp);
    t = Hashv4f(dot(ip, cHashA3.xy));
    return mix(mix(t.x, t.y, fp.x), mix(t.z, t.w, fp.x), fp.y);
}

float Fbmn(vec3 p, vec3 n)
{
    vec3 s;
    float a;
    s = vec3(0.);
    a = 1.;
    for(int i = 0; i < 5; i++)
    {
        s += a * vec3(Noisefv2(p.yz), Noisefv2(p.zx), Noisefv2(p.xy));
        a *= 0.5;
        p *= 2.;
    }
    return dot(s, abs(n));
}

vec3 VaryNf(vec3 p, vec3 n, float f)
{
    vec3 g;
    const vec3 e = vec3(0.1, 0., 0.);
    g = vec3(Fbmn(p + e.xyy, n), Fbmn(p + e.yxy, n), Fbmn(p + e.yyx, n)) -
        Fbmn(p, n);
    return normalize(n + f * (g - n * dot(n, g)));
}

---
if (uRayTrace == 1)
{
   float fl = -1. / uProj[3].z; // FOCAL LENGTH OF VIRTUAL CAMERA
   
   /* for (int i = 0 ; i < 4 ; i++)
   {
      light[i]  = vec4((uView * vec4(uL[i].xyz,0.)).xyz,uL[i].w);
      sphere[i] = vec4((uView * uS[i]).xyz,.25) - vec4(0.,0.,fl,0.);
   } */
   
   /* vec3 V = vec3(0.,10.,uP[0].z);
   vec3 W = normalize(vec3(2.*vUV.x-1.,1.-2.*vUV.y,-fl));
   //float tMin = 1000.;
   
   float rz = march(V,W);
   if ( rz < FAR )
   {
      vec3 pos = V + rz * W;
      vec2 e = vec2(-1., 1.)*0.01;
      vec3 nor = normalize(e.yxx*map(pos + e.yxx) + e.xxy*map(pos + e.xxy) + 
                  e.xyx*map(pos + e.xyx) + e.yyy*map(pos + e.yyy)) + water_normal(pos);
      vec3 ligt = normalize(vec3(-.2, 0.05, -0.2));
      float dif = clamp(dot( nor, ligt ), 0., 1.);
      float fre = pow(clamp(1.0+dot(nor,W),0.0,1.0), 3.);

      vec3 diffuse = vec3(.1);
      vec3 specular = vec3(.5);
      for (int l = 0 ; l < nL ; l++)
      {
         vec3 R = 2. * dot(nor, uLd[l]) * nor - uLd[l];
         color += uLc[l] * (diffuse * max(0.,dot(nor, uLd[l])) * uLi
               + specular * pow(max(0., dot(R, W)), 20.)) * uLi;
      }
      vec3 brdf = 2.*vec3(0.10,0.11,0.1);
      brdf += 1.9*dif*vec3(.8,1.,.05);
      color = color*brdf + fre*0.1*vec3(.7,.8,1.);
      color += .4 * turbulence(vAPos);
      color *= vec3(0.32) + vec3(-0.3, -0.3, 0.4);
   }
   else
   {
      float starXf = (atan(W.x, W.z) + 1.57) / 6.28;
      float starYf = (asin(W.y) / 1.57);
      int starX = int(starXf * 1000.0 * 16.0);
      int starY = int(starYf * 250.0 * 16.0);
      float starTest = float(7 + starX * starY * 13);
      float value = abs(mod(starTest, 5000.0));
      if( value >= 0.0 && value <= .5)
      {
         color = vec3(value * 0.5 + .5);
      }
      else
      {
         opacity = 0.2;
      }
      color = noise(vAPos) * u_sky_color;
   } */

   /* vec3 color_offset = vec3(clamp(-cos(Time), -0.2, 0.2), 0., clamp(-sin(Time), -0.2, 0.2));
   vec2 uv = gl_FragCoord.xy / resolution.xy;
   float waterdepth = 2.1;
   vec3 wfloor = vec3(0.0, -waterdepth, 0.0);
   vec3 wceil = vec3(0.0, 0.0, 0.0);
   vec3 orig = vec3(0.0, 2.0, 0.0);
   vec3 ray = getRay(uv);
   float hihit = intersectPlane(orig, ray, wceil, vec3(0.0, 1.0, 0.0));
   if(ray.y >= -0.01){
       vec3 C = getatm(ray) * 2.0 + sun(ray);
       //tonemapping
       color = aces_tonemap(C) + color_offset;   
       return;
   }
   float lohit = intersectPlane(orig, ray, wfloor, vec3(0.0, 1.0, 0.0));
   vec3 hipos = orig + ray * hihit;
   vec3 lopos = orig + ray * lohit;
   float dist = raymarchwater(orig, hipos, lopos, waterdepth);
   vec3 pos = orig + ray * dist;

   vec3 N = wave_normal(pos.xz, 0.001, waterdepth);
   vec2 velocity = N.xz * (1.0 - N.y);
   N = mix(vec3(0.0, 1.0, 0.0), N, 1.0 / (dist * dist * 0.01 + 1.0));
   vec3 R = reflect(ray, N);
   float fresnel = (0.04 + (1.0-0.04)*(pow(1.0 - max(0.0, dot(-N, ray)), 5.0)));
   
   vec3 C = fresnel * getatm(R) * 2.0 + fresnel * sun(R);
   //tonemapping
   color = aces_tonemap(C) + color_offset; */

    mat3 vuMat;
    vec2 mPtr;
    vec3 ro, rd, fpF, fpB, vd;
    vec2 canvas, uv, ori, ca, sa;
    float el, az, t, tt, a;
    canvas = resolution.xy;
    uv = 2. * gl_FragCoord.xy / canvas - 1.;
    uv.x *= canvas.x / canvas.y;
    tCur = time;
    mPtr = mouse;
    mPtr.xy = mPtr.xy / canvas - 0.5;
    t = 1. * tCur;
    az = 0.;
    el = 0.;
    az = az + 2. * pi * mPtr.x;
    el = el + 0.95 * pi * mPtr.y;
    /* if(mPtr.z > 0.)
    {
        az = az + 2. * pi * mPtr.x;
        el = el + 0.95 * pi * mPtr.y;
    }
    else
    {
        tt = mod(floor(0.05 * tCur), 4.);
        a = 0.45 * pi * SmoothBump(0.75, 0.95, 0.05, mod(0.05 * tCur, 1.));
        if(tt < 2.) el = (2. * tt - 1.) * a;
        else az = (2. * tt - 5.) * a;
    } */
    htWat = -0.5;
    for(int k = 0; k < 2; k++)
    {
        fpF = TrackPath(t + 3. + 3. * float(k) + 0.1);
        fpB = TrackPath(t + 3. + 3. * float(k) - 0.1);
        boatPos[k] = 0.5 * (fpF + fpB);
        boatPos[k].y = htWat + 0.01;
        vd = fpF - fpB;
        boatAng[k] = (length(vd.xz) > 0.) ? atan(vd.x, vd.z) : 0.5 * pi;
    }
    fpF = TrackPath(t + 0.1);
    fpB = TrackPath(t - 0.1);
    ro = 0.5 * (fpF + fpB);
    vd = fpF - fpB;
    ori = vec2(el, az + ((length(vd.xz) > 0.) ? atan(vd.x, vd.z) : 0.5 * pi));
    ca = cos(ori);
    sa = sin(ori);
    vuMat = mat3(ca.y, 0., -sa.y, 0., 1., 0., sa.y, 0., ca.y) *
        mat3(1., 0., 0., 0., ca.x, -sa.x, 0., sa.x, ca.x);
    rd = vuMat * normalize(vec3(uv, 2.));
    ltPos[0] = 0.5 * vuMat * vec3(0., 1., -1.);
    ltPos[1] = 0.5 * vuMat * vec3(0., -1., -1.);
    dstFar = 50.;
    color = ShowScene(ro, rd);
}

if(uTerrainTexture == 1)
{
  color = uBgColor * uLi;
  vec3 V = vec3(0.,.2,fl);
  vec3 W = normalize(vec3(vPos.xy, -fl));
  float tMin = 10000.;
  float tS = -V.y / W.y;
  color = water_surface(W, V + tS * W);
  color += .1 * turbulence(5. * vAPos);
}
      `);
    });
}