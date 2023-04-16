import * as cg from "../render/core/cg.js";
import {g2} from "../util/g2.js";
import * as global from "../global.js";
import {Gltf2Node} from "../render/nodes/gltf2.js";
import { controllerMatrix, buttonState, joyStickState } from "../render/core/controllerInput.js";
import * as keyboardInput from "../util/input_keyboard.js";
let leftTriggerPrev = false;
let rightTriggerPrev = false;

let prev_time = 0;

export const init = async model => {
    let screen = model.add('cube').scale(5, 5, -5);
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

---
if (uRayTrace == 1)
{
   float fl = -1. / uProj[3].z; // FOCAL LENGTH OF VIRTUAL CAMERA
   
   /* for (int i = 0 ; i < 4 ; i++)
   {
      light[i]  = vec4((uView * vec4(uL[i].xyz,0.)).xyz,uL[i].w);
      sphere[i] = vec4((uView * uS[i]).xyz,.25) - vec4(0.,0.,fl,0.);
   } */
   
   vec3 V = vec3(0.,10.,uP[0].z);
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
      /* float starXf = (atan(W.x, W.z) + 1.57) / 6.28;
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
      } */
      color = noise(vAPos) * u_sky_color;
   }
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