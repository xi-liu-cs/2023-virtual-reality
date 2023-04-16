import * as cg from "../render/core/cg.js";
import {g2} from "../util/g2.js";
import * as global from "../global.js";
import {Gltf2Node} from "../render/nodes/gltf2.js";
import { controllerMatrix, buttonState, joyStickState } from "../render/core/controllerInput.js";
import * as keyboardInput from "../util/input_keyboard.js";

export const init = async model =>
{
    let gltf = new Gltf2Node({url: './media/gltf/0/mountain.glb'});
    gltf.translation = [20, 0, -20];
    global.gltfRoot.addNode(gltf);

    let isAnimate = true, isBlending = true, isRubber = true, t = 0;
    model.setTable(false);
    model.setRoom(false);

    model.move(0, -0.5, 0);
    model.control('a', 'animate' , () => isAnimate  = ! isAnimate );
    model.control('b', 'blending', () => isBlending = ! isBlending);
    model.control('r', 'rubber'  , () => isRubber   = ! isRubber  );
 
    model.color(1,.5,.5);
 
    /* model.move(0,1.5,0).scale(.3);
    let shape1 = model.add('sphere');
    let shape2 = model.add('sphere');
    let elbow  = model.add('sphere');
    let muscle = model.add('sphere'); */

   let dim = 20,
   one_over_dim = 1 / dim,
   N = dim * dim * dim; /* let N = 10000 */
   let particles = model.add('particles').info(N).texture('media/textures/disk.jpg').flag('uTransparentTexture');

   let m = views[0]._viewMatrix, c = .5 * Math.cos(model.time), s = .5 * Math.sin(model.time);

   clay.defineMesh('terrain', clay.createGrid(100, 100));
   let terrain = model.add('terrain').color(.1,.3,1).opacity(.7); /* water */

   clay.defineMesh('terrain0', clay.createGrid(100, 100));
   let terrain0 = model.add('terrain0').color(.1,.3,1).opacity(.7); /* water */

   clay.defineMesh('terrain1', clay.createGrid(100, 100));
   let terrain1 = model.add('terrain1').color(.1,.3,1).opacity(.7); /* water */

   clay.defineMesh('terrain2', clay.createGrid(100, 100));
   let terrain2 = model.add('terrain2').opacity(1); /* island */

   clay.defineMesh('terrain3', clay.createGrid(100, 100));
   let terrain3 = model.add('terrain3').color(.1,.3,1).opacity(.7); /* water */

   clay.defineMesh('terrain4', clay.createGrid(100, 100));
   let terrain4 = model.add('terrain4').color(.1,.3,1).opacity(.7); /* water */

   clay.defineMesh('terrain5', clay.createGrid(100, 100));
   let terrain5 = model.add('terrain5').color(.1,.3,1).opacity(.7); /* water */

   model.customShader(`
   const int nL = 1;
   uniform int uObjTexture;
   uniform int uSkyTexture;
   uniform int uTerrainTexture;
   uniform int uTerrain2Texture;
   uniform int uRayTrace;
   uniform vec3 uBgColor; /* background color */
   uniform vec3 uLd[nL]; /* light direction */
   uniform vec3 uLc[nL]; /* light color */
   uniform float uLi; /* light intensity */
   uniform vec4 uC[4], uL[4], uS[4];
   vec4 light[4], sphere[4];
   vec3 u_sky_color = vec3(.2, .2, .7);
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

   ---

   if(uObjTexture == 1)
   {
      float n = 0.;
      for(int i = 0; i < 3; ++i)
         n += noise(vPos + vec3(i, i, i));
      color += .1 * vec3(1.0 - sin(uTime * n), 1.0 - cos(uTime * n), sin(uTime * 1.3 * n));
   }

   if(uSkyTexture == 1)
   {
      float height = 1.0 / vAPos.y;
      if(vAPos.y < 0.1)
         color = u_floor_color + .1 * pattern(3. * vAPos);
      else
         color = 0.5 * height * u_sky_color;
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

   if(uTerrain2Texture == 1)
   {
      color *= 1. + pattern(3. * vAPos);
   }

   if(uRayTrace == 1)
   {
      float fl = -1. / uProj[3].z; // FOCAL LENGTH OF VIRTUAL CAMERA
      for(int i = 0; i < 4; ++i)
      {
         light[i]  = vec4((uView * vec4(uL[i].xyz,0.)).xyz,uL[i].w);
         sphere[i] = vec4((uView * uS[i]).xyz,.25) - vec4(0.,0.,fl,0.);
      }
      vec3 V = vec3(0.);
      vec3 W = normalize(vec3(2.*vUV.x-1.,1.-2.*vUV.y,-fl));
      float tMin = 1000.;
      for(int i = 0; i < 4; ++i)
      {
         float t = ray_sphere(V, W, sphere[i]);
         if(t > 0. && t < tMin)
         {
            tMin = t;
            color = shade_sphere(t * W, sphere[i], uC[i]) + pattern(vAPos);
         }
      }
      if(tMin == 1000.)
         opacity = 0.;
   }
   `);
   let roomBackground = model.add('roomBackground');
   roomBackground.flag('uSkyTexture');
   /* model.add('sphere').scale(100,100,-100).flag('uSkyTexture'); */

   let data = [], V = [];

   function hash(i, j, k)
   {
        return (i * dim + j) * dim + k;
   }

   function rand(min, max)
   {
        return Math.random() * (max - min) + min;
   }

   for(let i = 0; i < dim; ++i)
   {
        for(let j = 0; j < dim; ++j)
        {
            for(let k = 0; k < dim; ++k)
            {
                //  p: position of the particle
                //  n: surface normal of the particle (defaults to facing the camera)
                //  s: size of the particle (defaults to 0.01)
                //  c: color of the particle (defaults to white)
                //  t: texture range [ uLo,vLo,uHi-uLo,vHi-vLo ] (defaults to [0,0,1,1])
                let r = rand(0, 5),
                x = i * one_over_dim, y = j * one_over_dim, z = k * one_over_dim,
                factor = 5;
                data.push({
                    // p: [ Math.random() - .5, Math.random() - .5, Math.random() - .5 ],
                    p: [factor * r * x, 0.3 * y, factor * r * z],
                    s: 0.1 * Math.random(),
                    c: [Math.random(), Math.random(), Math.random()],
                });
                V[hash(i, j, k)] = [0,0,0];
            }
        }
   }

   function uv_function(u, v)
   {
    return [2*u-1,2*v-1,.4 * cg.noise(3*u-model.time,3*v,model.time)];
   }

    /* for(let n = 0; n < N; ++n)
    {
        //  p: position of the particle
        //  n: surface normal of the particle (defaults to facing the camera)
        //  s: size of the particle (defaults to 0.01)
        //  c: color of the particle (defaults to white)
        //  t: texture range [ uLo,vLo,uHi-uLo,vHi-vLo ] (defaults to [0,0,1,1])

        let id = n / N;

        data.push({
        p: uv_function(id, id),
        s: .008 + .08 * Math.random(),
        c: [Math.random(), Math.random(), Math.random()],
        });
        V[n] = [0,0,0];
    } */

   console.log(data);

   /* let nu = 60, nv = 60;
   clay.defineMesh('terrain', clay.createGrid(nu, nv));
   let terrain = model.add('terrain').color(0,.5,1).opacity(.7); */

   model.move(0,1.5,0).scale(.6).animate(() => {
    /* model.blend(isBlending);
    model.melt(isAnimate && ! isRubber);
    t += isAnimate ? model.deltaTime : 0;
    shape1.identity().move(-1,0,0).scale(1.1,.28,.28);
    let bend = .7-.7*Math.cos(t);
    shape2.identity().turnZ(bend).move(.8,0,0).scale(1,.23,.23);
    elbow .identity().move(-.07,-.07,0).scale(.2);
    muscle.identity().move(-.8,.05+.15*bend,0).scale(.35,.2,.24); */

    model.setUniform('4fv','uL', [.5,.5,.5,1., -.5,-.5,-.5,.2, .7,-.7,0,.2, -.7,.7,0,.2]);
    model.setUniform('4fv','uS', [c,s,0,0, s,0,c,0, 0,c,s,0, -c,-s,0,0]);
    model.setUniform('4fv','uC', [c,s,0,2, s,0,c,2, 0,c,s,2, c,0,s,2]);
    let timeOfDay = 98;
    let ldX = 1. - timeOfDay / 50.;
    let ldY = Math.sin((timeOfDay / 100.) * 3.14159265);
    let lcG = .5 + .5 * Math.sin((timeOfDay / 100.) * 3.14159265);
    let lcB = Math.sin((timeOfDay / 100.) * 3.14159265);
    let ldData = [cg.normalize([.2,ldY,ldX])];
    let lIntensity = 1. + .5 * Math.sin((timeOfDay / 100.) * 3.14159265);
    let lcData = [1,lcG,lcB];
    model.setUniform('3fv', 'uLd', ldData.flat());
    model.setUniform('3fv', 'uLc', lcData);
    model.setUniform('1f', 'uLi', lIntensity);
    model.setUniform('3fv', 'uBgColor', [ .15,.2,.85 ]);

    particles.flag('uTerrainTexture');

    terrain.flag('uTerrainTexture');
    terrain.identity().move(0, 0.5, 0).turnX(-.5 * Math.PI).scale(2);
    terrain.setVertices((u, v) => {
       return [(2 * u - 1) * model.time, (2 * v - 1) * model.time, .1 * cg.noise(30 * u - .3 * c, 30 * v * c, .3 * c)];
    });
 
    terrain0.flag('uTerrainTexture');
    terrain0.identity().move(0, 0.5, 0).turnX(-.5 * Math.PI).scale(2);
    terrain0.setVertices((u, v) => {
       return [2 * u - 1, 2 * v - 1, .2 * cg.noise(3 * u - model.time, 3 * v, model.time)];
    });
 
    terrain1.flag('uTerrainTexture');
    terrain1.identity().move(0, 0.5, -7).turnX(-.5 * Math.PI).scale(2);
    terrain1.setVertices((u, v) => {
       return [(2 * u - 1) * model.time, (2 * v - 1) * model.time, .07 * cg.noise(2 * u - s, 2 * v * s, .5 * s)];
    });
 
    terrain2.flag('uTerrain2Texture');
    terrain2.identity().move(0, 0.5, -10).turnX(-.5 * Math.PI).scale(.5, .5, 5);
    terrain2.setVertices((u, v) => {
       return [2 * u - 1, 2 * v - 1, .07 * cg.noise(2 * u - 1, 2 * v - 1, .5)];
    });
 
    terrain3.flag('uTerrainTexture');
    terrain3.identity().move(-0.5, 0.5, -3).turnX(-.5 * Math.PI).scale(2);
    terrain3.setVertices((u, v) => {
       return [2 * u - 1, 2 * v - 1, .2 * cg.noise(2 * u - model.time, 2 * v, model.time)];
    });
 
    terrain4.flag('uTerrainTexture');
    terrain4.identity().move(0.5, 0.5, -4).turnX(-.5 * Math.PI).scale(2);
    terrain4.setVertices((u, v) => {
       return [2 * u - 1, 2 * v - 1, .2 * cg.noise(2.5 * u - model.time, 2.5 * v, model.time)];
    });
 
    terrain5.flag('uTerrainTexture');
    terrain5.identity().move(1, 0.5, -5).turnX(-.5 * Math.PI).scale(2);
    terrain5.setVertices((u, v) => {
       return [2 * u - 1, 2 * v - 1, .2 * cg.noise(2.7 * u - model.time, 2.7 * v, model.time)];
    });

    for (let n = 0 ; n < N ; n++)
    {
        for (let i = 0 ; i < 3 ; i++)
        {
            V[n][i] = Math.max(-.05, Math.min(.05, V[n][i] + (Math.random() - .5) * model.deltaTime));
            data[n].p[i] = data[n].p[i] + V[n][i] * Math.sin(model.deltaTime);
            /* data[n].p[i] = Math.max(-.5 , Math.min(.5 , data[n].p[i] + V[n][i] * model.deltaTime)); */
        }
    /* data[n].p = cg.scale(cg.normalize(data[n].p), .5); */
    }
    particles.setParticles(data /*, 'yaw' */);

    /* terrain.identity().turnX(-.3 * Math.PI);
    let axis = Math.sqrt(N);
    terrain.setVertices2((u,v,data) => {
        let i = u * nu + v * nv * axis;
        return [data[i].p[0], data[i].p[1], data[i].p[2]];
    }); */

          /* clay.animateWire(wire, .1, f);
      clay.animateWire(wire2, .1, f);
      clay.animateWire(wire3, .1, f);
      obj.hud().scale(.05,.05,.0001); */

      /* let left0 = buttonState.left[0].pressed,
      left1 = buttonState.left[1].pressed,
      left2 = buttonState.left[2].pressed,
      left3 = buttonState.left[3].pressed,
      left4 = buttonState.left[4].pressed,
      left5 = buttonState.left[5].pressed,
      left6 = buttonState.left[6].pressed,
      right0 = buttonState.right[0].pressed,
      right1 = buttonState.right[1].pressed,
      right2 = buttonState.right[2].pressed,
      right3 = buttonState.right[3].pressed,
      right4 = buttonState.right[4].pressed,
      right5 = buttonState.right[5].pressed,
      right6 = buttonState.right[6].pressed;
      if(left1) console.log('left1');
      if(left2) console.log('left2');
      if(left3) console.log('left3');
      if(left4) console.log('left4');
      if(left5) console.log('left5');
      if(left6) console.log('left6');
      if(right1) console.log('right1');
      if(right2) console.log('right2');
      if(right3) console.log('right3');
      if(right4) console.log('right4');
      if(right5) console.log('right5');
      if(right6) console.log('right6');
      console.log(buttonState.left);
      console.log(buttonState.right); */
      if(buttonState.left[0].pressed)
      {/* down */
         gltf.translation[1] -= 1.5;
         model.move(0,-1.5,0);
         model.scale(1.1, 1.1, 1.1);
      }
      if(buttonState.right[0].pressed)
      {/* up */
         gltf.translation[1] += 1.5;
         model.move(0,1.5,0);
         model.scale(0.9, 0.9, 0.9);
      }
      if(keyboardInput.keyIsDown(keyboardInput.KEY_W))
      {
         model.move(0,1.5,0);
      }
      if(keyboardInput.keyIsDown(keyboardInput.KEY_S))
      {
         model.move(0,-1.5,0);
      }
});
}