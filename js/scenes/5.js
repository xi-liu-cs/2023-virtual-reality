import * as cg from "../render/core/cg.js";
import {g2} from "../util/g2.js";

export const init = async model => {
   g2.textHeight(.1);
   model.setTable(false);

   let obj = model.add('cube').texture(() => {
      g2.setColor('black');
      g2.fillRect(0, 0, 1, 1);
      g2.setColor('white');
      g2.textHeight(.05);
      g2.fillText('text', .1, .45, 'center');
      g2.fillRect(.07, .35, .2, .05);
      g2.drawWidgets(obj);
   });
   obj.value = [.5,.5];
   g2.addWidget(obj, 'trackpad', .5, .6, '#ff8080', 'trackpad', value => obj.value = value);
   g2.addWidget(obj, 'button' , .5, .1, '#ffffff', 'button', value => {});

   let roomBackground = model.add('roomBackground');

   clay.defineMesh('terrain', clay.createGrid(100, 100));
   let terrain = model.add('terrain').color(0,.5,1).opacity(.7); /* water */

   clay.defineMesh('terrain0', clay.createGrid(100, 100));
   let terrain0 = model.add('terrain0').color(.2,.5,1).opacity(.5); /* water */

   clay.defineMesh('terrain1', clay.createGrid(100, 100));
   let terrain1 = model.add('terrain1').color(.3,.5,1).opacity(.5); /* water */

   clay.defineMesh('terrain2', clay.createGrid(100, 100));
   let terrain2 = model.add('terrain2').opacity(1); /* island */

   let isHUD = false;
   model.control('h', 'toggle HUD', () => isHUD = ! isHUD);

   let N = (t,a,b,c,d) => cg.noise(a * t, b * t, c * t + d * model.time);
   let f = t => [ N(t,3,3.3,3.6,.3), N(t,3.3,3.6,3,.25), N(t,3.6,3,3.3,.2) ];
   let wire = model.add(clay.wire(200,8)).move(.5,.3,0).scale(.5);
   let wire2 = model.add(clay.wire(200,8)).move(0,.3,0).scale(.5);
   let wire3 = model.add(clay.wire(200,8)).move(0,.3,0).scale(.5);

   model.animate(() => {
      let m = views[0]._viewMatrix, c = .5 * Math.cos(model.time), s = .5 * Math.sin(model.time);

      terrain2.flag('uRayTrace');
      model.setUniform('4fv','uL', [.5,.5,.5,1., -.5,-.5,-.5,.2, .7,-.7,0,.2, -.7,.7,0,.2]);
      model.setUniform('4fv','uS', [c,s,0,0, s,0,c,0, 0,c,s,0, -c,-s,0,0]);
      model.setUniform('4fv','uC', [c,s,0,2, s,0,c,2, 0,c,s,2, c,0,s,2]);

      model.customShader(`
      uniform int uObjTexture;
      uniform int uSkyTexture;
      uniform int uTerrainTexture;
      uniform int uTerrain2Texture;
      uniform int uRayTrace;
      uniform vec4 uC[4], uL[4], uS[4];
      vec4 light[4], sphere[4];
      vec3 u_sky_color = vec3(.2, .2, .7);
      vec3 u_floor_color = vec3(.7, .6, .5);

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
            color = 0.05 * height * u_sky_color;
      }

      if(uTerrainTexture == 1)
      {
         for(int l = 0; l < 4; l++)
         {
            vec3 N = normalize(vec3(0., 0., 1.));
            vec3 lDir = light[l].xyz;
            float lBrightness = light[l].w;
            vec3 R = 2. * N * dot(N, lDir) - lDir;
            color += lBrightness * (.9 * max(0., dot(N, lDir)) + vec3(pow(max(0., R.z), 10.)));
         }
         color *= .7 + turbulence(5. * vAPos);
      }

      if(uTerrain2Texture == 1)
      {
         color += vec3(.1, .3, .1);
         color *= .5 + pattern(3. * vAPos);
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

      obj.flag('uObjTexture');

      roomBackground.flag('uSkyTexture');

      function turbulence(p1, p2, p3)
      {
         let t = 0., f = 1.;
         for(let i = 0; i < 10; ++i)
         {
            t += Math.abs(cg.noise(f * p1, f * p2, f * p3)) / f;
            f *= 2.;
         }
         return t;
      }

      function pattern(v1, v2, v3)
      {
         let n = 10,
         res = 0., f = 1.;
         for(let i = 1; i < n; ++i)
         {
            res += cg.noise(f * v1, f * v2, f * v3) / i;
            f = f * i + i;
         }
         return res;
      }

      terrain.flag('uTerrainTexture');
      terrain.identity().move(0, 0.5, 0).turnX(-.5 * Math.PI).scale(5);
      terrain.setVertices((u, v) => {
         return [1000 * (2 * u - 1) * model.time, 1000 * (2 * v - 1) * model.time, .1 * cg.noise(30 * u - .3 * c, 30 * v * c, .3 * c)];
      });

      terrain0.flag('uTerrainTexture');
      terrain0.identity().move(0, 0.5, 0).turnX(-.5 * Math.PI).scale(1, 1, 10);
      terrain0.setVertices((u, v) => {
         return [2 * u - 1, 2 * v - 1, .2 * cg.noise(3 * u - model.time, 3 * v, model.time)];
      });

      terrain1.flag('uTerrainTexture');
      terrain1.identity().move(0, 0.5, 0).turnX(-.5 * Math.PI).scale(5);
      terrain1.setVertices((u, v) => {
         return [(2 * u - 1) * model.time, (2 * v - 1) * model.time, .07 * cg.noise(2 * u - s, 2 * v * s, .5 * s)];
      });

      terrain2.flag('uTerrain2Texture');
      terrain2.identity().move(0, 0.5, 0).turnX(-.5 * Math.PI).scale(1, 1, 10).color(.3, .5, .3);
      terrain2.setVertices((u, v) => {
         return [2 * u - 1, 2 * v - 1, .07 * cg.noise(2 * u, 2 * v, .5)];
      });

      clay.animateWire(wire, .1, f);
      clay.animateWire(wire2, .1, f);
      clay.animateWire(wire3, .1, f);
      wire.flag('uRayTrace');
      wire2.flag('uRayTrace');
      wire3.flag('uRayTrace');
      obj.hud().scale(.05,.05,.0001);
   });
}