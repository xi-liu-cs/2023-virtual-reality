import * as cg from "../render/core/cg.js";
import {g2} from "../util/g2.js";
import * as global from "../global.js";
import {Gltf2Node} from "../render/nodes/gltf2.js";

export const init = async model => {
   g2.textHeight(.1);
   model.setTable(false);

   let gltf = new Gltf2Node({url: './media/gltf/0/mountain.glb'});
   gltf.translation = [0, -1, -20];
   global.gltfRoot.addNode(gltf);

   /* let obj = model.add('cube').texture(() => {
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
   g2.addWidget(obj, 'button' , .5, .1, '#ffffff', 'button', value => {}); */

   let roomBackground = model.add('roomBackground');

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

   let isHUD = false;
   model.control('h', 'toggle HUD', () => isHUD = ! isHUD);

   model.animate(() => {
      let m = views[0]._viewMatrix, c = .5 * Math.cos(model.time), s = .5 * Math.sin(model.time);

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

      /* obj.flag('uObjTexture'); */

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


      /* clay.animateWire(wire, .1, f);
      clay.animateWire(wire2, .1, f);
      clay.animateWire(wire3, .1, f);
      obj.hud().scale(.05,.05,.0001); */
   });
}