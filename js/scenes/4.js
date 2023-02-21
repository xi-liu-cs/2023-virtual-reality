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

   clay.defineMesh('terrain2', clay.createGrid(100, 100));
   let terrain2 = model.add('terrain2').opacity(1); /* island */

   model.animate(() => {
      model.customShader(`
      uniform int uObjTexture;
      uniform int uSkyTexture;
      uniform int uTerrainTexture;
      uniform int uTerrain2Texture;
      vec3 u_sky_color = vec3(.7, .2, .7);

      float turbulence(vec3 p)
      {
         float t = 0., f = 1.;
         for (int i = 0 ; i < 10 ; ++i)
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
         color = 0.05 * height * u_sky_color;
      }
      if(uTerrainTexture == 1)
         color *= .7 + turbulence(5. * vAPos);
      if(uTerrain2Texture == 1)
      {
         color += vec3(.1, .3, .1);
         color *= .5 + pattern(3. * vAPos);
      }
      `);

      obj.flag('uObjTexture');

      roomBackground.flag('uSkyTexture');

      terrain.flag('uTerrainTexture');
      terrain.identity().move(0, 0.3, 0).turnX(-.5 * Math.PI).scale(5);
      terrain.setVertices((u, v) => {
         return [2 * u - 1, 2 * v - 1, .07 * cg.noise(3 * u - .5 * model.time, 3 * v, .5 * model.time)];
      });

      terrain2.flag('uTerrain2Texture');
      terrain2.identity().move(0, 0.5, 0).turnX(-.5 * Math.PI).scale(1, 1, 10).color(.3, .5, .3);
      terrain2.setVertices((u, v) => {
         return [2 * u - 1, 2 * v - 1, .07 * cg.noise(2 * u, 2 * v, .5)];
      });

      obj.hud().scale(.05,.05,.0001);
   });
}