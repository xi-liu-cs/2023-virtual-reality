import { g2 } from "../util/g2.js";

export const init = async model => {

   g2.textHeight(.1);
   let obj1 = model.add('cube').texture(() => g2.clock(0,0,1,1)); // HUD object.
   let obj2 = model.add('cube').texture(() => g2.clock(0,0,1,1)); // non-HUD object.
   // let box = model.add('cube');
   model.setTable(false);
   model.animate(() => {
      // console.log('flag', box.flag);
      // box.flag('uNoiseTexture');
      // model.customShader(`
      // uniform int uNoiseTexture;
      // ---
      // if(uNoiseTexture == 1)
      //    color *= .5 + noise(3. * vAPos);
      // `);
      // box.identity().move(0, 1.6, 0).turnY(model.time).turnX(model.time).scale(.1);
      obj1.hud().scale(.2,.2,.0001);
      // obj2.identity().move(0,1.5,-1).scale(.2,.2,.0001);
   });
}