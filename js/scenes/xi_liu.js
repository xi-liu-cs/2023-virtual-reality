import { controllerMatrix, buttonState, joyStickState } from "../render/core/controllerInput.js";

export const init = async model => 
{
   let c = [.8,.5,.3],
   a1 = model.add('tubeY').color(c),
   a2 = model.add('tubeY').color(c),
   a2_2 = model.add('sphere').color(c),
   a2_3 = model.add('sphere').color(c),
   a3 = model.add('sphere').color(c),
   a3_2 = model.add('sphere').color(c),
   a3_3 = a3.add('sphere').color(c),
   a3_4 = a3_2.add('sphere').color(c),
   a3_5 = a3_3.add('cube').color(c),
   a3_6 = a3_4.add('cube').color(c),
   a4 = a2.add('sphere').color(c),
   a4_2 = a2.add('sphere').color(c),
   a4_low = a4.add('sphere').color(c),
   a4_low_2 = a4_2.add('sphere').color(c),
   eye = model.add(),
   head = eye.add('sphere').color(c);

   model.move(0,1.5,0).scale(.3).animate(() => 
   {
      model.turnY(.1 * Math.cos(2.1*model.time));

      /* neck */
      a1.identity().move(0, .5, 0).scale(.2);

      /* body */
      a2.identity().move(0, 0, 0).scale(.2, .6, .2);
      a2_2.identity().move(0, -.2, 0).scale(.35, .6, .3);
      a2_3.identity().move(0, .2, 0).scale(.35, .6, .3);

      /* upper leg */
      a3.identity().move(-.2, -1, 0).turnX(.25 * Math.sin(2 * model.time)).scale(.13, .7, .13);
      a3_2.identity().move(.2, -1, 0).turnX(-.25 * Math.sin(2 * model.time)).scale(.13, .7, .13);

      /* lower leg */
      a3_3.identity().move(0, -1.5, 0).scale(.7, .7, .7);
      a3_4.identity().move(0, -1.5, 0).scale(.7, .7, .7);

      /* feet */
      a3_5.identity().move(-.2, -1, 1).scale(.9, .06, 2);
      a3_6.identity().move(.2, -1, 1).scale(.9, .06, 2);

      /* upper arm */
      a4.identity().move(-2, .17, 0).turnY(-Math.PI / 2.5).turnX(Math.PI / 2.5).turnY(.5 * Math.sin(2 * model.time)).scale(.4, .4, .7);
      a4_2.identity().move(2, .17, 0).turnY(Math.PI / 2.5).turnX(Math.PI / 2.5).turnY(.5 * Math.sin(2 * model.time)).scale(.4, .4, .7);

      /* lower arm */
      a4_low.identity().move(0, .17, 1.5).turnX(Math.PI / 3);
      a4_low_2.identity().move(0, .17, 1.5).turnX(Math.PI / 3);

      head.identity().move(0,1,0)
                    .turnX(Math.sin(2.1*model.time))
                    .turnY(Math.sin(model.time))
                    .turnY(.1 * Math.cos(2.1*model.time))
                    .scale(.3);
   });
}