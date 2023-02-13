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
   left_hand = a4_low.add('sphere').color(c),
   right_hand = a4_low_2.add('sphere').color(c),
   a = model.add(),
   head = a.add('sphere').color(c),
   leye = head.add('sphere'),
   reye = head.add('sphere'),
   leye_center = leye.add('sphere'),
   reye_center = reye.add('sphere'),
   nose = head.add('sphere'),
   mouth = head.add('sphere');

   model.move(0,1.5,0).scale(.3).animate(() => 
   {
      function clamp(min, max, x){ return Math.min(Math.max(x, min), max); }
      let c = [clamp(.5, 1, Math.cos(model.time)), clamp(.4, 1, Math.sin(model.time)), clamp(.3, 1, Math.cos(model.time))],
      tex = '../media/textures/blue.jpg';
      model.turnY(.1 * Math.cos(2.1*model.time));

      /* neck */
      a1.identity().move(0, .5, 0).scale(.2).color(c).texture(tex);

      /* body */
      a2.identity().move(0, 0, 0).scale(.2, .6, .2).color(c).texture(tex);
      a2_2.identity().move(0, -.2, 0).scale(.35, .6, .3).color(c).texture(tex);
      a2_3.identity().move(0, .2, 0).scale(.35, .6, .3).color(c).texture(tex);

      /* upper leg */
      a3.identity().move(-.2, -1, 0).turnX(.25 * Math.sin(2 * model.time)).scale(.13, .7, .13).color(c).texture(tex);
      a3_2.identity().move(.2, -1, 0).turnX(-.25 * Math.sin(2 * model.time)).scale(.13, .7, .13).color(c).texture(tex);

      /* lower leg */
      a3_3.identity().move(0, -1.5, 0).scale(.7, .7, .7).color(c).texture(tex);
      a3_4.identity().move(0, -1.5, 0).scale(.7, .7, .7).color(c).texture(tex);

      /* feet */
      a3_5.identity().move(-.2, -1, 1).scale(.9, .06, 2).color(c).texture(tex);
      a3_6.identity().move(.2, -1, 1).scale(.9, .06, 2).color(c).texture(tex);

      /* upper arm */
      a4.identity().move(-2, .17, 0).turnY(-Math.PI / 2.5).turnX(Math.PI / 2.5).turnY(.5 * Math.sin(2 * model.time)).scale(.4, .4, .7).color(c).texture(tex);
      a4_2.identity().move(2, .17, 0).turnY(Math.PI / 2.5).turnX(Math.PI / 2.5).turnY(.5 * Math.sin(2 * model.time)).scale(.4, .4, .7).color(c).texture(tex);

      /* lower arm */
      a4_low.identity().move(0, .17, 1.5).turnX(Math.PI / 3).color(c).texture(tex);
      a4_low_2.identity().move(0, .17, 1.5).turnX(Math.PI / 3).color(c).texture(tex);

      /* hand */
      left_hand.identity().move(0, .4, .5).turnZ(Math.PI / 4).scale(.4, .6, .9).color(c).texture(tex);
      right_hand.identity().move(0, .4, .5).turnZ(Math.PI / 4).scale(.4, .6, .9).color(c).texture(tex);

      head.identity().move(0,1,0)
                    .scale(.95, 1.1, .95)
                    .turnY(clamp(-.5, .5, Math.cos(model.time)))
                    .scale(.3)
                    .color(c)
                    .texture(tex);
      
      leye.identity().move(-.4, .3, .9).scale(.2, .1, .1).color(1, 1, 1);
      reye.identity().move(.4, .3, .9).scale(.2, .1, .1).color(1, 1, 1);
      leye_center.identity().scale(.7, 1.1, 1.1).color(0, 0, 0);
      reye_center.identity().scale(.7, 1.1, 1.1).color(0, 0, 0);
      nose.identity().move(0, 0, .1).scale(.4, .6, .97).color(.5, .1, .1);
      mouth.identity().move(0, -.4, .1).scale(.6, .23, .93).color(.5, .1, .1);
   });
}