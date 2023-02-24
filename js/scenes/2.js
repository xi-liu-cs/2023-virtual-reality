import * as cg from "../render/core/cg.js";
import { controllerMatrix, buttonState, joyStickState } from "../render/core/controllerInput.js";
import { lcb, rcb } from '../handle_scenes.js';

let leftTriggerPrev = false;
let rightTriggerPrev = false;
let right1pre = false;
let right3pre = false;

let MP = cg.mTranslate(0,1,.5);
let MP2 = cg.mTranslate(0,-1,.5);
let A = [0,0,0];
let MA = cg.mIdentity();
let MA1 = cg.mIdentity();
let MA2 = cg.mIdentity();
let tex = '../media/textures/planet.jpg';

export const init = async model => {

   // CREATE THE BOX.

   let box = model.add('tubeY').texture(tex);
   let obj;

   // FUNCTION TO RETURN TRUE IF A POINT IS INSIDE THE BOX, OTHERWISE FALSE.

   let isInBox = p => {

      // FIRST TRANSFORM THE POINT BY THE INVERSE OF THE BOX'S MATRIX.

      let q = cg.mTransform(cg.mInverse(box.getMatrix()), p);

      // THEN WE JUST NEED TO SEE IF THE RESULT IS INSIDE A UNIT CUBE.

      /* return q[0] >= -1 & q[0] <= 1 &&
             q[1] >= -1 & q[1] <= 1 &&
             q[2] >= -1 & q[2] <= 1 ; */
      
      return q[0] >= -2 & q[0] <= 2 &&
      q[1] >= -2 & q[1] <= 2 &&
      q[2] >= -2 & q[2] <= 2 ;
   }

   let pre_time = model.time;
   model.animate(() => {

      // FETCH THE MATRIXES FOR THE LEFT AND RIGHT CONTROLLER.

      let ml = controllerMatrix.left;
      let mr = controllerMatrix.right;

      // EXTRACT THE LOCATION OF EACH CONTROLLER FROM ITS MATRIX,
      // AND USE IT TO SEE WHETHER THAT CONTROLLER IS INSIDE THE BOX.

      let isLeftInBox  = isInBox(ml.slice(12,15));
      let isRightInBox = isInBox(mr.slice(12,15));

      // IF NEITHER CONTROLLER IS INSIDE THE BOX, COLOR THE BOX WHITE.

      if (! isLeftInBox && ! isRightInBox)
         box.color(1,1,1);

      // IF THE LEFT CONTROLLER IS INSIDE THE BOX

      if (isLeftInBox) {

         // COLOR THE BOX PINK.

         box.color(1,.5,.5);

         // IF THE LEFT TRIGGER IS SQUEEZED

         let leftTrigger = buttonState.left[0].pressed;
	 if (leftTrigger) {

            // COLOR THE BOX RED AND MOVE THE BOX.

            box.color(1,0,0);
            let B = ml.slice(12,15);
            if (! leftTriggerPrev)         // ON LEFT DOWN EVENT:
               A = B;                      // INITIALIZE PREVIOUS LOCATION.
            else
               MP = cg.mMultiply(cg.mTranslate(cg.subtract(B, A)), MP);

	    A = B;                         // REMEMBER PREVIOUS LOCATION.
         }
         leftTriggerPrev = leftTrigger;
      }

      // IF THE RIGHT CONTROLLER IS INSIDE THE BOX

      if (isRightInBox) {

         // COLOR THE BOX LIGHT BLUE.

         box.color(.5,.5,1);

	 // IF THE RIGHT TRIGGGER IS SQUEEZED

         let rightTrigger = buttonState.right[0].pressed;
         /* let right1 = buttonState.right[1].pressed,
         right2 = buttonState.right[2].pressed,
         right3 = buttonState.right[3].pressed,
         right4 = buttonState.right[4].pressed,
         right5 = buttonState.right[5].pressed,
         right6 = buttonState.right[6].pressed;
         if(right1) box.color(1,0,0); // lower left elliptical button
         if(right2) box.color(0,1,0);
         if(right3) box.color(0,0,0); // thumbstick
         if(right4) box.color(1,1,0); // a
         if(right5) box.color(0,1,1); // b
         if(right6) box.color(1,0,1); */
	 if (rightTrigger) {

	    // COLOR THE BOX BLUE AND MOVE AND ROTATE THE BOX.

            box.color(0,0,1);
            let MB = mr.slice();
            if (! rightTriggerPrev)        // ON RIGHT DOWN EVENT:
               MA = MB;                    // INITIALIZE PREVIOUS MATRIX.
            else
	       MP = cg.mMultiply(cg.mMultiply(MB, cg.mInverse(MA)), MP);

	    MA = MB;                       // REMEMBER PREVIOUS MATRIX.
         }
         rightTriggerPrev = rightTrigger;
   if(buttonState.right[1].pressed)
   {
      box.color(1, 0, 0);
      model.add('donut').texture(tex);
      obj = model.add('cube').texture(tex);
      let MB1 = mr.slice();
      if (!right1pre)                // ON RIGHT DOWN EVENT:
      {
         box.add('tubeX').texture(tex);
         MA1 = MB1;                  // INITIALIZE PREVIOUS MATRIX.
      }                    
      else
         MP2 = cg.mMultiply(cg.mMultiply(MB1, cg.mInverse(MA1)), MP2);
      MA1 = MB1;      
   }
   if(buttonState.right[3].pressed) // thumbstick
   {
      box.color(0, 0, 0);
      let MB2 = mr.slice();
      if (!right3pre)                // ON RIGHT DOWN EVENT:
      {
         box.add('tubeZ').texture(tex);
         MA2 = MB2;                  // INITIALIZE PREVIOUS MATRIX.
      }                    
      else
         MP = cg.mMultiply(cg.mMultiply(MB2, cg.mInverse(MA2)), MP);
      MA2 = MB2;        
   }
   right3pre = buttonState.right[3].pressed;
   if(buttonState.right[4].pressed) // a
   {
      box.color(1, 1, 0);
      MP = cg.mMultiply(MP, cg.mScale(1.1, 1.1, 1.1));
   }
   if(buttonState.right[5].pressed) // b
   {
      box.color(0, 1, 1);
      MP = cg.mMultiply(MP, cg.mScale(.9, .9, .9));
   }
      }

      // DISPLAY THE BOX.
      box.setMatrix(MP).scale(.2); 
      obj.setMatrix(MP2).scale(.2);
   });
}