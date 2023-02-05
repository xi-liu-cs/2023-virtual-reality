import * as cg from "../render/core/cg.js";
import { controllerMatrix, buttonState, joyStickState } from "../render/core/controllerInput.js";

let leftTriggerPrev = false;
let rightTriggerPrev = false;

let MP = cg.mTranslate(0,1,.5);
let A = [0,0,0];
let MA = cg.mIdentity();

export const init = async model => {

   // CREATE THE BOX.

   let box = model.add('sphere');

   // FUNCTION TO RETURN TRUE IF A POINT IS INSIDE THE BOX, OTHERWISE FALSE.

   let isInBox = p => {

      // FIRST TRANSFORM THE POINT BY THE INVERSE OF THE BOX'S MATRIX.

      let q = cg.mTransform(cg.mInverse(box.getMatrix()), p);

      // THEN WE JUST NEED TO SEE IF THE RESULT IS INSIDE A UNIT CUBE.

      return q[0] >= -1 & q[0] <= 1 &&
             q[1] >= -1 & q[1] <= 1 &&
             q[2] >= -1 & q[2] <= 1 ;
   }

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
      }

      // DISPLAY THE BOX.
   
      box.setMatrix(MP).scale(.1);
   });
}

// import 'https://cdn.jsdelivr.net/npm/@webxr-polyfill/polyfill@0.2.2/dist/polyfill.min.js';

// navigator.xr.requestDevice().then(device => {
//     device.supportsSession({ mode: 'immersive-vr' }).then(function (supported) {
//       if (supported) {
//         // Start an XR session
//         device.requestSession({ mode: 'immersive-vr' }).then(session => {
//           let xrSession = session;
//           // Get reference to XR controller
//           let controller = xrSession.inputSources[0].gamepad;
  
//           // Check for button events on the XR controller
//           controller.addEventListener("buttonchange", (event) => {
//             // Check if the button is being pressed
//             if (event.pressed) {
//               // Check which button is being pressed
//               switch (event.button) {
//                 case 0: // Button 0, move XR origin along the x-axis
//                   xrSession.updateWorldTrackingState({
//                     transform: new XRRigidTransform({
//                       position: new DOMPoint(xrSession.worldInformation.identity.position.x + 0.1, 0, 0)
//                     })
//                   });
//                   break;
//                 case 1: // Button 1, move XR origin along the y-axis
//                   xrSession.updateWorldTrackingState({
//                     transform: new XRRigidTransform({
//                       position: new DOMPoint(0, xrSession.worldInformation.identity.position.y + 0.1, 0)
//                     })
//                   });
//                   break;
//                 case 2: // Button 2, move XR origin along the z-axis
//                   xrSession.updateWorldTrackingState({
//                     transform: new XRRigidTransform({
//                       position: new DOMPoint(0, 0, xrSession.worldInformation.identity.position.z + 0.1)
//                     })
//                   });
//                   break;
//                 default:
//                   break;
//               }
//             }
//           });
//         });
//       }
//     });
//   });