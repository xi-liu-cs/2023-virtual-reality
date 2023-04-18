import * as cg from "../render/core/cg.js";
import {
    controllerMatrix,
    buttonState,
    joyStickState
} from "../render/core/controllerInput.js";
import {
    lcb,
    rcb
} from "../handle_scenes.js";
import {
    g2
} from "../util/g2.js";

/*********************************************************************************

This demo demonstrates objects with physics.

*********************************************************************************/

let leftTriggerPrev = false;
let rightTriggerPrev = false;

let MP = cg.mTranslate(0, 1, .5);
let A = [0, 1, .5];
let MA = cg.mIdentity();

export const init = async model => {
	model.setRoom(false);
	model.setTable(false);

    // CREATE THE BOX.
    let box = model.add('cube');
	let large_box = model.add('cube');

    let gravity = -9.8;
    let antigravity = false;

    let boxPos_x = 0;
    let boxPos_y = 0;
    let boxPos_z = 0;

    let boxVel_x = 0.03;
    let boxVel_y = 0.01;
    let boxVel_z = 0.015;

    let timestep = 0.05;

    let box_chart = model.add("cube").texture(() => {
        g2.setColor("blue");
        let values = [Math.abs(10 * boxVel_x), Math.abs(10 * boxVel_y), Math.abs(10 * boxVel_z)];
        g2.barChart(.3, .2, .5, .5, values, ["X speed", "Y speed", "Z speed"], ["red", "green", "blue"]);
    });
    g2.addWidget(box_chart, 'button', 0.1, 0.1, "#ff8080", "enable antigravity", () => {
        antigravity = !antigravity;
    });

    // FUNCTION TO RETURN TRUE IF A POINT IS INSIDE THE BOX, OTHERWISE FALSE.

    let isInBox = p => {

        // FIRST TRANSFORM THE POINT BY THE INVERSE OF THE BOX'S MATRIX.

        let q = cg.mTransform(cg.mInverse(box.getMatrix()), p);

        // THEN WE JUST NEED TO SEE IF THE RESULT IS INSIDE A UNIT CUBE.

        return q[0] >= -1 & q[0] <= 1 &&
            q[1] >= -1 & q[1] <= 1 &&
            q[2] >= -1 & q[2] <= 1;
    }


    let beamInBox = (con, b, size) => {
        let point = con == "l" ? lcb.projectOntoBeam(A) : rcb.projectOntoBeam(A);
        let diff = cg.subtract(point, box.getMatrix().slice(12, 15));
        return Math.abs(diff[0]) <= size && Math.abs(diff[1]) <= size && Math.abs(diff[2]) <= size;
        return false;
    }


    model.animate(() => {

        // display the bar chart
        box_chart.identity().move(0, 2, 0).scale(.7, .7, .0001);
        // enable physics, unless you're holding the object

        let physics_on = true;

        // FETCH THE MATRIXES FOR THE LEFT AND RIGHT CONTROLLER.

        let ml = controllerMatrix.left;
        let mr = controllerMatrix.right;

        // EXTRACT THE LOCATION OF EACH CONTROLLER FROM ITS MATRIX,
        // AND USE IT TO SEE WHETHER THAT CONTROLLER IS INSIDE THE BOX.

        let isLeftInBox = beamInBox("l", box, 0.2);
        let isRightInBox = beamInBox("r", box, 0.2);

        // IF NEITHER CONTROLLER['s beam] IS INSIDE THE BOX, COLOR THE BOX WHITE.
        if (!isLeftInBox && !isRightInBox)
            box.color(1, 1, 1);

        // IF THE LEFT CONTROLLER['s beam] IS INSIDE THE BOX

        if (isLeftInBox) {

            // COLOR THE BOX PINK.

            box.color(1, .5, .5);

            // IF THE LEFT TRIGGER IS SQUEEZED

            let leftTrigger = buttonState.left[0].pressed;
            if (leftTrigger) {
                physics_on = false;

                // COLOR THE BOX RED AND MOVE THE BOX.

                box.color(1, 0, 0);
                A = lcb.projectOntoBeam(A);
            }
            leftTriggerPrev = leftTrigger;
        }

        // IF THE RIGHT CONTROLLER['s beam] IS INSIDE THE BOX

        if (isRightInBox) {

            // COLOR THE BOX LIGHT BLUE.

            box.color(.5, .5, 1);

            // IF THE RIGHT TRIGGGER IS SQUEEZED

            let rightTrigger = buttonState.right[0].pressed;
            if (rightTrigger) {
                physics_on = false;

                // COLOR THE BOX BLUE AND MOVE AND ROTATE THE BOX.

                box.color(0, 0, 1);
                let MB = mr.slice();
                if (!rightTriggerPrev) // ON RIGHT DOWN EVENT:
                    MA = MB; // INITIALIZE PREVIOUS MATRIX.
                else
                    MP = cg.mMultiply(cg.mMultiply(MB, cg.mInverse(MA)), MP);

                MA = MB; // REMEMBER PREVIOUS MATRIX.
            }
            rightTriggerPrev = rightTrigger;
        }

        // do physics

        if (physics_on) {
            //let current_pos = MP.slice(12,15);
            let current_pos = [...A];
            let friction = 0.00001;

            boxVel_y += gravity * timestep * 0.001;
            current_pos[0] += boxVel_x;
            current_pos[1] += boxVel_y;
            current_pos[2] += boxVel_z;

            if (Math.abs(current_pos[0]) > 4) {
                current_pos[0] = 4 * Math.sign(current_pos[0]);
                boxVel_x *= -0.8;
            }
            if (Math.abs(current_pos[2]) > 3.8) {
                current_pos[2] = 3.8 * Math.sign(current_pos[2]);
                boxVel_z *= -0.8;
            }
            if (current_pos[1] < 0.2) {
                current_pos[1] = 0.2;
                boxVel_y *= -0.8;
                boxVel_x -= friction * Math.sign(boxVel_x);
                boxVel_z -= friction * Math.sign(boxVel_z);
            }
            if (current_pos[1] > 2) {
                current_pos[1] = 2;
                boxVel_y *= -0.8;
            }

            MP = cg.mMultiply(cg.mTranslate(cg.subtract(current_pos, A)), MP);
            A = current_pos;

        }

        //display box
        box.identity().move(A).scale(.2);
		large_box.identity().move(A).move(10,10,10);
        physics_on = true;

    });
}
