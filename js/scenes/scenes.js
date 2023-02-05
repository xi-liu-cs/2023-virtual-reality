import * as global from "../global.js";
import { Gltf2Node } from "../render/nodes/gltf2.js";

export default () => {
   global.scene().addNode(new Gltf2Node({
      url: "../media/gltf/0/1.glb"
      // url: "../media/gltf/60_fifth_ave/60_fifth_ave.gltf"
   })).name = "backGround";

   return {
      enableSceneReloading: true,
      scenes: [ 
         {name: "xi_liu" , path: "./1.js"},
         {name: '0', path: './0.js'}
      ]
   };
}