import * as global from "../global.js";
import { Gltf2Node } from "../render/nodes/gltf2.js";

export default () => {
   global.scene().addNode(new Gltf2Node({
      url: "../media/gltf/cave/cave.gltf"
   })).name = "backGround";

   global.scene().addNode(new Gltf2Node({
      url: "../media/gltf/0/human.glb"
   })).name = "human";

   return {
      enableSceneReloading: true,
      scenes: [ 
         {name: "xi_liu" , path: "./1.js"},
         {name: "0", path: "./0.js"}
      ]
   };
}