import * as cg from "../render/core/cg.js";
import { g2 } from "../util/g2.js";
import { controllerMatrix, buttonState, joyStickState } from "../render/core/controllerInput.js";
import { lcb, rcb } from '../handle_scenes.js';

let center = [0, 1.5, 0],
radius = 0.1,
tex = '../media/textures/planet.jpg';

function particle()
{
    let max_particles = 125,
    domain_width = 40,
    domain_height = 80,

    particle_mass = 1,
    isotropic_exponent = 20,
    base_density = 1,
    smoothing_length = 5,
    dynamic_viscosity = 0.5,
    damping_coefficient = -0.9,
    constant_force = [0, -0.1],
    
    time_step_length = 0.01,
    n_time_steps = 2_500,
    add_particles_every = 50,
    
    figure_size = [4, 6],
    plot_every = 6,
    scatter_dot_size = 2_000,
    
    domain_x_lim = [smoothing_length, domain_width - smoothing_length],
    domain_y_lim = [smoothing_length, domain_height - smoothing_length],

    n_particles = 1;
}

function color(i, j, t)
{
  let buf = new Array(16);
  for(let i = 0; i <= 9; ++i)
    buf[i] = String.fromCharCode('0'.charCodeAt(0) + i);
  for(let i = 10; i <= 15; ++i)
    buf[i] = String.fromCharCode('a'.charCodeAt(0) + i - 10);
  let buf_n = buf.length,
  format_n = '#000000'.length,
  str = '#' + (Math.abs(parseInt(i * 100 + j * 10 + t))).toString(),
  str_n = str.length;
  if(str_n < format_n)
  {
    let diff = format_n - str_n;
    for(let idx = 0; idx < diff; ++idx)
      str += buf[(i * j * idx) % buf_n].toString().slice(0, 1);
  }
  else if(str_n > format_n)
    str = str.slice(0, format_n);
  return str;
}

export const init = async model =>
{
   let box = model.add('cube').texture(tex);
   let obj1 = model.add('cube').texture(() =>
   {
        g2.setColor('black');
        g2.fillRect(0, 0, 1, 1);
        g2.setColor('white');
        g2.textHeight(.05);
        g2.fillText('text', .1, .45, 'center');
        g2.fillRect(.07, .35, .2, .05);
        g2.drawWidgets(obj1);
    });
    obj1.value = [.5, .5];

    let scale_factor = 1;

    g2.addWidget(obj1, 'trackpad', .5, .6, '#ff8080', 'trackpad', value => obj1.value = value);
    g2.addWidget(obj1, 'button' , .5, .1, '#ffffff', 'button', value => {model.add('smooth_octahedron').texture(tex);});

    let obj2 = model.add('cube').texture(() =>
    {
        g2.setColor('black');
        g2.textHeight(.1);
  
        g2.setColor('#ff8080');
        g2.lineWidth(.05);
        let path = [];
        for (let i = 0 ; i < 100 ; i++)
           path.push([.5 + .1 * Math.cos(.08 * i - 2 * model.time), .66 - .005 * i]);
        g2.drawPath(path);
    });

    let obj3 = model.add('cube').texture(() =>
    {
        g2.setColor('white');
        g2.fillRect(0, 0, .1, .1);
    });

    let a = [-1, 0, 0, -1, 0, 0], A = [ 1, 0, 0,  1, 0, 0],
    b = [ 0,-1, 0,  0,-1, 0], B = [ 0, 1, 0,  0, 1, 0],
    c = [ 0, 0,-1,  0, 0,-1], C = [ 0, 0, 1,  0, 0, 1];
    clay.defineMesh('smooth_octahedron', clay.trianglesMesh([
    a,b,C, a,B,c, A,b,c, A,B,C, a,C,B, A,c,B, A,C,b, a,c,b
    ]));
    let obj4 = model.add('smooth_octahedron').texture(tex);

   model.animate(() =>
   {
        let point = rcb.projectOntoBeam(center),
        diff = cg.subtract(point, center),
        hit = cg.norm(diff) < radius,
        rt = buttonState.right[0].pressed; /* right trigger */
        if(hit && rt) center = point;
        box.color(hit ? rt ? [1, 0, 0] : [1, .5, .5] : [1, 1, 1]);
        box.identity().move(center).scale(radius);
        obj1.identity().move(0,1.2,0).scale(.5,.5,.0001);
        obj4.identity().move(0, 1.2, 0).scale(.1, .1, .1).turnZ(model.time/2).turnY(model.time/2).turnZ(model.time);
   });
}