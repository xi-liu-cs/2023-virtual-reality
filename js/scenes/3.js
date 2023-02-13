import * as cg from "../render/core/cg.js";
import { g2 } from "../util/g2.js";
import { controllerMatrix, buttonState, joyStickState } from "../render/core/controllerInput.js";
import { lcb, rcb } from '../handle_scenes.js';

let center = [0, 1.5, 0],
radius = 0.1,
tex = '../media/textures/planet.jpg';

function is_hit(center, radius)
{
    let point = rcb.projectOntoBeam(center),
    diff = cg.subtract(point, center),
    hit = cg.norm(diff) < radius;
    return [hit, point];
}

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
    damping_coefficient = - 0.9,
    constant_force = [0, -0.1],
    
    time_step_length = 0.01,
    n_time_steps = 2_500,
    add_particles_every = 50,
    
    figure_size = [4, 6],
    plot_every = 6,
    scatter_dot_size = 2_000,
    
    domain_x_lim = [smoothing_length, domain_width - smoothing_length],
    domain_y_lim = [smoothing_length, domain_height - smoothing_length];

    n_particles = 1;
}

export const init = async model =>
{
   let box = model.add('cube').texture(tex);

   let a = model.add('cube').texture(() =>
   {
        g2.setColor('black');
        g2.fillRect(0, 0, .5, .5);
        g2.setColor('white');
        g2.textHeight(.05);
        g2.fillText('particle', .2, .45, 'center');
        g2.drawWidgets(a);
    });
    a.value = [.5, .5];
    g2.addWidget(a, 'button', .03, .04, 'white', 'number of particles', () => { });

   model.animate(() =>
   {
      let [hit, point] = is_hit(center, radius),
      rt = buttonState.right[0].pressed; /* right trigger */
      if(hit && rt) center = point;
      box.color(hit ? rt ? [1, 0, 0] : [1, .5, .5] : [1, 1, 1]);
      box.identity().move(center).scale(radius);
      a.identity().move(0,1.2,0).scale(.5,.5,.0001);
   });
}