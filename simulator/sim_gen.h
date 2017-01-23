#ifndef SIM_GEN_H
#define SIM_GEN_H
/*
 * construct a simulator
 * Feng Li, lifen@iupui.edu
 * Fist Created at Jan 18, 2017
 *
 * Notes:
 *  1. all attributes will be stored in flat format(continuous),
 *  2. slowest in the x and fastest in the z directions
 *  3. vel data will be aligned in format as: vx0, vy0, vz0, vx1, vy1, vz1
 */

/*
 * simulation process
 * input:
 *  timestep: current timestep
 *  dim, data layout 
 * output:
 *  vel_data & pressure
 *
 * Notes: vel will be (y*time,0,0)
 */
void update_attributes(int timeStep, int *dim,float *vel_data, float*p_data);



#endif
