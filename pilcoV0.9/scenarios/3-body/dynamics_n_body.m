%% dynamics_cp.m
% *Summary:* Implements ths ODE for simulating the cart-pole dynamics.
%
%    function dz = dynamics_cp(t, z, f)
%
%
% *Input arguments:*
%
%		t     current time step (called from ODE solver)
%   z     state                                                    [4 x 1]
%   f     (optional): force f(t)
%
% *Output arguments:*
%
%   dz    if 3 input arguments:      state derivative wrt time
%         if only 2 input arguments: total mechanical energy
%
%
% Note: It is assumed that the state variables are of the following order:
%       x:        [m]     position of cart
%       dx:       [m/s]   velocity of cart
%       dtheta:   [rad/s] angular velocity
%       theta:    [rad]   angle
%
%
% A detailed derivation of the dynamics can be found in:
%
% M.P. Deisenroth:
% Efficient Reinforcement Learning Using Gaussian Processes, Appendix C,
% KIT Scientific Publishing, 2010.
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-08
% This is where math goes

function dz = dynamics_cp(t,z)
%% Code

% Variables

G=6.67408e-11 %N-m2/kg2     Reference quantities
m_nd=1.989e+30 %kg      mass of the sun
r_nd=5.326e+12 %m       distance between stars in Alpha Centauri
v_nd=30000 %m/s     relative velocity of earth around the sun
t_nd=79.91*365*24*3600*0.51 %s  orbital period of Alpha Centauri#Net constants
K1=G*t_nd*m_nd/(r_nd**2*v_nd)
K2=v_nd*t_nd/r_nd

m1=1.0 
m2=1.0 
m3 = 1.0


if nargin==2
  r1=z(1, 1:3)
  r2=z(1, 4:6)
  r3=z(1, 7:9)
  v1=z(1, 10:12)
  v2=z(1, 13:15)
  v3=z(1, 16:18)

  r12=norm(r2-r1)
  r13=norm(r3-r1)
  r23=norm(r3-r2)

  dv1bydt=K1*m2*(r2-r1)/r12**3+K1*m3*(r3-r1)/r13^3
  dv2bydt=K1*m1*(r1-r2)/r12**3+K1*m3*(r3-r2)/r23^3
  dv3bydt=K1*m1*(r1-r3)/r13**3+K1*m2*(r2-r3)/r23^3
  dr1bydt=K2*v1
  dr2bydt=K2*v2
  dr3bydt=K2*v3
   
  r12_derivs = [dr1bydt dr2bydt]
  r_derivers = [r12_derivs dr3bydt]
  v12_derivs = [dv1bydt dv2bydt]
  v_derivs = [v12_derivs dv3bydt]
  dz = [r_derivs v_derivs]

end
