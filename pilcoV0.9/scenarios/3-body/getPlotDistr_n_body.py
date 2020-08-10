"""
%% getPlotDistr_cp.m
% *Summary:* Compute means and covariances of the Cartesian coordinates of
% the tips both the inner and outer pendulum assuming that the joint state
% $x$ of the cart-double-pendulum system is Gaussian, i.e., $x\sim N(m, s)$
%
%
%     function [M, S] = getPlotDistr_cp(m, s, ell)
%
%
%
% *Input arguments:*
%
%   m       mean of full state                                    [4 x 1]
%   s       covariance of full state                              [4 x 4]
%   ell     length of pendulum
%
%   Note: this code assumes that the following order of the state:
%          1: cart pos., 
%          2: cart vel., 
%          3: pendulum angular velocity,
%          4: pendulum angle
%
% *Output arguments:*
%
%   M      mean of tip of pendulum                               [2 x 1]
%   S      covariance of tip of pendulum                         [2 x 2]
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modification: 2013-03-27
%
%% High-Level Steps
% # Augment input distribution to complex angle representation
% # Compute means of tips of pendulums (in Cartesian coordinates)
% # Compute covariances of tips of pendulums (in Cartesian coordinates) 

"""

@mfunction("M, S")
def getPlotDistr_cp(m=None, s=None, ell=None):
    #% Code

    # 1. Augment input distribution to complex angle representation
    [m1, s1, c1] = gTrig(m, s, 4, ell)# map input distribution through sin/cos
    m1 = mcat([m, OMPCSEMI, m1])# mean of joint
    c1 = s * c1# cross-covariance between input and prediction
    s1 = mcat([s, c1, OMPCSEMI, c1.cT, s1])# covariance of joint

    # 2. Compute means of tips of pendulums (in Cartesian coordinates)
    M = mcat([m1(1) + m1(5), OMPCSEMI, -m1(6)])

    # 3. Compute covariances of tips of pendulums (in Cartesian coordinates) 
    s11 = s1(1, 1) + s1(5, 5) + s1(1, 5) + s1(5, 1)# x+l sin(theta)
    s22 = s1(6, 6)# -l*cos(theta)
    s12 = -(s1(1, 6) + s1(5, 6))# cov(x+l*sin(th), -l*cos(th)

    S = mcat([s11, s12, OMPCSEMI, s12.cT, s22])
        try:
        chol(S)

