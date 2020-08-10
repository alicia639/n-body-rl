# fix function calls
"""
%% settings_cp.m
% *Summary:* Script set up the cart-pole scenario
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-05-24
%
%% High-Level Steps
% # Define state and important indices
% # Set up scenario
% # Set up the plant structure
% # Set up the policy structure
% # Set up the cost structure
% # Set up the GP dynamics model structure
% # Parameters for policy optimization
% # Plotting verbosity
% # Some array initializations

%% Code
"""
# fix

rand(mstring('seed'), 1)
randn(mstring('seed'), 1)
format("short")
format("compact")

# include some paths
    try:
    rd = mstring('../../')
    addpath(mcat([rd, mstring('base')]), mcat([rd, mstring('util')]), mcat([rd, mstring('gp')]), mcat([rd, mstring('control')]), mcat([rd, mstring('loss')]))


"""
% 1. Define state and important indices

% 1a. Full state representation (including all augmentations)
%
%  1  x          cart position
%  2  v          cart velocity
%  3  dtheta     angular velocity
%  4  theta      angle of the pendulum
%  5  sin(theta) complex representation ...
%  6  cos(theta) of theta
%  7  u          force applied to cart
%

% 1b. Important indices
% odei  indicies for the ode solver
% augi  indicies for variables augmented to the ode variables
% dyno  indicies for the output from the dynamics model and indicies to loss
% angi  indicies for variables treated as angles (using sin/cos representation)
% dyni  indicies for inputs to the dynamics model
% poli  indicies for the inputs to the policy
% difi  indicies for training targets that are differences (rather than values)
"""

odei = mcat([1, 2, 3, 4])# varibles for the ode solver
augi = mcat([])# variables to be augmented
dyno = mcat([1, 2, 3, 4])# variables to be predicted (and known to loss)
angi = mcat([4])# angle variables
dyni = mcat([1, 2, 3, 5, 6])# variables that serve as inputs to the dynamics GP
poli = mcat([1, 2, 3, 5, 6])# variables that serve as inputs to the policy
difi = mcat([1, 2, 3, 4])# variables that are learned via differences


# 2. Set up the scenario
dt = 0.10# [s] sampling time
T = 4.0# [s] initial prediction horizon time
H = ceil(T / dt)# prediction steps (optimization horizon)
mu0 = mcat([0, 0, 0, 0]).cT# initial state mean
S0 = diag(mcat([0.1, 0.1, 0.1, 0.1]) **elpow** 2)# initial state covariance
N = 15# number controller optimizations
J = 1# initial J trajectories of length H
K = 1# no. of initial states for which we optimize
nc = 10# number of controller basis functions

# 3. Plant structure: fix this part
plant.dynamics = @# dynamics ode function
plant.noise = diag(ones(1, 4) * 0.01 **elpow** 2)# measurement noise
plant.dt = dt
plant.ctrl = @# controler is zero order hold
plant.odei = odei
plant.augi = augi
plant.angi = angi
plant.poli = poli
plant.dyno = dyno
plant.dyni = dyni
plant.difi = difi
plant.prop = @


# 4. Policy structure
policy.fcn = lambda policy, m, s: conCat(@, @, policy, m, s)# controller
# representation
policy.maxU = 10# max. amplitude of
# control
[mm, ss, cc] = gTrig(mu0, S0, plant.angi)# represent angles
mm = mcat([mu0, OMPCSEMI, mm])
cc = S0 * cc; print cc
ss = mcat([S0, cc, OMPCSEMI, cc.cT, ss])
# in complex plane
policy.p.inputs = gaussian(mm(poli), ss(poli, poli), nc).cT# init. location of
# basis functions
policy.p.targets = 0.1 * randn(nc, length(policy.maxU))# init. policy targets
# (close to zero)
policy.p.hyp = log(mcat([1, 1, 1, 0.7, 0.7, 1, 0.01])).cT# initialize policy
# hyper-parameters


# 5. Set up the cost structure
cost.fcn = @# cost function
cost.gamma = 1# discount factor
cost.p = 0.5# length of pendulum
cost.width = 0.25# cost function width
cost.expl = 0.0# exploration parameter (UCB)
cost.angle = plant.angi# index of angle (for cost function)
cost.target = mcat([0, 0, 0, pi]).cT# target state

# 6. Dynamics model structure


dynmodel.fcn = @# function for GP predictions
dynmodel.train = @# function to train dynamics model
dynmodel.induce = zeros(300, 0, 1)# shared inducing inputs (sparse GP)
trainOpt = mcat([300, 500])# defines the max. number of line searches
# when training the GP dynamics models
# trainOpt(1): full GP,
# trainOpt(2): sparse GP (FITC)


# 7. Parameters for policy optimization
opt.length = 150# max. number of line searches
opt.MFEPLS = 30# max. number of function evaluations
# per line search
opt.verbosity = 1# verbosity: specifies how much
# information is displayed during
# policy learning. Options: 0-3


# 8. Plotting verbosity
plotting.verbosity = 0# 0: no plots
# 1: some plots
# 2: all plots

# 9. Some initializations
x = mcat([])
y = mcat([]); print y

fantasy.mean = cell(1, N)
fantasy.std = cell(1, N); print fantasy.std

realCost = cell(1, N)
M = cell(N, 1); print M
Sigma = cell(N, 1)
