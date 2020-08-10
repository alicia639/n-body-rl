"""
%% cartPole_learn.m
% *Summary:* Script to learn a controller for the cart-pole swingup
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-27
%
%% High-Level Steps
% # Load parameters
% # Create J initial trajectories by applying random controls
% # Controlled learning (train dynamics model, policy learning, policy
% application)

%% Code
"""

# 1. Initialization
clear(mstring('all'), mstring('close'), mstring('all'))
settings_cp()# load scenario-specific settings
basename = mstring('cartPole_')# filename used for saving data

# 2. Initial J random rollouts
for jj in mslice[1:J]:
    [xx, yy, realCost(jj), latent(jj)] = rollout(gaussian(mu0, S0), struct(mstring('maxU'), policy.maxU), H, plant, cost)
    x = mcat([x, OMPCSEMI, xx])
    y = mcat([y, OMPCSEMI, yy]); print y
    # augment training sets for dynamics model
    if plotting.verbosity > 0:    # visualization of trajectory
        if not ishandle(1):
            figure(1)
        else:
            set(0, mstring('CurrentFigure'), 1)
        end
        clf(1)

        draw_rollout_cp()
    end

end

mu0Sim(odei, mslice[:]).lvalue = mu0
S0Sim(odei, odei).lvalue = S0; print S0Sim

mu0Sim = mu0Sim(dyno)
S0Sim = S0Sim(dyno, dyno); print S0Sim


# 3. Controlled learning (N iterations)
for j in mslice[1:N]:
    trainDynModel()# train (GP) dynamics model
    learnPolicy()# learn policy
    applyController()# apply controller to system
    disp(mcat([mstring('controlled trial # '), num2str(j)]))
    if plotting.verbosity > 0:    # visualization of trajectory
        if not ishandle(1):
            figure(1)
        else:
            set(0, mstring('CurrentFigure'), 1)
        end
        clf(1)

        draw_rollout_cp()
    end
end
