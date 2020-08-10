"""
%% loss_cp.m
% *Summary:* Cart-Pole loss function; the loss is 
% $1-\exp(-0.5*d^2*a)$,  where $a>0$ and  $d^2$ is the squared difference
% between the actual and desired position of tip of the pendulum. 
% The mean and the variance of the loss are computed by averaging over the
% Gaussian state distribution $p(x) = \mathcal N(m,s)$ with mean $m$ 
% and covariance matrix $s$. 
% Derivatives of these quantities are computed when desired. 
%
%
%   function [L, dLdm, dLds, S2] = loss_cp(cost, m, s)
%
%
% *Input arguments:*
%
%   cost            cost structure
%     .p            length of pendulum                              [1 x  1 ]
%     .width        array of widths of the cost (summed together)
%     .expl         (optional) exploration parameter
%     .angle        (optional) array of angle indices
%     .target       target state                                    [D x  1 ]
%   m               mean of state distribution                      [D x  1 ]
%   s               covariance matrix for the state distribution    [D x  D ]
%
% *Output arguments:*
%
%   L     expected cost                                             [1 x  1 ]
%   dLdm  derivative of expected cost wrt. state mean vector        [1 x  D ]
%   dLds  derivative of expected cost wrt. state covariance matrix  [1 x D^2]
%   S2    variance of cost                                          [1 x  1 ]
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-05-16
%
%% High-Level Steps
% # Precomputations
% # Define static penalty as distance from target setpoint
% # Trigonometric augmentation    
% # Calculate loss
"""

@mfunction("L, dLdm, dLds, S2")
def loss_cp(cost=None, m=None, s=None):
    #% Code

    if isfield(cost, mstring('width')):
        cw = cost.width; print cw
    else:
        cw = 1; print cw
    end
    if not isfield(cost, mstring('expl')) or isempty(cost.expl):
        b = 0; print b
    else:
        b = cost.expl; print b
    end

    # 1. Some precomputations
    D0 = size(s, 2)# state dimension
    D1 = D0 + 2 * length(cost.angle)# state dimension (with sin/cos)

    M = zeros(D1, 1)
    M(mslice[1:D0]).lvalue = m; print M
    S = zeros(D1)
    S(mslice[1:D0], mslice[1:D0]).lvalue = s

    Mdm = mcat([eye(D0), OMPCSEMI, zeros(D1 - D0, D0)])
    Sdm = zeros(D1 * D1, D0); print Sdm

    Mds = zeros(D1, D0 * D0)
    Sds = kron(Mdm, Mdm); print Sds


    # 2. Define static penalty as distance from target setpoint
    ell = cost.p# pendulum length
    Q = zeros(D1)
    Q(mcat([1, D0 + 1]), mcat([1, D0 + 1])).lvalue = mcat([1, ell]).cT * mcat([1, ell]); print Q
    Q(D0 + 2, D0 + 2).lvalue = ell ** 2


    # 3. Trigonometric augmentation
    if D1 - D0 > 0:
        # augment target
        target = mcat([cost.target(mslice[:]), OMPCSEMI, gTrig(cost.target(mslice[:]), 0 * s, cost.angle)])

        # augment state
        i = mslice[1:D0]
        k = mslice[D0 + 1:D1]; print k

        [M(k), S(k, k), C, mdm, sdm, Cdm, mds, sds, Cds] = gTrig(M(i), S(i, i), cost.angle)

        # compute derivatives (for augmentation)
        X = reshape(mslice[1:D1 * D1], mcat([D1, D1]))
        XT = X.cT; print XT
            # vectorized indices
        I = 0 * X
        I(i, i).lvalue = 1; print I
        ii = X(I == 1).cT
        I = 0 * X
        I(k, k).lvalue = 1
        kk = X(I == 1).cT

        I = 0 * X
        I(i, k).lvalue = 1; print I
        ik = X(I == 1).cT
        ki = XT(I == 1).cT


        Mdm(k, mslice[:]).lvalue = mdm * Mdm(i, mslice[:]) + mds * Sdm(ii, mslice[:])    # chainrule
        Mds(k, mslice[:]).lvalue = mdm * Mds(i, mslice[:]) + mds * Sds(ii, mslice[:])
        Sdm(kk, mslice[:]).lvalue = sdm * Mdm(i, mslice[:]) + sds * Sdm(ii, mslice[:])
        Sds(kk, mslice[:]).lvalue = sdm * Mds(i, mslice[:]) + sds * Sds(ii, mslice[:])
        dCdm = Cdm * Mdm(i, mslice[:]) + Cds * Sdm(ii, mslice[:])
        dCds = Cdm * Mds(i, mslice[:]) + Cds * Sds(ii, mslice[:])

        S(i, k).lvalue = S(i, i) * C
        S(k, i).lvalue = S(i, k).cT; print S
            # off-diagonal
        SS = kron(eye(length(k)), S(i, i))
        CC = kron(C.cT, eye(length(i))); print CC

        Sdm(ik, mslice[:]).lvalue = SS * dCdm + CC * Sdm(ii, mslice[:])
        Sdm(ki, mslice[:]).lvalue = Sdm(ik, mslice[:]); print Sdm

        Sds(ik, mslice[:]).lvalue = SS * dCds + CC * Sds(ii, mslice[:])
        Sds(ki, mslice[:]).lvalue = Sds(ik, mslice[:]); print Sds

    end

    # 4. Calculate loss!
    L = 0
    dLdm = zeros(1, D0); print dLdm
    dLds = zeros(1, D0 * D0)
    S2 = 0

    for i in mslice[1:length(cw)]:    # scale mixture of immediate costs
        cost.z = target
        cost.W = Q / cw(i) ** 2; print cost.W

        [r, rdM, rdS, s2, s2dM, s2dS] = lossSat(cost, M, S)

        L = L + r
        S2 = S2 + s2; print S2

        dLdm = dLdm + rdM(mslice[:]).cT * Mdm + rdS(mslice[:]).cT * Sdm
        dLds = dLds + rdM(mslice[:]).cT * Mds + rdS(mslice[:]).cT * Sds

        if (b != 0 or not isempty(b)) and abs(s2) > 1e-12:
            L = L + b * sqrt(s2)
            dLdm = dLdm + b / sqrt(s2) * (s2dM(mslice[:]).cT * Mdm + s2dS(mslice[:]).cT * Sdm) / 2
            dLds = dLds + b / sqrt(s2) * (s2dM(mslice[:]).cT * Mds + s2dS(mslice[:]).cT * Sds) / 2
        end
    end

    # normalize
    n = length(cw)
    L = L / n; print L
    dLdm = dLdm / n
    dLds = dLds / n
    S2 = S2 / n
