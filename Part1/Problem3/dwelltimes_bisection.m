% -------------------------------------------------------------------------
% 4CM20 Hybrid Systems and Control 2018-2019
% Dwell time optimization by bisection
% -------------------------------------------------------------------------

% This file allows the optimization of the dwell time of a given switched
% system. Since the parameter that is to be optimized appears in a non-linear
% manner in the matrix inequality, the problem to be solved is no longer an
% an LMI. Generally, the solvers that were supplied will not allow you to
% solve this problem directly. However using bisection methods and
% iteratively solving the LMI problem for a fixed parameter leads to an
% optimized result where the optimization problem can be formulated as an
% LMI for each iteration.
%
% Code written by: Mark Haring, 2-12-2008
%
% Modifications:
%       v1: David Rijlaarsdam, 03-12-2008
%       v2: David Rijlaarsdam, 30-11-2009
%       v3: Niek Borgers,      09-11-2016
%       v4: Daniel Deenen,     05-12-2017

%% Decalarations
clear all
close all
clc

% System matrices
A1 = [-1 10; -100 -1];
A2 = [-1 100; -10 -1];

% Parameters
alpha = 0.01;       % \alpha of exercise
beta = 100;         % \beta of exercise
I = eye(2);

% Loop parameters
lb = 1;             % initial value lower search bound
ub = 10;            % initial value upper search bound
dsigma = 1e-8;      % maximal difference in sigma between to different iteration steps
ii = 100;           % maximal number of iteration steps

% Set initial loop values
lbf = NaN;      % lower bound feasible? NaN: not determined yet, 0: no, 1: yes
ubf = NaN;      % upper bound feasible? NaN: not determined yet, 0: no, 1: yes

sigmal = [];    % list to store sigmas of iteration loop

%% Iteration loop

if lb >= ub     % lower bound is larger than upper bound
    disp('Lower bound larger than upper bound. Change values.')
elseif lb <=0   % lower bound is negative
    disp('Sigma must be positive. Change lower bound.')
else
    for i = 1:ii
        % Determine feasibility of upper and lower bound
        if isnan(lbf)
            sigma = lb;      % determine feasibility lower bound
        elseif isnan(ubf)
            sigma = ub;     % determine feasibility upper bound
            
            % Expand range if both values are feasibile or not feasible
        elseif lbf == 0 && ubf == 0  	% both values are not feasible
            sigma = lb/2;                    % try lower value, reducing the lower bound by a factor 2
        elseif lbf == 1 && ubf == 1     % both values are feasible
            sigma = 2*ub;                  % try higher value, by doubling the interval length and keeping the lower bound fixed
            
            % Lower bound is feasible, upper bound is not feasible
        else
            sigma = (ub-lb)/2 + lb;      % bisecting the interval
        end
        
        % LMIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %               FORMULATE THE LMI PROBLEM HERE
        %
        % Defining the LMI problem and solve using SeDuMi:
        % Use notation P1var, P2var for the matrices to be solved and L for
        % your final set of constraints to follow the syntax of this file.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        options = sdpsettings('solver','sedumi');
        P1var = sdpvar(2,2);
        P2var = sdpvar(2,2);
        Lf1 = [A1'*P1var+P1var*A1+sigma*P1var <= 0];
        Lf2 = [A2'*P2var+P2var*A2+sigma*P2var <= 0];
        cP1 = [P1var >= alpha*I, P1var <= beta*I];
        cP2 = [P2var >= alpha*I, P2var <= beta*I];
        L = Lf1 + Lf2 +cP1 +cP2;
        diagnostics = optimize(L,[],options); % solve the LMI problem
        
        % Check result
        % pres : Primal constraint residuals
        % dres : Dual constraint residuals
        [pres,~] = check(L);
        
        % Results initial feasibility check
        if sigma == lb
            lbf = min(pres) > 0 && diagnostics.problem == 0;    % lower bound feasible? 0: no, 1: yes
        elseif sigma == ub
            ubf = min(pres) > 0 && diagnostics.problem == 0;    % upper bound feasible? 0: no, 1: yes
            
            % Resuls feasibility check if range is broadened
        elseif sigma > ub
            lb = ub;                % set new lower bound
            ub = sigma;             % set new upper bound
            ubf = min(pres) > 0 && diagnostics.problem == 0;    % determine feasibility upper bound
        elseif sigma < lb
            ub = lb;                % set new upper bound
            lb = sigma;             % set new upper bound
            lbf = min(pres) > 0 && diagnostics.problem == 0;    % determine feasibility lower bound
            
            % Determine feasibility after bisection
        else
            if min(pres) > 0 && diagnostics.problem == 0
                lb = sigma;     % new lower bound
            else
                ub = sigma;     % new upper bound
            end
            
            % Check convergence criterion?
            if abs(sigma-sigmal(i-1)) < dsigma
                break;              % break loop if accurate enough
            end
        end
        
        % Add sigma to history list
        sigmal(i) = sigma;
        
        % Store feasible solution
        if min(pres) > 0 && diagnostics.problem == 0
            feasible.P1var = P1var;
            feasible.P2var = P2var;
            feasible.sigma = sigma;
        end
        
    end
    
    % Plot iteration history
    figure(1)
    plot(sigmal,'ro-','LineWidth',2)
    xlabel('iteration step')
    ylabel('\sigma')
    title('Iteration results')
    grid on
    axis tight
    
    % Retrieve output
    P1    = value(feasible.P1var)
    P2    = value(feasible.P2var)
    sigma = feasible.sigma
end


%% Calculate Td
c1 = sqrt(max(eig(P1))/min(eig(P1)));
c2 = sqrt(max(eig(P2))/min(eig(P2)));
c = max(c1,c2);
lambda = sigma/2;
Td = log(c)/lambda;