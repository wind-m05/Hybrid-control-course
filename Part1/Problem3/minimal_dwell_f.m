%% Decalarations
clear all
close all
clc

% System matrices
A1 = [-7 ,4   ,6;
       8 ,-47 ,-60;
       0 ,36  ,45];
A2 = [7 , 16,  8;
      3 ,  5, -1;
      1 ,  6, 11];
% Parameters
alpha = 0.01;       % \alpha of exercise
beta = 100;         % \beta of exercise
I = eye(3);

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
        options = sdpsettings('solver','mosek','verbose',0);
%         P1var = sdpvar(3,3);
        P2var = sdpvar(3,3);
%         Lf1 = [A1'*P1var+P1var*A1+sigma*P1var <= 0];
        Lf2 = [A2'*P2var+P2var*A2+sigma*P2var >= 0];
%         cP1 = [P1var >= alpha*I, P1var <= beta*I];
        cP2 = [P2var <= -alpha*I, P2var >= -beta*I];
        L = Lf2 +cP2;
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
%             feasible.P1var = P1var;
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
%     P1    = value(feasible.P1var)
    P2    = value(feasible.P2var)
    sigma = feasible.sigma
end


%% Calculate Td
% mu12 = max(eig(P1))/min(eig(P2));
% mu21 = max(eig(P2))/min(eig(P1));
% mu = (mu12+mu21)/2;
% lambda0 = sigma;
% Td = log(mu)/lambda0;