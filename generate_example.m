%% generate_example - generate illustrative example
% Description:
%   Generates an illustrative two-resource allocation example. Used for the
%   generation of the illustrative example in [1].
% Outputs:
%   - network.mat (n, M, u_min, u_bar, u_max, weight_min, weight_max,
%   fw, P_go, T, l, C, d0, w_star, l_star, C_star)
% Assumptions and limitations:
%   - Sensitivity distribution is uniform (easy to extend to a general 
%   distribution, though)
% Other m-files required: none
% MAT-files required: none
% Toolboxes required: 
%   - YALMIP [2]
% Authors:
%   Leonardo Pedroso, Andrea Agazzi, W.P.M.H. (Maurice) Heemels, Mauro
%   Salazar
% Revision history:
%   01/03/2024 - Leonardo Pedroso
%       * Initial implementation
% References: 
%   [1] L. Pedroso, A. Agazzi, W.P.M.H. Heemels and M. Salazar, "Fair 
%   Artificial Currency Incentives in Repeated Weighted Congestion Games:
%   Equity vs. Equality", 2024. (submitted)
%   [2] Lofberg, Johan. "YALMIP: A toolbox for modeling and optimization in
%   MATLAB." In 2004 IEEE international conference on robotics and 
%   automation, pp. 284-289. IEEE, 2004.

%% Initialization
clear;
rng(1); % Reproducibility

%% Parameters
% Number of resources
n = 2; 
% Population
M = 1000;
% Urgency is an uniform distribution
u_min = 0;
u_bar = 1;
u_max = 2;
% Weight is a truncated normal distribution
weight_min = 1/2;
weight_bar = 1;
weight_max = 3/2;
weight_std = 0.15;
fw = @(w)  normpdf(w,weight_bar,weight_std).*(w>=weight_min).*(w<=weight_max)/...
    (normcdf(weight_max,weight_bar,weight_std)-normcdf(weight_min,weight_bar,weight_std));
% Behaviour
P_home = 0.05; % probability staying at home
P_go = 1 - P_home; % probability of travelling 
T = 4; % individual decision window
% Latency function (BPR)
alpha = 0.15;
beta = 4;
d0 = [1; 2];
kappa = [1/2; 2/3];
l = @(w) d0.*(1+alpha.*(w./kappa).^beta);
% Cost function
C = @(w) l(w)'*w;

%% Compute system optimum
% YALMIP to parse 
% Variables
w_star = sdpvar(n,1);
% Constraints
constr = [ones(1,n)*w_star == P_go;...
          0<= w_star <= 1];
opts = sdpsettings('verbose',0);
% Optimize
optimize(constr,C(w_star),opts);
% Get value
w_star = value(w_star);
% Compute ordering
% l_1(w_1^star) < ... < l_n(w_n^star)
[~,srt_idx] = sort(l(w_star));
d0 = d0(srt_idx);
kappa = kappa(srt_idx);
l = @(w) d0.*(1+alpha.*(w./kappa).^beta);
C = @(w) l(w)'*w;
w_star = w_star(srt_idx);

% SO
fprintf("System opt. flows (x_star):\t");
fprintf("%g ",w_star);
fprintf("\n");
% Latency corresponding to socially social optimum flows
l_star = l(w_star);
fprintf("System opt. disc. (l_star):\t");
fprintf("%g ",l_star);
fprintf("\n");
% Cost corresponding to socially social optimum flows 
C_star = C(w_star);
fprintf("System opt. cost:\t\t");
fprintf("%g ",C_star);
fprintf("\n");

%% Save illustrative network
save('example.mat','n','M','u_min','u_bar','u_max','weight_min',...
    'weight_max','fw','P_go','T','l','C','d0','w_star','l_star','C_star');
clear;