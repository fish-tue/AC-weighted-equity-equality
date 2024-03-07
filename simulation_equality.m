%% simulation_equality - simulate AC mechanism designed for equality
% Description:
%   Simulation of the daily Nash equilibrium and the decisions of each
%   user, when the prices are designed for equality [1]
% Outputs:
%   simulation_eql.mat (simulation output)
% Assumptions and limitations:
%   - Sensitivity distribution is uniform (easy to extend to a general 
%   distribution, though)
% Other m-files required:
%   - n_arcs_individual.m
%   - gamma_n_arcs_individual.m
% MAT-files required: example.mat, prices.mat
% Toolboxes required: none
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

%% Initialization
clear;
rng(1); % Reproducibility

%% Load example and prices 
load('example.mat','n','M','u_min','u_bar','u_max','weight_min',...
    'weight_max','fw','P_go','T','l','C','d0','w_star','l_star','C_star');
load('prices.mat','pEql','InEqt_Eql','InEql_Eql','S');
p = pEql;

%% Simulation parameters
% Nash iterations
it_nash_max = 50;
it_nash_epsl = 1e-3;
% Simulation
T_sim = 500; 
k_init_max = 10*S;

%% Simulation
% Simulation - Initialization
k = zeros(M,T_sim+1); % AC 
weight = zeros(M,1); % weight
s = zeros(M,T_sim+1); % sensitivity
y = zeros(M,T_sim+1);
x = zeros(n,T_sim+1); % aggregate choices
it_nash = zeros(T_sim+1,1); % number of iterations for the nash equilibrium
% Weight initialization
% Inverse transform sampling of weight distribution
u_pop = rand(M,1);
opts = optimset('Display','off');
for j = 1:M
    fun = @(x) integral(@(x) fw(x),0,x)-u_pop(j);
    weight(j,1) = fsolve(fun,(weight_min+weight_max)/2,opts);
end
w_sum = sum(weight);
% AC initialization
k(:,1) = unifrnd(k_init_max/2,k_init_max,M,1); % uniform initial karma distribution
% Simulation loop
for t = 1:T_sim+1
    % Pick daily sensitivities
    s(:,t) = unifrnd(0,u_max,M,1);
    % Agents decisions
    y_go = binornd(1,P_go,M,1); % mask of travelling agents    
    % Init Nash equilibrium iterations
    if t > 1      
        x(:,t) = x(:,t-1);
    else
        x(:,t) = w_star;
    end
    % Until convergence of Nash equilibrium is reached
    I = eye(n);
    y(:,t) = ones(M,1);
    while true
        % Previous iteration's flows
        x_prev = x(:,t);
        % Agents decisions
        for i = randperm(M)
            if ~y_go(i)
                y(i,t) = 0;
                continue; 
            end % non-traveling agent
            y(i,t) = n_arcs_individual(...
               l(x(:,t)+weight(i)*(-I(:,y(i,t))/M+ones(n,1)/M)),T,p(weight(i)),k(i,t),0,s(i,t),u_min,u_bar,u_max);
            % Aggregate behaviour
            for j = 1:n
                x(j,t) = sum(weight(y(:,t)==j))/w_sum;
            end
        end
        % Catch infeasibility
        if sum(isnan(y(:,t)))
            error("Caught infeasibility!");
        end
        % Catch Nash equilibrium iterations not converging
        it_nash(t) = it_nash(t)+1;
        if norm(x_prev-x(:,t))<1e-3/M, break; end
        if it_nash(t) > it_nash_max
            warning("Nash iterations did not converge.");
            break; 
        end
    end
    % AC dynamics
    if t < T_sim+1
        for i = 1:M
            if ~y_go(i)
                k(i,t+1) = k(i,t);
            else
                p_aux = p(weight(i));
                k(i,t+1) = k(i,t)-p_aux(y(i,t));
            end
        end
    end  
end

%% Plot simulation evolution

% Colorblind-safe color palette (Wong, B. Points of view: Color blindness.
% Nat Methods 8,441 (2011). https://doi.org/10.1038/nmeth.1618)
color.black = [0 0 0]/255;
color.orange = [230 159 0]/255;
color.cyan = [86 180 233]/255;
color.green = [0 158 115]/255;
color.yellow = [240 228 66]/255;
color.blue = [0 114 178]/255;
color.red = [213 94 0]/255;
color.pink = [204 121 167]/255;

%% Decision Evolution
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
for j = n:-1:1
    aux_x = (0:T_sim)';
    aux_y = sum(x(1:j,1:T_sim+1),1)';
    [aux_x,aux_y] = stairs_vector(aux_x,aux_y);
    if j==1
        area(aux_x,aux_y,'LineWidth',2,'FaceColor',color.blue);
    elseif j ==2
        area(aux_x,aux_y,'LineWidth',2,'FaceColor',color.red);
    end
end
for j = n:-1:1
    if j == 1
        plot([0 T_sim],[sum(w_star(1:j)) sum(w_star(1:j))],'--','Color','black','LineWidth',3);
    else
        plot([0 T_sim],[sum(w_star(1:j)) sum(w_star(1:j))],'--','Color','black','LineWidth',3,'HandleVisibility','off');
    end
end
legend({'$\mathbf{w}^{A_t}_2$','$\mathbf{w}^{A_t}_1$','$\,\mathbf{w}^\star$'},'Location','northeast','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylim([0 1]);
xlim ([0 100]);
% Save figure to .fig and .png formats
% savefig('./figures/decision_Eql.fig');
% set(gcf,'renderer','Painters');
% exportgraphics(gcf,'./figures/decision_Eql.png','Resolution',300);
% exportgraphics(gcf,'./figures/decision_Eql.pdf');
hold off;

% AC level
k_mean = mean(k(:,1:T_sim+1),1);
k_max = max(k(:,1:T_sim+1));
k_min = min(k(:,1:T_sim+1));
k_var = sqrt(var(k(:,1:T_sim+1)));
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
aux_x2 = [0:T_sim fliplr(0:T_sim)]';
aux_y2 = [k_max fliplr(k_min)]';
[aux_x2,aux_y2] = stairs_vector(aux_x2,aux_y2);
fill(aux_x2,aux_y2*1e-2,'k',...
    'LineWidth',2,'FaceColor',color.blue,'FaceAlpha',0.2,'EdgeAlpha',0);
aux_x1 = [0:T_sim fliplr(0:T_sim)]';
aux_y1 = [k_mean+k_var fliplr(k_mean-k_var)]';
aux_y1(aux_y1<0) = 0;
[aux_x1,aux_y1] = stairs_vector(aux_x1,aux_y1);
fill(aux_x1,aux_y1*1e-2,'k',...
    'LineWidth',2,'FaceColor',color.blue,'FaceAlpha',0.4,'EdgeAlpha',0);
stairs(0:T_sim,k_mean*1e-2,'LineWidth',2,'Color','black');
stairs(aux_x1,aux_y1*1e-2,'LineWidth',1,'Color',[0.4 0.4 0.4]);
stairs(aux_x2,aux_y2*1e-2,'LineWidth',1,'Color',[0.4 0.4 0.4]);
legend({' $\max$/$\min$ $\{K_t\}$',' $\mathrm{E}[K_t]\pm\sqrt{\mathrm{Var}[K_t]}$',' $\mathrm{E}[K_t]$'},...
    'Location','northeast','Interpreter','latex');
ylabel('AC level $\times 10^{-2}$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
xlim ([0 100]);
ylim ([0 1.5]);
% Save figure to .fig and .png formats
% savefig('./figures/AC_Eql.fig');
% set(gcf,'renderer','Painters');
% exportgraphics(gcf,'./figures/AC_Eql.png','Resolution',300);
hold off;

% System's cost
cost_soc = zeros(T_sim+1,1);
cost_soc_rel_opt = zeros(T_sim+1,1);
cost_soc_opt = C(w_star);
for t = 1:T_sim+1
    cost_soc(t) = C(x(:,t));
    cost_soc_rel_opt(t) = (cost_soc(t)-cost_soc_opt)/cost_soc_opt;
end
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
stairs(0:T_sim,cost_soc_rel_opt*100,'LineWidth',2.5,'Color',color.black);
ylabel('$\Delta$ Societal cost $(\%)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylim([-3 35]);
hold off;
xlim ([0 100]);
% Save figure to .fig and .png formats
% savefig('./figures/cost_Eql.fig');
% set(gcf,'renderer','Painters');
% exportgraphics(gcf,'./figures/cost_Eql.png','Resolution',300);
hold off;

% Sensitivity to latency
delta_s = 100*mean(s(:,1:T_sim+1)-u_bar)/u_bar;
delta_d = zeros(T_sim+1,1);
I = eye(n);
for t = 1:T_sim+1
    aux1 = 0;
    aux2 = 0;
    for i = 1:M
        if ~y(i,t), continue; end
        aux1 = aux1 + s(i,t)*l(x(:,t))'*I(:,y(i,t)) - u_bar.*l(x(:,t))'*I(:,y(i,t));
        aux2 = aux2 + u_bar.*l(x(:,t))'*I(:,y(i,t));
    end
    delta_d(t) = 100*aux1/aux2;
end
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex');
stairs(0:T_sim,delta_s,':','LineWidth',2.5,'Color','black');
stairs(0:T_sim,delta_d,'LineWidth',2.5,'Color','black');
legend({'$\Delta \bar{u}(t)$','$\Delta \bar{l}(t)$'},'Interpreter','latex','Location','southwest');
ylabel('Urgency, latency [$\%$]','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
hold off;
xlim ([0 100]);
ylim ([-15 5]);
% Save figure to .fig and .png formats
% savefig('./figures/sensitivity_Eql.fig');
% set(gcf,'renderer','Painters');
% exportgraphics(gcf,'./figures/sensitivity_Eql.png','Resolution',300);
hold off;

% Equity vs Equality
unf_l_equity = zeros(M,T_sim+1);
unf_l_equality = zeros(M,T_sim+1);
E_L = zeros(T_sim+1,1);
E_L_per_weight = zeros(T_sim+1,1);
InEqt = zeros(T_sim+1,1);
InEql = zeros(T_sim+1,1);
part_times_n = zeros(M,T_sim+1);
for i = 1:M
    for t = 1:T_sim+1
        part_times_n(i,t) = sum(y(i,1:t)~=0);
        if t ~= 1
            if y(i,t)
                unf_l_equity(i,t) = unf_l_equity(i,t-1)*((part_times_n(i,t)-1)/part_times_n(i,t))+l(x(:,t))'*I(:,y(i,t))*(1/part_times_n(i,t));
                unf_l_equality(i,t) = unf_l_equality(i,t-1)*((part_times_n(i,t)-1)/part_times_n(i,t))+(1/weight(i))*l(x(:,t))'*I(:,y(i,t))*(1/part_times_n(i,t));
            else
                unf_l_equity(i,t) = unf_l_equity(i,t-1);
                unf_l_equality(i,t) = unf_l_equality(i,t-1);
            end
        else
            if y(i,t)
                unf_l_equity(i,t) = l(x(:,t))'*I(:,y(i,t));
                unf_l_equality(i,t) = (1/weight(i))*l(x(:,t))'*I(:,y(i,t));
            else
                unf_l_equity(i,t) = 0;
                unf_l_equality(i,t) = 0;
            end
        end
    end
end
for t = 1:T_sim+1
    E_L(t) = mean(unf_l_equity(:,t));
    E_L_per_weight(t) = mean(unf_l_equality(:,t));
    InEqt(t) = sqrt(var(unf_l_equity(:,t)));
    InEql(t) = sqrt(var(unf_l_equality(:,t)));
end
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
aux_x1 = [0:T_sim fliplr(0:T_sim)]';
aux_y1 = [E_L'+InEqt' fliplr(E_L'-InEqt')]';
aux_y1(aux_y1<0) = 0;
[aux_x1,aux_y1] = stairs_vector(aux_x1,aux_y1);
fill(aux_x1,aux_y1,'k',...
    'LineWidth',2,'FaceColor',color.green,'FaceAlpha',0.5,'EdgeAlpha',0);
aux_x2 = [0:T_sim fliplr(0:T_sim)]';
aux_y2 = [E_L_per_weight'+InEql' fliplr(E_L_per_weight'-InEql')]';
aux_y2(aux_y2<0) = 0;
[aux_x2,aux_y2] = stairs_vector(aux_x2,aux_y2);
fill(aux_x2,aux_y2,'k',...
    'LineWidth',2,'FaceColor',color.orange,'FaceAlpha',0.5,'EdgeAlpha',0);
stairs(0:T_sim,E_L,'LineWidth',2.5,'Color',color.green);
stairs(0:T_sim,E_L_per_weight,'LineWidth',2.5,'Color',color.orange);
stairs(aux_x1,aux_y1,'LineWidth',1,'Color',color.green);
stairs(aux_x2,aux_y2,'LineWidth',1,'Color',color.orange);
legend({'$\mathrm{E}[L_t] \pm \sqrt{\mathrm{Var}[L_t]}$', '$\mathrm{E}[L_t/W] \pm \sqrt{\mathrm{Var}[L_t/W]}$'},'Interpreter','latex','Location','northeast');
xlabel('$t$','Interpreter','latex');
ylim([1.3 2.5]);
% Save figure to .fig and .png formats
% savefig('./figures/L_Eql.fig');
% set(gcf,'renderer','Painters');
% exportgraphics(gcf,'./figures/L_Eql.png','Resolution',300);
% exportgraphics(gcf,'./figures/L_Eql.pdf');
hold off;

figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex');
% set(gca,'XScale','log')
stairs(0:T_sim,InEqt,'LineWidth',2.5,'Color',color.green);
stairs(0:T_sim,InEql,'LineWidth',2.5,'Color',color.orange);
plot([0 1 T_sim],InEqt_Eql*ones(1,3),'--','LineWidth',2,'Color',color.green);
plot([0 1 T_sim],InEql_Eql*ones(1,3),'--','LineWidth',2,'Color',color.orange);
legend({'$\mathrm{InEqt}_t$', '$\mathrm{InEql}_t$'},'Interpreter','latex','Location','northeast');
xlabel('$t$','Interpreter','latex');
ylim([0 0.6]);
% Save figure to .fig and .png formats
% savefig('./figures/Unf_Eql.fig');
% set(gcf,'renderer','Painters');
% exportgraphics(gcf,'./figures/Unf_Eql.png','Resolution',300);
hold off;

%% Save simulation results
save('simulation_eql.mat')
clear;

%% Auxiliary functions
% To use the equivalent of stairs with area and fill functions
function [x,y] = stairs_vector(x,y)
    x = [x';x'];
    y = [y';y'];
    y = y(:);
    x = x([2:end end])';
end