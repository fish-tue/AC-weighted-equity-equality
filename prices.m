%% prices - design AC prices for equity and for equality
% Description:
%   Design AC prices for equity and equality according to [1]
% Outputs:
%   prices.mat (prices design for equality and equity and expected fairness
%   metrics at the limit)
% Assumptions and limitations: none
% Other m-files required: none
% MAT-files required: example.mat
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

%% Load example 
load('example.mat','n','M','u_min','u_bar','u_max','weight_min',...
    'weight_max','fw','P_go','T','l','C','d0','w_star','l_star','C_star');

%% Parameters
S = 10;

%% Equity
% Prices
pEqt_cte = zeros(2,1);
pEqt_cte(1) = S*(w_star(2)/w_star(1));
pEqt_cte(2) = -S;
fprintf("Prices (Equity):\t\tp1 = %g \t p2 = %g\n",pEqt_cte(1),pEqt_cte(2));
% Expected unfairness metrics for equity design
L_Eqt_w = @(w) ones(size(w)).*(1/P_go)*(w_star(1)*l_star(1) + w_star(2)*l_star(2));
L_Eqt_per_weight_w  = @(w) L_Eqt_w(w)./w;
E_L_Eqt_w = integral(@(w) fw(w).*L_Eqt_w(w),0,inf);
E_L_Eqt_per_weight_w = integral(@(w) fw(w).*L_Eqt_per_weight_w(w),0,inf);
InEqt_Eqt = sqrt(integral(@(w) fw(w).*(L_Eqt_w(w)-E_L_Eqt_w).^2,0,inf));
InEql_Eqt = sqrt(integral(@(w) fw(w).*(L_Eqt_per_weight_w(w)-E_L_Eqt_per_weight_w).^2,0,inf));
pEqt = @(w) pEqt_cte;
fprintf("InEqt (Design for InEqt):\t%g\n",InEqt_Eqt);
fprintf("InEql (Design for InEqt):\t%g\n",InEql_Eqt);


%% Equality
% Parameters (s = w/theta)
ksi = P_go*l_star(2)/(l_star(2)-l_star(1));
sw_min = 1-(P_go-w_star(1))/(ksi-w_star(1));
sw_max = 1+(w_star(1))/(ksi-w_star(1));
% Synthesis (s = w/theta)
n1_s = @(s) min(max(-s*(ksi-w_star(1))+ksi,0),P_go);
expected_x1_pop = @(theta) integral(@(w) fw(w).*n1_s(w/theta).*w,0,+inf)/integral(@(w) fw(w).*w,0,+inf);
theta_star = fsolve(@(theta) w_star(1)-expected_x1_pop(theta),1,optimoptions('fsolve','Display','off'));
pEql = @(w) p_fnc(w,theta_star,P_go,n1_s,w_star,S);
% Expected unfairness metrics for equality design
L_Eql_w = @(w) (1/P_go)*(n1_s(w./theta_star)*l_star(1) + (P_go-n1_s(w./theta_star))*l_star(2));
L_Eql_per_weight_w  = @(w) L_Eql_w(w)./w;
E_L_Eql_w = integral(@(w) fw(w).*L_Eql_w(w),0,inf);
E_L_Eql_per_weight_w = integral(@(w) fw(w).*L_Eql_per_weight_w(w),0,inf);
InEqt_Eql = sqrt(integral(@(w) fw(w).*(L_Eql_w(w)-E_L_Eql_w).^2,0,inf));
InEql_Eql = sqrt(integral(@(w) fw(w).*(L_Eql_per_weight_w(w)-E_L_Eql_per_weight_w).^2,0,inf));
fprintf("InEqt (Design for InEql):\t%g\n",InEqt_Eql);
fprintf("InEql (Design for InEql):\t%g\n",InEql_Eql);

%% Plots

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

% Plot n_1(w,\theta^\star) and f_W(w)
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex');
N_points = 1000;
w_points = weight_min:((weight_max-weight_min)/(N_points -1)):weight_max;
n_aux = zeros(2,N_points);
for i = 1:N_points 
    n_aux(1,i) = n1_s(w_points(i)/theta_star);
    n_aux(2,i) = fw(w_points(i));
end
yyaxis left
ylim([0 1]);
plot(w_points,n_aux(1,:),'LineWidth',2,'Color',color.blue);
plot([theta_star theta_star],[0 1],'--','LineWidth',2,'Color','k');
yyaxis right
plot(w_points,n_aux(2,:),'LineWidth',2,'Color',color.red);
legend({'$n_1(w,\theta^\star)$','$\theta_\star$','$f_W(w)$'},'Interpreter','latex','Location','northeast');
xlabel('$w$','Interpreter','latex');
% Save figure to .fig and .svg formats
savefig('./figures/n1.fig');
set(gcf,'renderer','Painters');
exportgraphics(gcf,'figures/n1.png','Resolution',300);
hold off;

% Plot prices
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex');
N_points = 1000;
w_points = weight_min:((weight_max-weight_min)/(N_points -1)):weight_max;
p_aux = zeros(2,N_points);
for i = 1:N_points 
    p_aux(:,i) = pEql(w_points(i));
end
plot(w_points,p_aux(1,:),'LineWidth',2,'Color',color.orange);
plot(w_points,p_aux(2,:),'LineWidth',2,'Color',color.green);
plot([theta_star theta_star],[-S S],'--','LineWidth',2,'Color','k');
plot([weight_min weight_max],pEqt_cte(1)*ones(1,2),'--','LineWidth',2,'Color',color.orange);
plot([weight_min weight_max],pEqt_cte(2)*ones(1,2),'--','LineWidth',2,'Color',color.green);
legend({'$\mathbf{p}_1^{\mathrm{Eql}}$', '$\mathbf{p}_2^{\mathrm{Eql}}$','$\theta_\star$','$\mathbf{p}_1^{\mathrm{Eqt}}$', '$\mathbf{p}_2^{\mathrm{Eqt}}$'},'Interpreter','latex','Location','northeast','NumColumns',2);
xlabel('$w$','Interpreter','latex');
% Save figure to .fig and .png formats
savefig('./figures/prices.fig');
set(gcf,'renderer','Painters');
exportgraphics(gcf,'figures/prices.png','Resolution',300);
hold off;

%% Save prices
save('prices.mat','pEqt','pEql','InEqt_Eql','InEql_Eql','InEqt_Eqt',...
    'InEql_Eqt','S');
clear;

%% Auxiliary functions
function p = p_fnc(w,theta,P_go,n1_s,w_star,S)   
    if w/theta>1
        p = S*(w_star(2)/w_star(1))*[1 -n1_s(w/theta)/(P_go-n1_s(w/theta))]';
    else
        p = S*[(P_go-n1_s(w/theta))/n1_s(w/theta) -1]';
    end
end