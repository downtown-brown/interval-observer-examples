% Author: Scott Brown (sab007@ucsd.edu)
% 2022

clear all
set(groot,'defaultAxesTickLabelInterpreter','latex');  
rng(1234)
for i = 1:4
    try
        clf(i)
    end
end
try
    clf(100)
end
theta = 0;
A = [cos(theta) -sin(theta) 1 0;sin(theta) cos(theta) 0 1;0 0 1 0; 0 0 0 1];
% A = randn(2)/2;
N = 3;
n = size(A,2);
m = [1 1 1];

L_s = [0 0 1;
       0 0 1;
       1 1 0];

neigh = cell(N,1);
for i = 1:N
    neigh{i} = find(L_s(i,:));
    neigh{i} = neigh{i}(neigh{i}~=i);
end
gamma = 1;
c{1} = [1 0 0 0];
c{2} = [0 1 0 0];
c{3} = [1 1 0 0];

for i = 1:N
        cvx_begin sdp
            variable L_c(n,m)
            variable T_c(n,n)
            variable G_c(n,m)
            variable E(n,n)
            A_res = T_c*A-L_c*c{i};
            minimize sum(sum(E))
            subject to
                T_c == eye(n) - G_c*c{i}
                -vec(E) <= vec(A_res)
                vec(E) >= vec(A_res)
                vec(E) >= 0
        cvx_end

        L{i} = L_c;
        T_{i} = T_c;
        N_{i} = G_c;

end

for i = 1:N
    L_p{i} = L{i}.*(L{i}>0);
    L_m{i} = L_p{i} - L{i};
    N_p{i} = N_{i}.*(N_{i}>0);
    N_m{i} = N_p{i} - N_{i};
    T_p{i} = T_{i}.*(T_{i}>0);
    T_m{i} = T_p{i} - T_{i};
    ALH{i} = T_{i}*A-L{i}*c{i};
    ALH_p{i} = ALH{i}.*(ALH{i}>0);
    ALH_m{i} = ALH_p{i} - ALH{i};
end

c_c = vertcat(c{:})
cvx_begin sdp
    variable L_c(n,3)
    variable T_c(n,n)
    variable N_c(n,3)
    variable E(n,n)
    A_res = T_c*A-L_c*c_c;
    minimize sum(sum(E))
    subject to
        T_c == eye(n) - N_c*c_c
        -vec(E) <= vec(A_res)
        vec(E) >= vec(A_res)
        vec(E) >= 0
cvx_end


L_p_c = L_c.*(L_c>0);
L_m_c = L_p_c - L_c;
N_p_c = N_c.*(N_c>0);
N_m_c = N_p_c - N_c;
T_p_c = T_c.*(T_c>0);
T_m_c = T_p_c - T_c;
ALH_c = T_c*A-L_c*c_c;

T = 16;

d_w = [1, 1, .2, .2]';
d_v = cell(N,1);
d_v{1} = 1;
d_v{2} = 2;
d_v{3} = 3;


%%
x_0_l = [-20 -15 -0.5 0]';
x_0_h = [10 25 2 3]';
dx_0 = x_0_h - x_0_l;

x(:,1) = [-14, 3, 1.7, 1.4];

x_h_c(:,1) = ones(n,1);
x_l_c(:,1) = zeros(n,1);

w = d_w.*[sin((1:T)*pi/5);
          cos((1:T)*pi/5);
          rand(1,T)-0.5;
          rand(1,T)-0.5;]
w_h = [1, 1, .1, .1]';
w_l = -[1, 1, .1, .1]';

% w(3,:) = cos((1:T)*0.05)/2
% w(4,:) = sin((1:T)*0.05)/2

for i = 1:N
    v{i} = d_v{i}.*(rand(m(i),T) - 0.5);
    v_h{i} = d_v{i}*0.5;
    v_l{i} = d_v{i}*-0.5;
end

v_c = vertcat(v{:});
v_h_c = vertcat(v_h{:});
v_l_c = vertcat(v_l{:});

x_h_c(:,1) = x_0_h;
x_l_c(:,1) = x_0_l;

for i = 1:N
    y{i}(:,1) = c{i}*x(:,1);% + v{i}(:,1);
    x_h{i}(:,1) = x_0_h;
    x_l{i}(:,1) = x_0_l;
    x_h_loc{i}(:,1) = x_0_h;
    x_l_loc{i}(:,1) = x_0_l;
    x_h_best{i}(:,1) = x_0_h;
    x_l_best{i}(:,1) = x_0_l;
end

for t = 1:T-1
    x(:,t+1) = A*x(:,t) + w(:,t);
    for i = 1:N
        y{i}(:,t+1) = c{i}*x(:,t+1) + v{i}(:,t+1);
    end

    y_c = vertcat(y{:});
    x_h_c(:,t+1) = ALH_c*x_h_c(:,t) + L_c*y_c(:,t) + N_c*y_c(:,t+1) + T_p_c*w_h - T_m_c*w_l + (L_p_c + N_p_c)*v_h_c - (L_m_c + N_m_c)*v_l_c;
    x_l_c(:,t+1) = ALH_c*x_l_c(:,t) + L_c*y_c(:,t) + N_c*y_c(:,t+1) + T_p_c*w_l - T_m_c*w_h + (L_p_c + N_p_c)*v_l_c - (L_m_c + N_m_c)*v_h_c;

    
    for i = 1:N
        x_h{i}(:,t+1) = ALH_p{i}*x_h{i}(:,t) - ALH_m{i}*x_l{i}(:,t) + L{i}*y{i}(:,t) + N_{i}*y{i}(:,t+1) + T_p{i}*w_h - T_m{i}*w_l + (L_p{i} + N_p{i})*v_h{i} - (L_m{i} + N_m{i})*v_l{i};
        x_l{i}(:,t+1) = ALH_p{i}*x_l{i}(:,t) - ALH_m{i}*x_h{i}(:,t) + L{i}*y{i}(:,t) + N_{i}*y{i}(:,t+1) + T_p{i}*w_l - T_m{i}*w_h + (L_p{i} + N_p{i})*v_l{i} - (L_m{i} + N_m{i})*v_h{i};
    
        x_h_loc{i}(:,t+1) = ALH_p{i}*x_h_loc{i}(:,t) - ALH_m{i}*x_l_loc{i}(:,t) + L{i}*y{i}(:,t) + N_{i}*y{i}(:,t+1) + T_p{i}*w_h - T_m{i}*w_l + (L_p{i} + N_p{i})*v_h{i} - (L_m{i} + N_m{i})*v_l{i};
        x_l_loc{i}(:,t+1) = ALH_p{i}*x_l_loc{i}(:,t) - ALH_m{i}*x_h_loc{i}(:,t) + L{i}*y{i}(:,t) + N_{i}*y{i}(:,t+1) + T_p{i}*w_l - T_m{i}*w_h + (L_p{i} + N_p{i})*v_l{i} - (L_m{i} + N_m{i})*v_h{i};

    end
    
    for i = 1:N
        x_h_best{i}(:,t+1) = x_h{i}(:,t+1);
        x_l_best{i}(:,t+1) = x_l{i}(:,t+1);
        for j = neigh{i}
            x_h_best{i}(:,t+1) = min(x_h_best{i}(:,t+1), x_h{j}(:,t+1));
            x_l_best{i}(:,t+1) = max(x_l_best{i}(:,t+1), x_l{j}(:,t+1));
        end
    end
    for i = 1:N
        x_h{i}(:,t+1) = x_h_best{i}(:,t+1);
        x_l{i}(:,t+1) = x_l_best{i}(:,t+1);
    end
end

e = [x_h{1}-x;-x_l{1}+x];

%%
ylblx = {'$x_1$','$x_2$','$x_3$','$x_4$'};
ylbld = {'$\delta_1 = \overline{x}_1-\underline{x}_1$','$\delta_2 = \overline{x}_2-\underline{x}_2$','$\delta_{3} = \overline{x}_{3}-\underline{x}_{3}$'};

limsx = {[-100, 100], [-44 67],[-15 17], [-14,14]};
limsd = {[0, 160], [0, 120],[0,40], [0,40]};
figure(1)
    tiledlayout(3,2,'TileSpacing','tight','Padding','tight')
    for kk = 1:3
        nexttile;
        plot(0:T-1,x(kk,:)','k','linewidth',1)
        hold on
        plot([NaN,NaN],'w')
        plot(0:T-1,x_h_loc{1}(kk,:)','rx--','MarkerSize',4,'linewidth',1)
        plot(0:T-1,x_l_loc{1}(kk,:)','bx--','MarkerSize',4,'linewidth',1)
        plot(0:T-1,x_h_loc{2}(kk,:)','r+-.','MarkerSize',4,'linewidth',1)
        plot(0:T-1,x_l_loc{2}(kk,:)','b+-.','MarkerSize',4,'linewidth',1)
        plot(0:T-1,x_h_loc{3}(kk,:)','r*--','MarkerSize',4,'linewidth',1)
        plot(0:T-1,x_l_loc{3}(kk,:)','b*--','MarkerSize',4,'linewidth',1)
        plot(0:T-1,x_h{1}(kk,:)','^--','MarkerSize',4,'linewidth',1,'color',[0 .5 0])
        plot(0:T-1,x_l{1}(kk,:)','v--','MarkerSize',4,'linewidth',1,'color',[0 .5 0])
        hold off
        ylabel(ylblx{kk}, 'interpreter','latex')
        ylim(limsx{kk})
        xlim([0,7])
        set(gca,'FontSize',10.5)
        if kk == 3
            xlabel('Time, $k$','interpreter','latex')
        end
        if kk == 1
        leg = legend('$x_i$','',...
            '$\overline{x}_i^1$','$\underline{x}_i^1$',...
            '$\overline{x}_i^2$','$\underline{x}_i^2$',...
            '$\overline{x}_i^3$','$\underline{x}_i^3$',...
            '$\overline{x}_i$', '$\underline{x}_i$',...
            'interpreter','latex',...
            'location','northoutside',...
            'Orientation','horizontal',...
            'numcolumns',5,...
            'FontSize',10)
            h = get(gca,'Children');
            set(gca,'Children',h([2,4,6,8,9,1,3,5,7,10]))
        leg.ItemTokenSize = [10,10];
        end
        set(gca,'FontSize',10)
        
        nexttile
        
        plot(0:T-1,x_h_loc{1}(kk,:)'- x_l_loc{1}(kk,:)','m-.','MarkerSize',4,'linewidth',1)
        hold on
        plot(0:T-1,x_h_loc{2}(kk,:)' - x_l_loc{2}(kk,:)','+--','color',[0 1/2 0],'MarkerSize',4,'linewidth',1)
        plot(0:T-1,x_h_loc{3}(kk,:)' - x_l_loc{3}(kk,:)','*--','color',[0 1/2 1/2],'MarkerSize',4,'linewidth',1)
        plot(0:T-1,x_h{1}(kk,:)' - x_l{1}(kk,:)','dk--','MarkerSize',4,'linewidth',1)
        hold off
        ylabel(ylbld{kk}, 'interpreter','latex')
        ylim(limsd{kk})
        xlim([0,7])
%         xlabel('$k$','interpreter','latex')
        
        if kk == 1
        leg = legend(...
            '$\delta_i^1$',...
            '$\delta_i^2$',...
            '$\delta_i^3$',...
            '$\delta_i$',...
            'interpreter','latex',...
            'location','northoutside',...
            'Orientation','horizontal',...
            'FontSize',10);
        leg.ItemTokenSize = [10,10];
        end
        set(gca,'FontSize',10)
    end
   
xlabel('Time, $k$','interpreter','latex')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4 5]);
saveas(gcf,sprintf('img/ex1.eps', kk),'epsc')
saveas(gcf,sprintf('img/ex1.png', kk))
saveas(gcf,sprintf('img/ex1.fig', kk))