% Author: Scott Brown (sab007@ucsd.edu)
% 2022
%
% This code requires CVX
%
% Run this script after running dist_interval_compr_theirs.m, for
% comparison

clear all
load theirs
dt_theirs = 0.01;

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

T = 30;
dt = 0.5;

N = 6;
n = 12;
m = 1;

A = kron(eye(N),[0,1;0,0]);
B = eye(n);

sys_c = ss(A,B,eye(n),zeros(n));
sys_d = c2d(sys_c,dt,'foh');

A = sys_d.A;
B = sys_d.B;

w = kron(ones(4,1),[sin(0:dt:T)/2; zeros(1,T/dt+1); cos(0:dt:T)/5]);
w_i = Interval(5*kron(ones(4,1),[-1/2; 0; -1/5]), kron(ones(4,1),[1/2; 0; 1/5]));

v_tmp = kron(ones(2,1),[sin(0:dt:T).^2/5; zeros(1,T/dt+1); cos(2*(0:dt:T))/3]);
v_h_tmp = kron(ones(2,1),[1/5; 0; 1/3]);
v_l_tmp = kron(ones(2,1),[0; 0; -1/3]);

for i = 1:N
    v{i} = v_tmp(i,:);
    v_i{i} = Interval(v_l_tmp(i), v_h_tmp(i));
end

L_s = eye(N) - diag(ones(N-1,1),-1);
L_s(1,N) = -1;

neigh = cell(N,1);
for i = 1:N
    neigh{i} = find(L_s(i,:));
    neigh{i} = neigh{i}(neigh{i}~=i);
end
gamma = 1;
C{1} = kron([1 0  0 0 0 0],[1,0]);
C{2} = kron([0 1  0 0 0 0],[1,0]);
C{3} = kron([0 0 1 0 0 0],[1,0]);
C{4} = kron([0 0  0 1 0 0],[1,0]);
C{5} = kron([0 0  0 0 1 0],[1,0]);
C{6} = kron([0 0  0 0 0 1],[1,0]);

%%

for i = 1:N
        cvx_begin
            variable L_c(n,m)
            variable T_c(n,n)
            variable G_c(n,m)
            variable E(n,n)
            A_res = T_c*A-L_c*C{i};
            minimize sum(sum(E))
            subject to
                T_c == eye(n) - G_c*C{i}
                -vec(E) <= vec(A_res)
                vec(E) >= vec(A_res)
                vec(E) >= 0
        cvx_end

        L{i} = L_c;
        T_{i} = T_c;
        N_{i} = G_c;
end

for i = 1:N
    for j = 1:n
        ALH{i} = T_{i}*A-L{i}*C{i};
    end
end

%% Centralized gains
C_c = vertcat(C{:});
cvx_begin
    variable L_c(n,N)
    variable T_c(n,n)
    variable N_c(n,N)
    variable E(n,n)
    A_res = T_c*A-L_c*C_c;
    minimize sum(sum(E))
    subject to
        T_c == eye(n) - N_c*C_c
        -vec(E) <= vec(A_res)
        vec(E) >= vec(A_res)
        vec(E) >= 0
cvx_end

ALH_c = T_c*A-L_c*C_c;
v_c = Interval(v_l_tmp, v_h_tmp);

%%
% Use the same x as the other example
for kk = 1:n
    x(kk,:) = downsample(x_t(kk,:),dt/dt_theirs);
end

alpha = [3.5508; 4.6249; 2.6457; 4.8559; 2.7363; 3.0322; 2.3639; 2.6104; 0.7610; 3.2712; 0.6522; 2.2707];
x_0_l = x(:,1) - alpha;
x_0_h = x(:,1) + alpha;


for i = 1:N
    y{i}(:,1) = C{i}*x(:,1) + v{i}(:,1);
    y_c(i,1) = y{i}(:,1);
    x_i{i} = Interval(x_0_l, x_0_h);
    x_c = Interval(x_0_l, x_0_h);
    x_i_best{i} = Interval(x_0_l, x_0_h);
end

err = zeros(5,1);

for d = [5:-1:1, 5]
for t = 1:T/dt-1
    for i = 1:N
        y{i}(:,t+1) = C{i}*x(:,t+1) + v{i}(:,t+1);
        y_c(i,t+1) = y{i}(:,t+1);
    end

    % Local update
    for i = 1:N
        x_i{i}(t+1) = ALH{i}*x_i{i}(t) + L{i}*y{i}(:,t) + N_{i}*y{i}(:,t+1) + (T_{i}*B)*w_i - L{i}*v_i{i} - N_{i}*v_i{i}; 
    end
    
    x_c(t+1) = ALH_c*x_c(t) + L_c*y_c(:,t) + N_c*y_c(:,t+1) + (T_c*B)*w_i - L_c*v_c - N_c*v_c;

    
    % Network update
    for jj = 1:d
        for i = 1:N
            x_i_best{i}(t+1) = x_i{i}(t+1);
            for j = neigh{i}
                x_i_best{i}(t+1) = intersect(x_i_best{i}(t+1), x_i{j}(t+1));
            end
        end
    
        for i = 1:N
                x_i{i}(t+1) = x_i_best{i}(t+1);
        end
    end
end

for i = 1:N
    tmp = horzcat(x_i{i}.h) - horzcat(x_i{i}.l);
    err(d) = err(d) + sum(sum(tmp.^2));
end
err(d) = sqrt(err(d)/(T/dt)/n/6);
end

%%
states = [3,8,10];
% states = 1:n;
ylblx = {'$x_3$','$x_8$','$x_{10}$'};
ylbld = {'$\delta_3 = \overline{x}_3-\underline{x}_3$','$\delta_7 = \overline{x}_7-\underline{x}_7$','$\delta_{10} = \overline{x}_{10}-\underline{x}_{10}$'};

limsx = {[-12,12],[-12,12], [-6,6]};
limsd = {[0, 25], [0, 25],[0,13]};

t = 0:dt:T;
t(end) = [];

figure(1)
tiledlayout(3,2,'TileSpacing','tight','Padding','tight')
% tiledlayout(n,2,'TileSpacing','tight','Padding','tight')
j = 1;
for kk = states
    nexttile;
    plot(t,x(kk,:)','k','linewidth',1)
    hold on
    plot([NaN,NaN],'w')
    for i = 1:N
        x_h_tmp = horzcat(x_i{i}.h);
        x_l_tmp = horzcat(x_i{i}.l);
        plot(t,x_h_tmp(kk,:)','r')

        plot(t_t,x_h_t{i}(kk,:)','b--')
        plot(t,x_l_tmp(kk,:)','r')
        plot(t_t,x_l_t{i}(kk,:)','b--')

    end
    x_h_tmp = horzcat(x_c.h);
    x_l_tmp = horzcat(x_c.l);

%     plot(t,x_h_tmp(kk,:)','k')        
%     plot(t,x_l_tmp(kk,:)','k')

    xlim([0,10])
    ylabel(ylblx{j},'interpreter','latex')
    leg = legend(ylblx{j},'','DIO','[20]','interpreter','latex','location','southeast','orientation','horizontal','FontSize',9,'box','off','numcolumns',2)
    leg.ItemTokenSize = [10,10]
    if kk == 10
        xlabel('Time, $k$','interpreter','latex')
    end
    ylim(limsx{j})
    set(gca,'FontSize',10)
    h = get(gca,'Children');
    set(gca,'Children',h([1:22,23,26,24,25]))
    nexttile;
    xlim([0,10])
    for i = 1:N
        x_h_tmp = horzcat(x_i{i}.h);
        x_l_tmp = horzcat(x_i{i}.l);
        plot(t,x_h_tmp(kk,:)' - x_l_tmp(kk,:)','r')
        hold on
        plot(t_t,x_h_t{i}(kk,:)' - x_l_t{i}(kk,:)','b--')
        legend('DIO','[20]','interpreter','latex','location','northeast','box','off')
    end
%     x_h_tmp = horzcat(x_c.h);
%     x_l_tmp = horzcat(x_c.l);
%     plot(t,x_h_tmp(kk,:)' - x_l_tmp(kk,:)','k')

    xlim([0,10])
    ylim(limsd{j})
    ylabel(ylbld{j},'interpreter','latex')
    set(gca,'FontSize',10)
    j = j + 1;
end
   
xlabel('Time, $k$','interpreter','latex')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4 4.5]);
saveas(gcf,sprintf('img/ex2.eps', kk),'epsc')
saveas(gcf,sprintf('img/ex2.png', kk))
saveas(gcf,sprintf('img/ex2.fig', kk))

%%
figure(1000)
plot(err,'ko-')
hold on

err_c = sqrt(mean(mean((horzcat(x_c.h) - horzcat(x_c.l)).^2)))
plot(1:5,err_c*ones(5,1),'r--')


err_t = 0;
for i = 1:N
    err_t = err_t + sum(sum((x_h_t{i} - x_l_t{i}).^2));
end
err_t = sqrt(err_t/(T/0.01)/n/6)
plot(1:5,err_t*ones(5,1),'b-.')
xlabel('$d$','interpreter','latex')
ylabel('RMS error','interpreter','latex')
legend('DIO','Centralized DIO','[21]', 'interpreter','latex')
leg = legend('DIO','Centralized DIO','[21]','interpreter','latex','location','northeast','orientation','horizontal','FontSize',9,'box','off','numcolumns',3)
leg.ItemTokenSize = [10,10]

ax = gca;
ax.XTick = 1:5;
ylim([0,7.5])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4 1.5]);
saveas(gcf,sprintf('img/e.eps', kk),'epsc')
saveas(gcf,sprintf('img/e.png', kk))
saveas(gcf,sprintf('img/e.fig', kk))

