% This code implements the observer example from
% https://ieeexplore.ieee.org/document/9713996
% 
% Author: Scott Brown (sab007@ucsd.edu)
% 2022
%
% Run this code first, before running the other, to plot the comparison

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

T = 30;
dt = 0.01;

N = 6;
n = 12;
m = 1;

Dn = [eye(n) -eye(n)];

A = kron(eye(N),[0,1;0,0]);

L = eye(N) - diag(ones(N-1,1),-1);
L(1,N) = -1;
W = abs(L);

neigh = cell(N,1);
for i = 1:N
    neigh{i} = find(L(i,:));
    neigh{i} = neigh{i}(neigh{i}~=i);
end
gamma = 1;
C{1} = kron([1 0  0 0 0 0],[1,0]);
C{2} = kron([0 1  0 0 0 0],[1,0]);
C{3} = kron([0 0 1 0 0 0],[1,0]);
C{4} = kron([0 0  0 1 0 0],[1,0]);
C{5} = kron([0 0  0 0 1 0],[1,0]);
C{6} = kron([0 0  0 0 0 1],[1,0]);

Kd{1} = [-2; -1];
Kd{2} = [-2; -1];
Kd{3} = [-6; -5];
Kd{4} = [-4.2426; -3.5355];
Kd{5} = [-6; -5];
Kd{6} = [-6; -5];

for i = 1:N
%     [~, ~, ~, T_{i}] = obsvf(A,eye(n),C{i});
%     U{i} = T_{i}(:,3:end);
%     D{i} = T_{i}(:,1:2);
    U{i} = eye(n);
    D{i} = U{i}(:,2*(i-1)+1:2*i);
    U{i}(:,2*(i-1)+1:2*i) = [];
    T_{i} = [D{i} U{i}];
    F{i} = [zeros(2,n); U{i}'];
    
    K{i} = T_{i}*[Kd{i}; zeros(n-2,1)];
    K_p{i} = K{i}.*(K{i}>0);
    K_m{i} = K_p{i} - K{i};
    
    At{i} = T_{i}'*A*T_{i};
    Ad{i} = At{i}(1:2,1:2);
    Ar{i} = At{i}(3:end,1:2);
    Au{i} = At{i}(3:end,3:end);
    
    Ct{i} = C{i}*T_{i};
    Cd{i} = Ct{i}(:,1:2);
    
    [S{i}, J{i}] = jordan(Ad{i}+Kd{i}*Cd{i});
    
    St{i} = blkdiag(S{i},eye(10));
    
    Ac{i} = [J{i} zeros(2,10);Ar{i}*S{i} Au{i}];
    B{i} = [D{i}*S{i} U{i}];
end

c = 90;


%%
x0 = kron(ones(N,1),[0.7032; 0.0457]);
x(:,1) = x0;

alpha = [3.5508; 4.6249; 2.6457; 4.8559; 2.7363; 0.0322; 2.3639; 2.6104; 0.7610; 3.2712; 0.6522; 2.2707];
x_0_l = x(:,1) - alpha;
x_0_h = x(:,1) + alpha;

w = kron(ones(4,1),[sin(0:dt:T)/2; zeros(1,T/dt+1); cos(0:dt:T)/5]);
w_l = kron(ones(4,1),[-1/2; 0; -1/5]);
w_h = kron(ones(4,1),[1/2; 0; 1/5]);

v_tmp = kron(ones(2,1),[sin(0:dt:T).^2/5; zeros(1,T/dt+1); cos(2*(0:dt:T))/3]);
v_h_tmp = kron(ones(2,1),[1/5; 0; 1/3]);
v_l_tmp = kron(ones(2,1),[0; 0; -1/3]);
for i = 1:N
    v{i} = v_tmp(i,:);
    v_l{i} = v_l_tmp(i);
    v_h{i} = v_h_tmp(i);
end

for i = 1:N
    y{i}(:,1) = C{i}*x(:,1) + v{i}(:,1);
    Z_l{i} = mp(inv(St{i})*inv(T_{i}))*G(x_0_l,x_0_l);
    Z_h{i} = mp(inv(St{i})*inv(T_{i}))*G(x_0_l,x_0_h);
end

%%
% 
% clear x
% Z_h_d = @(i,Z_l,Z_h,y,t) ...
%     Mz(Ac{i} - c*L(i,i)*F{i}*B{i})*Z_h(idx(i)) ...
%       + c*mp(F{i})*W(i,neigh{i})*mp(B{neigh{i}}*Dn)*G(Z_l(idx(neigh{i})),Z_h(idx(neigh{i}))) ...
%       + mp(-inv(St{i})*inv(T_{i})*K{i})*p(y(i)) ...
%       + mp(inv(St{i})*inv(T_{i}))*G(w_l,w_h) ...
%       + mp(inv(St{i})*inv(T_{i}))*G(K_p{i}*v_l{i}-K_m{i}*v_h{i},K_p{i}*v_h{i}-K_m{i}*v_l{i});
%   
% Z_l_d = @(i,Z_l,Z_h,y,t) ...
%     Mz(Ac{i} - c*L(i,i)*F{i}*B{i})*Z_l(idx(i)) ...
%       + c*mp(F{i})*W(i,neigh{i})*mp(B{neigh{i}}*Dn)*G(Z_l(idx(neigh{i})),Z_l(idx(neigh{i}))) ...
%       + mp(-inv(St{i})*inv(T_{i})*K{i})*p(y(i)) ...
%       + mp(inv(St{i})*inv(T_{i}))*G(w_l,w_l) ...
%       + mp(inv(St{i})*inv(T_{i}))*G(K_p{i}*v_l{i}-K_m{i}*v_h{i},K_p{i}*v_l{i}-K_m{i}*v_h{i});
% 
% x_d = @(x,t) A*x + w(t);
% 
% zd = @(t,x) [x_d(x(1:n),t);
%     Z_l_d(1,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_l_d(2,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_l_d(3,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_l_d(4,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_l_d(5,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_l_d(6,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_h_d(1,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_h_d(2,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_h_d(3,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_h_d(4,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_h_d(5,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t);
%     Z_h_d(6,x(n+1:13*n),x(13*n+1:end),vertcat(C{:})*x(1:n)+v(t),t)];
% 
% tmp = ode45(zd,[0,30],[x0;vertcat(Z_l_0{:});vertcat(Z_h_0{:})])
% 
% %%
% for i = 1:N
%     t_t = tmp.x
%     x = tmp.y(1:n,:)
%     Z_l{i} = tmp.y(n+idx(i),:)
%     Z_h{i} = tmp.y(13*n+idx(i),:)
% end

%%
for t = 1:T/dt-1
    x(:,t+1) = x(:,t) + dt*(A*x(:,t) + w(:,t));
    for i = 1:N
        y{i}(:,t) = C{i}*x(:,t) + v{i}(:,t);
    end
    
    for i = 1:N
        tmp = 0;
        for j = neigh{i}
            tmp = tmp + W(i,j)*mp(B{j}*Dn)*G(Z_l{j}(:,t),Z_l{j}(:,t));
        end
        Z_l{i}(:,t+1) = Z_l{i}(:,t) + dt*(Mz(Ac{i}-c*L(i,i)*F{i}*B{i})*Z_l{i}(:,t) + c*mp(F{i})*tmp + mp(-inv(St{i})*inv(T_{i})*K{i})*p(y{i}(:,t)) ...
            + mp(inv(St{i})*inv(T_{i}))*G(w_l,w_l) + mp(inv(St{i})*inv(T_{i}))*G(K_p{i}*v_l{i}-K_m{i}*v_h{i},K_p{i}*v_l{i}-K_m{i}*v_h{i}));
        
        tmp = 0;
        for j = neigh{i}
            tmp = tmp + W(i,j)*mp(B{j}*Dn)*G(Z_l{j}(:,t),Z_h{j}(:,t));
        end
        Z_h{i}(:,t+1) = Z_h{i}(:,t) + dt*(Mz(Ac{i}-c*L(i,i)*F{i}*B{i})*Z_h{i}(:,t) + c*mp(F{i})*tmp + mp(-inv(St{i})*inv(T_{i})*K{i})*p(y{i}(:,t)) ...
            + mp(inv(St{i})*inv(T_{i}))*G(w_l,w_h) + mp(inv(St{i})*inv(T_{i}))*G(K_p{i}*v_l{i}-K_m{i}*v_h{i},K_p{i}*v_h{i}-K_m{i}*v_l{i}));
    end
    
end

%%
for i = 1:N
    x_h_t{i} = [eye(n) zeros(n)]*mp(T_{i}*St{i})*Z_h{i} - [zeros(n) eye(n)]*mp(T_{i}*St{i})*Z_l{i};
    x_l_t{i} = [eye(n) zeros(n)]*mp(T_{i}*St{i})*Z_l{i} - [zeros(n) eye(n)]*mp(T_{i}*St{i})*Z_h{i};
%     x_h{i} = Dn*mp(T_{i}*St{i})*Z_h{i};
%     x_l{i} = Dn*mp(T_{i}*St{i})*Z_l{i};
end
t_t = 0:dt:T;
t_t(end) = [];
x_t = x;
save('theirs.mat','x_h_t','x_l_t','t_t','x_t')
%%
ylbl = {'$x_1$','$x_2$','$x_3$','$x_4$'};
lims = {[-100, 100], [-44 67],[-15 17], [-14,14]};

figure(1)
    tiledlayout(3,4,'TileSpacing','tight','Padding','tight')
    for kk = 1:n
        nexttile;
        plot(x(kk,:)','k','linewidth',1)
        hold on
        for i = 1:N
            plot(real(x_h_t{i}(kk,:))','r--')
            plot(real(x_l_t{i}(kk,:))','b--')
        end
        ylim([-10,10])
    end
%%
Mz(Ac{i} - c*L(i,i)*F{i}*B{i})
%%
function res = mp(A)
    Ap = A.*(A>0);
    Am = Ap-A;
    res = [Ap Am;Am Ap];
end

function res = Mz(A)
    dA = diag(diag(A));
    Ap = (A-dA).*((A-dA)>0);
    Am = Ap - (A-dA);
    res = [dA + Ap, Am;Am dA+Ap];
end

function res = G(a,x)
    am = -a.*(a<0);
    res = [x+am;am];
end

function res = p(x)
    xp = x.*(x>0);
    xm = xp-x;
    res = [xp;xm];
end

function res = idx(i)
    n = 12;
    res = 2*n*(i-1)+1:2*n*i;
end

% function res = v(t)
%     res = kron(ones(2,1),[sin(t).^2/5; 0; cos(2*(t))/3])/10;
% end
% 
% function res = w(t)
%     res = kron(ones(4,1),[sin(t)/2; 0; cos(t)/5])/10;
% end