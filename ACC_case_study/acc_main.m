clc;
x0 = [0 20 100]; % ego vehicle state: x(position), v(speed), z(distance with respect to the preceding vehicle)
g = 9.81; m = 1650;
f0 = 0.1; f1 = 5; f2 = 0.25; 
vd = 24; % desired speed
eps = 10;
gamma = 1;
ca = 0.4; cd = 0.4;%cd = 0.375;
psc = 1;
c = 10;
set(0,'DefaultTextInterpreter','latex')
warning('off');
global u v0;
v0 = 13.89;  % the speed of the preceding vehicle (constant)


mode = 0;  % control model: 0 - without minimum braking distance, 1-with minimum braking distance.
use_conversion = 0;  % use the conversion method or not
time1 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%quadprog
for i = 1:300
    i
Fr = f0 + f1*x0(2) + f2*x0(2)^2;
phi0 = -2*(x0(2) - vd)*Fr/m + eps*(x0(2) - vd)^2;
phi1 = 2*(x0(2) - vd)/m;

k = 2;
if mode == 0  % HOCBF b = z - c without minimum braking distance
    % LfB =v0 - x0(2) +1.8*Fr/m + k*(x0(3) - 1.8*x0(2));   % 1st degree CBF
    % LgB = 1.8/m;
    k = 1;   % HOCBF
    b = x0(3) - c;
    Lfb = v0 - x0(2);
    Lf2b = Fr/m;  % - 
    LgLfb = 1/m;
    LfB = Lf2b + 2*k*Lfb + k^2*b;
    LgB = LgLfb;
else  % CBF b = z - 1.8v with minimum braking distance
    LfB = v0 - x0(2) - ((v0 - x0(2))/(cd*g) - 1.8)*Fr/m + x0(3) - 1.8*x0(2) - 0.5*(v0 - x0(2))^2/(cd*g);
    LgB = -((v0 - x0(2))/(cd*g) - 1.8)/m;
end


A = [phi1 -1;LgB 0;1 0; -1 0];
b = [-phi0;LfB;ca*m*g; cd*m*g];
H = [2/(m^2) 0; 0 2*psc];
F = [-2*Fr/(m^2); 0];
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off');  % interior-point-convex or active-set
tic
[u,fval,exitflag,output,lambda] = ...
   quadprog(H,F,A,b,[],[],[],[],[0;0],options);
temp = toc;
time1 = [time1,temp];
t=[0 0.1];
[tt,xx]=ode45('acc',t,x0);
x0 = [xx(end, 1) xx(end, 2) xx(end, 3)];
result1(i,:) = [0.1*i xx(end,2) u(1) x0(3)-c]; %x0(3)-1.8*x0(2)
end


x0 = [0 20 100];
time2 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EIQP
for i = 1:300
    i
Fr = f0 + f1*x0(2) + f2*x0(2)^2;   % CLF for desired speed
phi0 = -2*(x0(2) - vd)*Fr/m + eps*(x0(2) - vd)^2;
phi1 = 2*(x0(2) - vd)/m;

k = 2;
if mode == 0  % HOCBF b = z - c, without minimum braking distance
    % LfB =v0 - x0(2) +1.8*Fr/m + k*(x0(3) - 1.8*x0(2));   % 1st degree CBF
    % LgB = 1.8/m;
    k = 1;
    b = x0(3) - c;
    Lfb = v0 - x0(2);
    Lf2b = Fr/m;  % - 
    LgLfb = 1/m;
    LfB = Lf2b + 2*k*Lfb + k^2*b;
    LgB = LgLfb;
else  % CBF b = z - 1.8v with minimum braking distance
    LfB = v0 - x0(2) - ((v0 - x0(2))/(cd*g) - 1.8)*Fr/m + x0(3) - 1.8*x0(2) - 0.5*(v0 - x0(2))^2/(cd*g);
    LgB = -((v0 - x0(2))/(cd*g) - 1.8)/m;
end

%%%%%% note the standard QP form in Matlab
%min u^THu + F^Tu, s.t., Au <= b
%%%%%%

if use_conversion == 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%with conversion
%%%%%%%%%%%% standard QP form in Matlab: min u^THu + F^Tu, s.t., Au <= b
A = [phi1 -1;LgB 0;1 0; -1 0];
b = [-phi0;LfB;ca*m*g; cd*m*g];
H = [2/(m^2) 0; 0 2*psc];
F = [-2*Fr/(m^2); 0];

%%%%%%%%%%%%%%%%%%%Convert to the QP form: min u^THu + F^Tu, s.t., Au >= b, u>=0
H_bar= [H, -H; -H, H];
F_bar = [F; -F];
A_bar = [-A, A];
b_bar = -b;

%solve QP with EIQP
epsilon = 1e-16; 
tic
% [u, status] = EIQP(H,F,A,b,epsilon);
[u, status] = EIQP(H_bar,F_bar,A_bar,b_bar,epsilon);
u = u(1:2) - u(3:4);
temp = toc;

else%%%%%%%%%%%%%%%%%%%%without conversion (set u <- u - lb), directly write in the form of: 
%%%%%%%%%%%% min u^THu + F^Tu, s.t., Au >= b, u>=0
%%%%%%%%%%%%where the upper bound of control is put in Au >= b, the lower
%%%%%%%%%%%%control bound is specified as u >= 0 (as we have set u <- u - lb)
lb = [-cd*m*g; 0];
A = [phi1 -1;LgB 0;1 0];
b = [-phi0;LfB;ca*m*g] - A*lb;  %% as u <- u - lb
H = [2/(m^2) 0; 0 2*psc];
F = [-2*Fr/(m^2); 0];
F = F + H*lb;                   %% as u <- u - lb
epsilon = 1e-16;  % 1e-18
tic
[u, status] = EIQP(H,F,-A,-b,epsilon);
u = u + lb;     %% get the true control, as u <- u - lb
temp = toc;
end

time2 = [time2, temp];
t=[0 0.1];
[tt,xx]=ode45('acc',t,x0);
x0 = [xx(end, 1) xx(end, 2) xx(end, 3)];
result2(i,:) = [0.1*i xx(end,2) u(1) x0(3)-c];
end

fprintf('QuadProg time mean, std: %f, %f\n', mean(time1), std(time1))
fprintf('EIQP time mean, std: %f, %f\n', mean(time2), std(time2))


figure(1)
subplot(3, 1, 1)

hold on
plot(result1(:, 1),result1(:, 2),'b-', result2(:, 1),result2(:, 2),'r--',  'LineWidth', 2)
plot([0,30],[vd,vd], 'k--',[0,30],[13.89,13.89], 'k--')
xlabel('$t/s$','fontsize',15)
ylabel('$v(t)$', 'fontsize',15)
lg = legend('Quadprog', 'EIQP');
set(lg, 'fontsize', 10)
axis([0 30 10 26]); 
grid on
subplot(3, 1, 2)
plot([0,30],[ca,ca], 'k--',[0,30],[-cd,-cd], 'k--')
hold on
plot(result1(:, 1),result1(:, 3)/(m*g), 'b-', result2(:, 1),result2(:, 3)/(m*g), 'r--',  'LineWidth', 2)
axis([0 30 -0.75 0.75]); 
xlabel('$t/s$','fontsize',15)
ylabel('$u(t)/(mg)$', 'fontsize',15)
grid on
subplot(3, 1, 3)
plot(result1(:, 1),result1(:, 4), 'b-', result2(:, 1),result2(:, 4), 'r--', [0,30],[0,0], 'k--', 'LineWidth', 2)
axis([0 30 -5 75]); 
grid on
xlabel('$t/s$','fontsize',15)
ylabel('$b(x(t))$', 'fontsize',15)


 