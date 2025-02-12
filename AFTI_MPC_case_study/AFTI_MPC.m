%%
clear,clc
Ts = 0.05;
NSim = 200;
Np_cases = [5;10;15;20;25]; % prediction horizon
Execution_time_EIQP = zeros(length(Np_cases),1);
MPC_cost = zeros(length(Np_cases),1);
MPC_violations = zeros(length(Np_cases),1);

% Define the continuous-time model
Ac = [-.0151 -60.5651 0 -32.174;
      -.0001 -1.3411 .9929 0;
      .00018 43.2541 -.86939 0;
      0      0       1      0];
Bc = [-2.516 -13.136;
      -.1689 -.2514;
      -17.251 -1.5766;
      0        0];
Cc = [0 1 0 0;
      0 0 0 1];
Dc = [0 0;
      0 0];
sys = ss(Ac,Bc,Cc,Dc);
model = c2d(sys,Ts);
A = model.A;
B = model.B;
C = model.C;
[nx,nu] = size(B);

% MPC-to-QP condense construction
Wy = 10*eye(2);
Wdu = 0.1*eye(2);
u_min = -25*ones(nu,1);
u_max = 25*ones(nu,1);
y_min = [-0.5; -100];
y_max = [0.5; 100];

%% obtain Execution_time_EIQP and MPC_cost
for j=1:length(Np_cases)
    Np = Np_cases(j);
    Ai = A;
    AA = Ai;
    for i=2:Np
        Ai = A*Ai;
        AA = [AA; Ai];
    end
    AiB = B; BB = kron(eye(Np),AiB);
    for i=1:Np-1
        AiB = model.A * AiB;
        BB = BB + kron(diag(ones(Np-i,1),-i),AiB);
    end
    QQ = blkdiag(kron(eye(Np-1),C'*Wy*C),C'*Wy*C);
    
    M = kron(eye(Np),eye(nu)) + kron(diag(ones(Np-1,1),-1),-eye(nu));
    RR = kron(eye(Np),Wdu);
    H = BB'*QQ*BB + M'*RR*M;
    
    Gu = kron(eye(Np),C)*BB;
    Gx = kron(eye(Np),C)*AA;
    %% Closed-loop simulation
    ref = [0; 10];
    x = zeros(nx,1);
    u = zeros(nu,1);
    for k=1:NSim
        if(mod(k,100)==0)
            ref = [0; 0];
        end
        %%
        gi = A*x;
        gg = gi;
        fi = -u;
        ff = fi;
        for i=2:Np
            gi = A*gi;
            gg = [gg; gi];
            fi = zeros(nu,1);
            ff = [ff; fi];
        end
        h = BB'*(QQ*gg - repmat(C'*Wy*ref,Np,1)) + M'*RR*ff;
        lb = repmat(u_min,Np,1);
        ub = repmat(u_max,Np,1);
        Glb = repmat(y_min,Np,1)- Gx*x;
        Gub = repmat(y_max,Np,1)-Gx*x;
        G = [Gu;-Gu];
        g = [Gub;-Glb];
        h_bar = H*lb + h;
        G_bar = [G;eye(Np*nu)];
        g_bar = [g-G*lb; ub-lb];
        tic
        [U,status] = EIQP(H,h_bar,-G_bar,-g_bar,1e-8);
        time_EIQP = toc;
        Execution_time_EIQP(j) = Execution_time_EIQP(j) + time_EIQP;
        U = U + lb;
        MPC_cost(j) = MPC_cost(j) +0.5*(U(1:nu)-u)'*Wdu*(U(1:nu)-u);
        u = U(1:nu);
        MPC_violations(j) = MPC_violations(j) + norm(max(u-u_max,0))+norm(max(u_min-u,0));
        x = model.A*x + model.B*u;
        y = model.C*x;
        MPC_cost(j) = MPC_cost(j) + 0.5*(y-ref)'*Wy*(y-ref);
        MPC_violations(j) = MPC_violations(j) + norm(max(y-y_max,0))+norm(max(y_min-y,0));
    end
    Execution_time_EIQP(j) = Execution_time_EIQP(j)/NSim;
    MPC_cost(j) = MPC_cost(j)/NSim;
    MPC_violations(j) = MPC_violations(j)/NSim;
    
end
Execution_time_EIQP
MPC_cost
MPC_violations

%% plotting
Np = 10;
Ai = A;
AA = Ai;
for i=2:Np
    Ai = A*Ai;
    AA = [AA; Ai];
end
AiB = B; BB = kron(eye(Np),AiB);
for i=1:Np-1
    AiB = model.A * AiB;
    BB = BB + kron(diag(ones(Np-i,1),-i),AiB);
end
QQ = blkdiag(kron(eye(Np-1),C'*Wy*C),C'*Wy*C);
    
M = kron(eye(Np),eye(nu)) + kron(diag(ones(Np-1,1),-1),-eye(nu));
RR = kron(eye(Np),Wdu);
H = BB'*QQ*BB + M'*RR*M;

Gu = kron(eye(Np),C)*BB;
Gx = kron(eye(Np),C)*AA;
% Closed-loop simulation
ref = [0; 10];
x = zeros(nx,1);
u = zeros(nu,1);
Ref_Hist = ref;
U_Hist = [];
Y_Hist = model.C*x;
for k=1:NSim
    if(mod(k,100)==0)
        ref = [0; 0];
    end
    Ref_Hist = [Ref_Hist, ref];
    %%
    gi = A*x;
    gg = gi;
    fi = -u;
    ff = fi;
    for i=2:Np
        gi = A*gi;
        gg = [gg; gi];
        fi = zeros(nu,1);
        ff = [ff; fi];
    end
    h = BB'*(QQ*gg - repmat(C'*Wy*ref,Np,1)) + M'*RR*ff;
    lb = repmat(u_min,Np,1);
    ub = repmat(u_max,Np,1);
    Glb = repmat(y_min,Np,1)- Gx*x;
    Gub = repmat(y_max,Np,1)-Gx*x;
    G = [Gu;-Gu];
    g = [Gub;-Glb];
    h_bar = H*lb + h;
    G_bar = [G;eye(Np*nu)];
    g_bar = [g-G*lb; ub-lb];
    [U,status] = EIQP(H,h_bar,-G_bar,-g_bar,1e-8);
    U = U + lb;
    u = U(1:nu);
    x = model.A*x + model.B*u;
    y = model.C*x;
    U_Hist = [U_Hist, u];
    Y_Hist = [Y_Hist, y];
end
TSim = 0:Ts:NSim*Ts;
figure()
subplot(411)
plot(TSim,Y_Hist(2,:),'r','linewidth',1.5),ylabel('y2')
hold on 
stairs(TSim,Ref_Hist(2,:),'b','linewidth',1.5)
subplot(412)
plot(TSim,Y_Hist(1,:),'r','linewidth',1.5),ylabel('y1'),
subplot(413)
plot(TSim,[U_Hist(1,:),NaN],'r','linewidth',1.5),ylim([-25,25]),ylabel('u1'),yticks([-25,0,25])
subplot(414)
plot(TSim,[U_Hist(2,:),NaN],'r','linewidth',1.5),ylim([-25,25]), xlabel('time'),ylabel('u2'),yticks([-25,0,25])
