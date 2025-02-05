%%
clear, close, clc
rng(20250131)
cr = [1e1,1e2,5e2,1e3,5e3,1e4,5e4,1e5,5e5,1e6];
num_cr = length(cr);

InfeasibiltiyDetectingRate_quadprog = zeros(num_cr,1);
InfeasibiltiyDetectingRate_OSQP_default = zeros(num_cr,1);
InfeasibiltiyDetectingRate_OSQP_maxiter5e5 = zeros(num_cr,1);
InfeasibiltiyDetectingRate_EIQP = zeros(num_cr,1);

%%
num_exp = 100;
for i=1:num_cr
    cr_i = cr(i);
    m = 20;
    n = 4*m;
    for j=1:num_exp
        v = [cr_i; (cr_i-1)*rand(m-2,1)+1; 1];
        [U,~] = qr(rand(m,m));
        Q = U*diag(v)*U';
        c = ones(m,1);
        A = full(2*sprand(n,m,0.15));
        b = ones(n,1);
        % To make the QP infeasible, we add the contradictory constraint
        % A = [A; -A(1,:)];
        % b = [b; -b(1)-1];
        A = [A; -A(1,:);-A(2,:)];
        b = [b; -b(1)-1;-b(2)-1];
        %% quadprog
        options = optimoptions('quadprog','Display','off');
        [z, ~, existflag] = quadprog((Q+Q')/2,c,A,b,[],[],[],[],[],options);
        if existflag==-2
            InfeasibiltiyDetectingRate_quadprog(i) = InfeasibiltiyDetectingRate_quadprog(i) + 1;
        end
        %% OSQP_default, please download OSQP from https://osqp.org/
        prob = osqp;
        prob.setup(Q, c, A, [], b,'verbose', false);
        res = prob.solve();
        if strcmp(res.info.status,'primal infeasible')
            InfeasibiltiyDetectingRate_OSQP_default(i) = InfeasibiltiyDetectingRate_OSQP_default(i) + 1;
        end
        %% OSQP_maxiter5e5, please download OSQP from https://osqp.org/
        prob_maxiter5e5 = osqp;
        prob_maxiter5e5.setup(Q, c, A, [], b,'verbose', false,'max_iter',5e5,'eps_abs',1e-10,'eps_rel',1e-10,'polish',true);
        res_maxiter5e5 = prob_maxiter5e5.solve();
        if strcmp(res_maxiter5e5.info.status,'primal infeasible')
            InfeasibiltiyDetectingRate_OSQP_maxiter5e5(i) = InfeasibiltiyDetectingRate_OSQP_maxiter5e5(i) + 1;
        end

        %% EIQP: min 1/2 z'*Q*z+c'*z, s.t. Az>=b, z>=0
        Q_bar = [Q,-Q;-Q, Q];
        c_bar = [c; -c];
        A_bar = [-A, A];
        b_bar = -b;
        [z_EIQP,status_EIQP] = EIQP(Q_bar,c_bar,A_bar,b_bar,1e-6);
        if status_EIQP=="Optimal"
            disp("Infeasibility detecting failed")
        else
            InfeasibiltiyDetectingRate_EIQP(i) = InfeasibiltiyDetectingRate_EIQP(i) + 1;
        end
    end
    InfeasibiltiyDetectingRate_quadprog(i) = InfeasibiltiyDetectingRate_quadprog(i)/num_exp;
    InfeasibiltiyDetectingRate_OSQP_default(i) = InfeasibiltiyDetectingRate_OSQP_default(i)/num_exp;
    InfeasibiltiyDetectingRate_OSQP_maxiter5e5(i) = InfeasibiltiyDetectingRate_OSQP_maxiter5e5(i)/num_exp;
    InfeasibiltiyDetectingRate_EIQP(i) = InfeasibiltiyDetectingRate_EIQP(i)/num_exp;
end

%% plot
figure(1)
semilogx(cr, InfeasibiltiyDetectingRate_quadprog, '--o', 'LineWidth', 2.0);
hold on
semilogx(cr,InfeasibiltiyDetectingRate_OSQP_default,'--o','LineWidth', 2.0)
semilogx(cr,InfeasibiltiyDetectingRate_OSQP_maxiter5e5,'--o','LineWidth', 2.0)
semilogx(cr,InfeasibiltiyDetectingRate_EIQP,'--o','LineWidth',2.0)
legend('quadprog','OSQP-default','OSQP-5e5','EIQP')
xlabel('Condition number of Q'),ylabel('Infeasibility Detecting Rate')
grid on
