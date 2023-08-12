%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Charge transfer by discrete breathers in 1D model
% written by Janis Bajars, February 2023
% Damped Gauss-Newton algorithm to obtain exact stationary polarobreather
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clearvars
clc
set(0,'DefaultAxesFontSize',16)
set(groot,'DefaultAxesTickLabelInterpreter','latex'); 
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

disp('Start of the computation!')
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add path
addpath('Functions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load approximate stationary polarobreather simulation data
load('SavedData/Approx_SimData_Stat_G04_E0shift_N64Tau001.mat',...
    'parm','U','P','A','B')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redefine some parameter values in structure <parm> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of time steps
parm.Nsteps = 1200;
% Final simulation time
parm.Tend = parm.h*parm.Nsteps;
% Time grid points
parm.t = 0:parm.h:parm.Tend;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define simulation variables in structure <vars>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Define_Variables;

% Pick an initial seed (solution) such that p~0
Nn = 165; % gamma=0.4
vars.u = U(:,Nn);
vars.a = A(:,Nn);
% Either one of the options for momenta p
vars.p = zeros(parm.N,1);
% vars.p = P(:,Nn);
vars.b = B(:,Nn);
vars.b(end) = 0;

% % Set displacements to zero: optional
% vars.u(1:parm.N/4) = 0;
% vars.u(3*parm.N/4:end) = 0;

% Make sure that the charge constraint is satisfied
C2 = sum(vars.a.^2 + vars.b.^2);
vars.a = vars.a/sqrt(C2)*parm.tau_f;
vars.b = vars.b/sqrt(C2)*parm.tau_f;

% Save initial values
u_n = vars.u;
a_n = vars.a;
p_n = vars.p;
b_n = vars.b;

% Initial E0 value
E0 = parm.E0;
format long
disp(E0)
format short

% Max number of iterations
Iter_max = 100;
% Save iteration and error values in structure <parm>
parm.iter_max = Iter_max;
parm.iter = zeros(Iter_max,1);
parm.error_p = zeros(Iter_max,1);
parm.error_b = zeros(Iter_max,1);
parm.error_C2 = zeros(Iter_max,1);
parm.error_f = zeros(Iter_max,1);
parm.values_E0 = zeros(Iter_max,1);
parm.mu = zeros(Iter_max,1);

% Errors 
Error_f = Inf;
Error_C2 = Inf;
% Tolerance
Tol = 1e-14;
% Numerical Jacobian computation
delta = 1e-6;
df = zeros(4*parm.N,4*parm.N+1);
% Number of constraints
parm.M = parm.N + 2;
% Jacobian matrix of linear constraints for p,
% total charge probability and b_N=0
dg = [zeros(parm.N,parm.N) zeros(parm.N,parm.N)... 
    eye(parm.N,parm.N) zeros(parm.N,parm.N+1)];
dg = [dg; zeros(2,4*parm.N+1)];
dg(end,end-1) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time integration to find f and the objective function F_old 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:parm.Nsteps
    vars = Method_Im_PQDWimDQP(parm,vars,parm.h); 
end
U_n = vars.u;
A_n = vars.a;
P_n = vars.p;
B_n = vars.b;
% Compute function value f 
f = [U_n-u_n; A_n-a_n; P_n-p_n; B_n-b_n];
% Objective function
F_old = (f'*f)/2;

iter = 1;
while Error_f > Tol || Error_C2 > Tol

    disp(iter)
    
% % % % %     % Optional: plot solutions at different iterations
% % % % %     figure(iter)
% % % % %     set(gcf, 'Position',  [600, 150, 1200, 800])
% % % % %         
% % % % %     subplot(221)
% % % % %     hold on
% % % % %     plot(parm.x,u_n,'-k','linewidth',1.5)
% % % % %     plot(parm.x,U_n,'--r','linewidth',1.5)
% % % % %     axis on
% % % % %     box on
% % % % %     grid on
% % % % %     xlabel('$n$')
% % % % %     ylabel('$u_n$')
% % % % %     axis([0 parm.N-1 -0.15 0.1501])
% % % % %         
% % % % %     subplot(222)
% % % % %     hold on
% % % % %     plot(parm.x,a_n,'-k','linewidth',1.5)
% % % % %     plot(parm.x,A_n,'--r','linewidth',1.5)
% % % % %     plot(parm.x,b_n,'-b','linewidth',1.5)
% % % % %     plot(parm.x,B_n,'--c','linewidth',1.5)
% % % % %     axis on
% % % % %     box on
% % % % %     grid on
% % % % %     xlabel('$n$')
% % % % %     ylabel('$a_n$, $b_n$')
% % % % %     axis([0 parm.N-1 -0.04 0.04])
% % % % %         
% % % % %     subplot(223)
% % % % %     hold on
% % % % %     plot(parm.x,p_n,'-k','linewidth',1.5)
% % % % %     plot(parm.x,P_n,'--r','linewidth',1.5)
% % % % %     axis on
% % % % %     box on
% % % % %     grid on
% % % % %     xlabel('$n$')
% % % % %     ylabel('$p_n$')
% % % % %     axis([0 parm.N-1 -0.15 0.1501])
% % % % % 
% % % % %     subplot(224)
% % % % %     hold on
% % % % %     plot(parm.x,a_n.^2+b_n.^2,'-k','linewidth',1.5)
% % % % %     plot(parm.x,A_n.^2+B_n.^2,'--r','linewidth',1.5)
% % % % %     axis on
% % % % %     box on
% % % % %     grid on
% % % % %     xlabel('$n$')
% % % % %     ylabel('$|c_n|^2$')
% % % % %     axis([0 parm.N-1 -0.0001 0.001])
% % % % %     drawnow
       
    % Numerical Jacobian computation 
    X = [u_n; a_n; p_n; b_n; E0];
       
    for j=1:4*parm.N+1
        
        % Modify one of the entries of X
        XX = X;        
        temp = X(j);
        h = delta*abs(temp);
        if h < 1e-8
            h = delta;
        end
        XX(j) = temp+h;
        h = XX(j)-temp;
        
        % Initial conditions for time integration
        vars.u = XX(1:parm.N);
        vars.a = XX(parm.N+1:2*parm.N);    
        vars.p = XX(2*parm.N+1:3*parm.N);
        vars.b = XX(3*parm.N+1:4*parm.N);
        parm.E0 = XX(end);
        parm.E0_tau = parm.E0/parm.tau;
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Time integration with new u, p, a, b, and E0 values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for n=1:parm.Nsteps
            vars = Method_Im_PQDWimDQP(parm,vars,parm.h); 
        end
                 
        % Compute function value g
        g = [vars.u-XX(1:parm.N);...
            vars.a-XX(parm.N+1:2*parm.N);...
            vars.p-XX(2*parm.N+1:3*parm.N);...
            vars.b-XX(3*parm.N+1:4*parm.N)];
        
        % Numerical Jacobian entry
        df(:,j) = (g-f)/h;
        
    end    
    
    % Compute constraint Jacobian matrix dg
    dg(parm.M-1,parm.N+1:2*parm.N) = a_n';
    dg(parm.M-1,3*parm.N+1:4*parm.N) = b_n';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Damped Gauss-Newton algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    % Regularization parameter mu
    if iter==1
        mu = 1e-6*max(diag(df'*df));
    end

    % Vector G = df'*f
    G = df'*f;

    % System matrix M
    M = [df'*df+mu*eye(4*parm.N+1,4*parm.N+1) dg'; dg zeros(parm.M,parm.M)];
    % Right-hand side
    Y = -[G; p_n; (a_n'*a_n + b_n'*b_n)/2-parm.tau; b_n(end)];
    % Solve the linear system of equations
    dX = M\Y;
    dX = dX(1:4*parm.N+1);
    % Update variables
    X = X + dX;

    % Initial conditions for the next iteration and time integration
    vars.u = X(1:parm.N);
    vars.a = X(parm.N+1:2*parm.N);    
    vars.p = X(2*parm.N+1:3*parm.N);
    vars.b = X(3*parm.N+1:4*parm.N);
    parm.E0 = X(end); 
    parm.E0_tau = parm.E0/parm.tau;

    % Save initial values
    u_n = vars.u;
    a_n = vars.a;
    p_n = vars.p;
    b_n = vars.b;

    % Value of E0
    format long
    E0 = parm.E0;
    disp(E0)
    format short
    parm.values_E0(iter,1) = E0; 
    
    disp('_______________________________________________________________')

    % Errors
    Error_p = max(abs(vars.p));
    parm.error_p(iter,1) = Error_p;
    disp(Error_p)

    Error_b = abs(vars.b(end));
    parm.error_b(iter,1) = Error_b;
    disp(Error_b)

    Error_C2 = abs(sum(vars.a.^2 + vars.b.^2)/2/parm.tau-1);
    parm.error_C2(iter,1) = Error_C2;
    disp(Error_C2)
               
    Error_f = max(abs(f));
    parm.error_f(iter,1) = Error_f;
    disp(Error_f)

    parm.mu(iter,1) = mu;
    disp(mu)

    disp('---------------------------------------------------------------')
            
    parm.iter(iter,1) = iter; 
        
    iter = iter + 1;
    if iter > Iter_max
        disp('Reached the maximum number of iterations!')        
        break
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time integration to find f and the objective function F_new 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n=1:parm.Nsteps
        vars = Method_Im_PQDWimDQP(parm,vars,parm.h); 
    end
    U_n = vars.u;
    A_n = vars.a;
    P_n = vars.p;
    B_n = vars.b;
    % Compute function value f 
    f = [U_n-u_n; A_n-a_n; P_n-p_n; B_n-b_n];
    % Objective function
    F_new = (f'*f)/2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the gain ratio rho and update the mu value
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho = 2*(F_old-F_new)/(dX'*(mu*dX-G));
    F_old = F_new;
    if rho < 0.25
        mu = 2*mu;
    elseif rho > 0.75
        mu = mu/3;
    end
        
end

% Save simulation data 
save('SavedData/ExactSol_Stat_G04_E0shift_N64Tau001.mat','parm','vars')

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove path
rmpath('Functions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
