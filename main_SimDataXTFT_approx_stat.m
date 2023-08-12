%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Charge transfer by discrete breathers in 1D model
% written by Janis Bajars, February 2023
% Code to obtain numerical data of approximate stationary polarobreather
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute parameter values in structure <parm> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameter_Values_SimDataXTFT_approx_stat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define simulation variables (u,a,p,b) in structure <vars>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Define_Variables;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions: stationary breather + charge
L = round(parm.N/2);
% Lattice
gamma = 0.4;
vars.p(L-1) = -1*gamma;
vars.p(L) = 2*gamma;
vars.p(L+1) = -2*gamma;
vars.p(L+2) = 1*gamma;
% Charge
vars.a(L) = 1;
vars.b(L) = -1;

% Make sure that the charge constraint is satisfied
C2 = sum(vars.a.^2 + vars.b.^2);
vars.a = vars.a/sqrt(C2)*parm.tau_f;
vars.b = vars.b/sqrt(C2)*parm.tau_f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time integration: initial evolution of Nsteps0 steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for n=1:parm.Nsteps0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Semi-implicit splitting method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    vars = Method_Im_PQDWimDQP(parm,vars,parm.h);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices for saving data in time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = zeros(parm.N,parm.Nsteps+1); 
U(:,1) = vars.u;
Q = zeros(parm.N,parm.Nsteps+1); 
Q(:,1) = parm.x + vars.u;
P = zeros(parm.N,parm.Nsteps+1); 
P(:,1) = vars.p;
A = zeros(parm.N,parm.Nsteps+1); 
A(:,1) = vars.a;
B = zeros(parm.N,parm.Nsteps+1); 
B(:,1) = vars.b;
H = zeros(parm.Nsteps+1,1); 
H(1,1) = Comp_Hamiltonian(parm,vars);
C2 = zeros(parm.Nsteps+1,1); 
C2(1,1) = Comp_Charge_Total(parm,vars);
C2_probability = zeros(parm.N,parm.Nsteps+1); 
C2_probability(:,1) = Comp_Charge_Probability(parm,vars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time integration: evolution of Nsteps steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:parm.Nsteps
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Semi-implicit splitting method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    vars = Method_Im_PQDWimDQP(parm,vars,parm.h);
    
    % Save data in time
    U(:,n+1) = vars.u;
    Q(:,n+1) = parm.x + vars.u;
    P(:,n+1) = vars.p;
    A(:,n+1) = vars.a;
    B(:,n+1) = vars.b;
    H(n+1,1) = Comp_Hamiltonian(parm,vars);
    C2(n+1,1) = Comp_Charge_Total(parm,vars);
    C2_probability(:,n+1) = Comp_Charge_Probability(parm,vars);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save simulation data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('SavedData/Approx_SimData_Stat_G04_E0shift_N64Tau001.mat',...
    'parm','U','P','A','B')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results and errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
set(gcf, 'Position',  [100, 550, 700, 400])
colormap(flipud(hot))
contourf(parm.t,parm.x,C2_probability,'edgecolor','none')
shading flat
caxis([0 1])
c=colorbar;
set(c,'TickLabelInterpreter','latex')
hold on
plot(parm.t,Q','-k','linewidth',1)
axis on
box on
xlabel('$t$')
ylabel('$q_n^0+1.5u_n$','interpreter','latex')
axis([0 parm.Tend 0 parm.N-1])
legend('Charge probability $P_n$')
set(legend,'location','northwest')

figure(2)
set(gcf, 'Position',  [900, 550, 700, 400])
ErrorH = abs((H-H(1))/H(1));
plot(parm.t,ErrorH,'-r','linewidth',1.5)
axis on
box on
grid on
xlabel('$t$')
ylabel('Relative Hamiltonian error')

figure(3)
set(gcf, 'Position',  [500, 50, 700, 400])
ErrorC = abs((C2-C2(1))/C2(1));
plot(parm.t,ErrorC,'-b','linewidth',1.5)
axis on
box on
grid on
xlabel('$t$')
ylabel('Relative probability error')

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove path
rmpath('Functions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
