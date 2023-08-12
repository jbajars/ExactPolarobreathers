%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute parameter values and save in structure <parm>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time step
parm.h = 0.01;
% Number of time steps
parm.Nsteps = 10000+1;
% Initial evolution
parm.Nsteps0 = 10000;
% Final simulation time
parm.Tend = parm.h*parm.Nsteps;
% Time grid points
parm.t = 0:parm.h:parm.Tend;

% Number of particles
parm.N = 64;
% Particle equilibrium positions 
parm.x = (0:1:parm.N-1)';
% Forward shift index vector 
parm.idx_p = [2:parm.N 1];
% Backward shift index vector
parm.idx_m = [parm.N 1:parm.N-1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice energy potentials' constants
parm.U0 = 1;
parm.V0 = 0.05;

% Scaled Planck constant
parm.tau = 0.001;
% Rescaling factor
parm.tau_f = sqrt(2*parm.tau);

% Alpha value
parm.alpha = 15;
% J_0/tau value
parm.J0_tau = 0.5*exp(parm.alpha);
% J_0 value
parm.J0 = parm.J0_tau*parm.tau;

% Charge energy parameter values  
parm.E0 = 0;
parm.E0 = 0.375*(2*pi*parm.tau);  % Stationary - gamma=0.4
parm.E0_tau = parm.E0/parm.tau;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Index vectors for implicit midpoint method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parm.I = [1:2*parm.N 1:parm.N 1:parm.N 2:parm.N 1 parm.N+1:2*parm.N...
    parm.N+1:2*parm.N parm.N+2:2*parm.N parm.N+1];
parm.J = [1:2*parm.N parm.N+1:2*parm.N parm.N+2:2*parm.N...
    parm.N+1 parm.N+1:2*parm.N 1:parm.N 2:parm.N 1 1:parm.N];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Index vectors to construct system matrix of charge equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parm.EI = [1:parm.N 1:parm.N 2:parm.N 1];
parm.EJ = [1:parm.N 2:parm.N 1 1:parm.N];
