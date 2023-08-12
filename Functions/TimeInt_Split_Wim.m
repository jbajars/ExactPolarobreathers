function vars = TimeInt_Split_Wim(parm,vars,h)

% Implicit midpoint integration of piece W with time step h
a = vars.a;
b = vars.b;

% Compute new a and b values
vars = TimeInt_ImMidpStepWithoutEn_ChargeEqs(parm,vars,h);
A = vars.a;
B = vars.b;

% Compute charge-lattice force G
vars.a = (A+a)/2;
vars.b = (B+b)/2;
G_W = Comp_ChargeLattice_GW_Force(parm,vars);

% Update p with time step h
vars.p = vars.p + h*G_W;

% New a and b values
vars.a = A;
vars.b = B;

end