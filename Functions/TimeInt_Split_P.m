function vars = TimeInt_Split_P(parm,vars,h)

% Exact integration of piece P with time step h

% Compute forces acting on lattice particles
G = Comp_Lattice_OnsiteForce(parm,vars);
F = Comp_Lattice_InteractionForce(parm,vars);

% Update p with time step h
vars.p = vars.p + h*(F + G);

end