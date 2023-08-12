function vars = TimeInt_Split_Q(~,vars,h)

% Exact integration of piece Q with time step h
vars.u = vars.u + h*vars.p;

end