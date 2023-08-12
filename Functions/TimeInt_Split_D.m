function vars = TimeInt_Split_D(parm,vars,h)

% Exact integration of piece D with time step h

% Exponent of En_tau
Exp_En_tau = exp(-parm.E0_tau*h*1i); 

% Exact integration 
z = vars.a + 1i*vars.b;
z = Exp_En_tau.*z;
vars.a = real(z); 
vars.b = imag(z);

end