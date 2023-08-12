%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Charge transfer by discrete breathers in 1D model
% written by Janis Bajars, February 2023
% Plot spectra of approximate and exact polarobreathers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clearvars
clc
set(0,'DefaultAxesFontSize',16)
set(groot,'DefaultAxesTickLabelInterpreter','latex'); 
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

disp('Plot spectra!')
tic;

% Load colormap
load ColorMapCM.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load simulation data and plot XTFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution. Exact: true, Approximate: false
Sol = true;
% Motion. Moving: true, Stationary: false
Mot = true;

if Sol
    if Mot
        % Exact moving polarobreather solution
        load('SavedData/Exact_SimData_Mov_G06_E0shift_N64Tau001.mat')
    else
        % Exact stationary polarobreather solution
        load('SavedData/Exact_SimData_Stat_G04_E0shift_N64Tau001.mat')
    end
else
    if Mot
        % Approximate moving polarobreather solution
        load('SavedData/Approx_SimData_Mov_G06_E0shift_N64Tau001.mat')
    else
        % Approximate stationary polarobreather solution
        load('SavedData/Approx_SimData_Stat_G04_E0shift_N64Tau001.mat')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice XTFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
set(gcf, 'Position',  [300, 500, 1300, 400])
t = tiledlayout(1,3,'TileSpacing','Compact');
nexttile
hold on

% Lattice dispersion relation
omega2 = (2*pi*sqrt(parm.U0))^2;
c2 = 72*parm.V0;
q = -pi:0.01:pi;
dr = omega2 + 4*c2*sin(q/2).^2;
plot(q/pi,sqrt(dr)/2/pi,'-','color',[0 0 0]+0.75,'linewidth',1.5)
plot(q/pi,-sqrt(dr)/2/pi,'-','color',[0 0 0]+0.75,'linewidth',1.5)

% Charge dispersion relation
q = -pi:0.01:pi;
wq = -2*parm.J0_tau*exp(-parm.alpha)*cos(q);
plot(q/pi,wq/2/pi,'--','color',[0 0 0]+0.75,'linewidth',1.5)

% XTFT of u
[~, nt] = size(U);
dt = parm.h;
ff = (-nt/2:nt/2-1)/nt/dt; 
qq = 2*(-parm.N/2:parm.N/2-1)/parm.N;
Uxt = fft2(U);
Uxt = fftshift(Uxt); 
Uxt = abs(Uxt).^2;   
fft_norm = max(max(Uxt));
Uxt = Uxt./fft_norm;

% Plot XTFT
colormap(cm)
contour(qq,-ff,Uxt',100)
xlabel('$q/\pi$')
ylabel('$\omega/2\pi$')
title('Normalized $|XTFT|^2$ of $u_n$','fontsize',16)
axis([-1 1 -1.5 1.5])
xticks([-1 -0.5 0 0.5 1])
yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
xticklabels({'$-$1','$-$0.5','0','0.5','1'})
yticklabels({'$-$1.5','$-$1','$-$0.5','0','0.5','1','1.5'})
axis on
box on
grid on
caxis([0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Charge probability XTFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
hold on

% Lattice dispersion relation
omega2 = (2*pi*sqrt(parm.U0))^2;
c2 = 72*parm.V0;
q = -pi:0.01:pi;
dr = omega2 + 4*c2*sin(q/2).^2;
plot(q/pi,sqrt(dr)/2/pi,'-','color',[0 0 0]+0.75,'linewidth',1.5)
plot(q/pi,-sqrt(dr)/2/pi,'-','color',[0 0 0]+0.75,'linewidth',1.5)

% Charge dispersion relation
q = -pi:0.01:pi;
wq = -2*parm.J0_tau*exp(-parm.alpha)*cos(q);
plot(q/pi,wq/2/pi,'--','color',[0 0 0]+0.75,'linewidth',1.5)

% XTFT of |c_n|^2
[~, nt] = size(A);
dt = parm.h;
ff = (-nt/2:nt/2-1)/nt/dt; 
kk = 2*(-parm.N/2:parm.N/2-1)/parm.N;
Cxt = fft2(abs((A+1i*B)/sqrt(2*parm.tau)).^2);
Cxt = fftshift(Cxt);
Cxt = abs(Cxt).^2;   
fft_norm = max(max(Cxt));
Cxt = Cxt./fft_norm;

% Plot XTFT
colormap(cm)
contour(kk,-ff,Cxt',100)
xlabel('$q/\pi$')
title('Normalized $|XTFT|^2$ of $|c_n|^2$','fontsize',16)
axis([-1 1 -1.5 1.5])
xticks([-1 -0.5 0 0.5 1])
yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
xticklabels({'$-$1','$-$0.5','0','0.5','1'})
yticklabels({'$-$1.5','$-$1','$-$0.5','0','0.5','1','1.5'})
axis on
box on
grid on
caxis([0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Charge amplitude XTFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
hold on

% Lattice dispersion relation
omega2 = (2*pi*sqrt(parm.U0))^2;
c2 = 72*parm.V0;
q = -pi:0.01:pi;
dr = omega2 + 4*c2*sin(q/2).^2;
plot(q/pi,sqrt(dr)/2/pi,'-','color',[0 0 0]+0.75,'linewidth',1.5)
plot(q/pi,-sqrt(dr)/2/pi,'-','color',[0 0 0]+0.75,'linewidth',1.5)

% Charge dispersion relation
q = -pi:0.01:pi;
wq = -2*parm.J0_tau*exp(-parm.alpha)*cos(q);
plot(q/pi,wq/2/pi,'--','color',[0 0 0]+0.75,'linewidth',1.5)

% XTFT of c_n
[~, nt] = size(A);
dt = parm.h;
ff = (-nt/2:nt/2-1)/nt/dt; 
kk = 2*(-parm.N/2:parm.N/2-1)/parm.N;
Cxt = fft2((A+1i*B)/sqrt(2*parm.tau));
Cxt = fftshift(Cxt);
Cxt = abs(Cxt).^2;   
fft_norm = max(max(Cxt));
Cxt = Cxt./fft_norm;

% Plot XTFT
colormap(cm)
contour(kk,-ff-parm.E0_tau/2/pi,Cxt',100)
xlabel('$q/\pi$')
title('Normalized $|XTFT|^2$ of $c_n$','fontsize',16)
axis([-1 1 -1.5 1.5])
xticks([-1 -0.5 0 0.5 1])
yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
xticklabels({'$-$1','$-$0.5','0','0.5','1'})
yticklabels({'$-$1.5','$-$1','$-$0.5','0','0.5','1','1.5'})
axis on
box on
grid on
caxis([0 1])

% Add colorbar
figure(1)
c=colorbar;
set(c,'TickLabelInterpreter','latex')

% Export and save figure
if Sol
    if Mot
        % Exact moving polarobreather solution
        exportgraphics(gcf,'Figures/mov_exact_spectrum.png',...
            'Resolution',300)
    else
        % Exact stationary polarobreather solution
        exportgraphics(gcf,'Figures/stat_exact_spectrum.png',...
            'Resolution',300)
    end
else
    if Mot
        % Approximate moving polarobreather solution
        exportgraphics(gcf,'Figures/mov_approx_spectrum.png',...
            'Resolution',300)
    else
        % Approximate stationary polarobreather solution
        exportgraphics(gcf,'Figures/stat_approx_spectrum.png',...
            'Resolution',300)
    end
end

toc;