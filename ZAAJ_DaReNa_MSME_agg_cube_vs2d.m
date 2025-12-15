%MaCaReNa 2.0 Simulation Code

%Yuri
%100x smaller delta T
%3 clusters
%3 protons
%steps 1.5e5
%matrix multiplication
%changing rp by tau
addpath(genpath('/Users/abigailfrey/Dropbox (MIT)/He-Abigail/matlab/Toolbox/'));
addpath(genpath('B:\Dropbox\Data\Matlab\Map'));
addpath(genpath('C:\Dropbox\Data\Matlab\Map'));
addpath(genpath('D:\Dropbox\Data\Matlab\Map'));
addpath 'C:\Users\polto\Desktop\MatLab';
addpath(genpath('C:\Users\Study Room\Dropbox\HW-Lab-Share\Data\Matlab\Map'));%for S1-257
addpath(genpath('C:\Users\nanolab2\Dropbox\Zuzu-HW\Modeling_Simulation\DaReNa')); %zuzu pc

clear all
close all
rng('shuffle');          % Reshuffle the seed
% pause(rand * 0.5);       % Optional random pause
% result = [];             % For any necessary reinitialization

folder = uigetdir;
date = '/agg/cube/20241014';
mkdir(folder, date);    % Change path 
path = strcat(folder,date);

% Write simulation conditions
fid = fopen(strcat(path,'\simulation_conditions.txt'), 'wt' );

header = ['simulation conditions\n'];
simulation_conditions1 = ['dnp = 5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60 nm\n'];
simulation_conditions2 = ['no liposome, C2AB added, 20 arenas and 50 spin-trajectories\n '];

fprintf(fid, simulation_conditions1);
fprintf(fid, simulation_conditions2);
fclose(fid);

% proton_folder = '/proton_trajectories';
% mkdir(path, proton_folder);
% path_proton = strcat(path,proton_folder);
path_proton = 'C:\Dropbox\Protocol_Data_Analysis\MRI modeling simulation\MaCaReNa simulation\free\20210918\proton_trajectories';

% max_saved_p = 500;

%% Tau D Loop

% tau_D_matrix = [10^(-5)*10^(-3); 1/4*10^(-4)*10^(-3); 1/2*10^(-4)*10^(-3);
% 3/4*10^(-4)*10^(-3); 10^(-4)*10^(-3); 1/4*10^(-3)*10^(-3);
% 1/2*10^(-3)*10^(-3); 3/4*10^(-3)*10^(-3); 10^(-3)*10^(-3); 1/4*10^(-2)*10^(-3);
% 1/2*10^(-2)*10^(-3); 10^(-2)*10^(-3); 1/4*10^(-1)*10^(-3)];  % 13 points, tau_D = (rp)^2/D

D = 2.5e-9; %(m^2/s)

dp_input = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60];  %diameter of NP in nanometers
rp_input = dp_input./2; %radius of NP in nanometers

tau_D_matrix = ((rp_input./1e9).^2)./D;

% TE_matrix = [200e-3; 100e-3; 50e-3; 50e-3; 30e-3; 20e-3; 20e-3; 20e-3; 20e-3; 40e-3; 40e-3; 100e-3; 200e-3];  
% TE_matrix = [40e-3; 40e-3; 20e-3; 20e-3; 20e-3; 20e-3; 20e-3; 20e-3; 20e-3; 20e-3; 40e-3; 40e-3; 40e-3];

max_h = length(dp_input);

TE_matrix = 10e-3.*ones(max_h,1);  % TE must be larger than T2 (1/R2)
R2_matrix = ones(max_h,1);
tau_D_matrix_used = ones(max_h,1);

for h = 1:max_h

loop_num_h = sprintf('h loop number = %d', h);
disp(loop_num_h);

%% Variables

tau_D = tau_D_matrix(h);
tau_D_matrix_used (h,1) = tau_D;
% tau_cp = 0.5e-3; %s
D = 2.5e-9; %(m^2/s)
rp = sqrt(tau_D*D)*1e9; %nm
PT_ZDS = 0.8;  %nm
C2AB = 6.5;  %nm, diameter for the calcium-sensitive protein
Bm3h = 6.0; %nm, diameter for the dopamine-sensitive protein
Ca = 0.231;
r_liposome = 50; %nm
dnp = 2*rp; %nm % interparticle distance from center of particle
edge = 100*dnp;
num_of_np_units = 1;
particles_number = 9*num_of_np_units;
sig_digits = 5;

% dist = PT_ZDS + C2AB + dnp;
dist = Bm3h + dnp; % IF C2ABs serve as bridges for NPs

wr = 2.36e7; %rad/s % rms angular frequency for magnetite (2.36e7) or maghemite (2.36e7/92*76 = 1.95e7)
gamma = 2.6752e8; % rad s^-1 T^-1 % gyromagnetic ratio of a proton

TE_max = TE_matrix(h); %s
tau_cp0 = TE_max/100;
tau_cp = max(tau_cp0,100*tau_D);   %tau_cp should satisfy the condition: tau_cp >> tau_D

% delta_t0 = (edge*1e-9)^2/(24*D)/((particles_number)^2/3); %s
delta_t0 = (dnp*1e-9)^2/(6*D);
% delta_t0 = (2*dnp*1e-9)^2/(6*D);        % 4 times faster than normal by doing 4 times bigger walking steps, but the accuracy is smaller
delta_t = min(delta_t0, 0.2*tau_cp); %s

max_s = 20;        % number of different cluster positions
max_p = 50;     %number of protons
steps = round(TE_max/delta_t);      %number of steps
TE_points = floor(TE_max/tau_cp/2);
I_norm = zeros(max_s,TE_points);
phi_all = zeros(max_p,TE_points);


%% Make random particles

for s=1:max_s
    
% particles_number = round((182.9/(rp*2+PZ_ZDS+C2AB/2))^3*(4*pi/3/8));     %cannot divide one radius by another directly. 182.9 is the diameter of aggregate given by DLS. But this number is for close-packed cube, so it has to be corrected by the volume ratio of (sphere/cube) who has the same radius.
if h <= 4
    loop_num_s = sprintf('s loop number = %d', s);
    disp(loop_num_s);
end

% size of box for 0.1 mM Fe concentration = (840.64^3*(6/27))^(1/3) = 509.1814 nm egdes

%edge = 509.1814;

%% Make Liposome/Particle Coordinates

M = [edge/2 edge/2 edge/2; edge/2 edge/2 edge*3/2; edge/2 edge*3/2 edge*3/2; edge/2 edge*3/2 edge/2; edge*3/2 edge/2 edge/2; edge*3/2 edge/2 edge*3/2; edge*3/2 edge*3/2 edge*3/2; 
    edge*3/2 edge*3/2 edge/2; edge*3/2 edge/2 edge*5/2; edge*3/2 edge*3/2 edge*5/2; edge/2 edge/2 edge*5/2; edge/2 edge*3/2 edge*5/2; edge/2 edge*5/2 edge*3/2; edge/2 edge*5/2 edge/2;
    edge*3/2 edge*5/2 edge*3/2; edge*3/2 edge*5/2 edge/2; edge*3/2 edge*5/2 edge*5/2; edge/2 edge*5/2 edge*5/2; edge*5/2 edge/2 edge/2; edge*5/2 edge/2 edge*3/2; edge*5/2 edge*3/2 edge*3/2;
    edge*5/2 edge*3/2 edge/2; edge*5/2 edge/2 edge*5/2; edge*5/2 edge*3/2 edge*5/2; edge*5/2 edge*5/2 edge*3/2; edge*5/2 edge*5/2 edge/2; edge*5/2 edge*5/2 edge*5/2];  % 27*3 matrix 


%% Plot Cube

%% Plot NPs 

particles = zeros(particles_number,3);

x1_np = rand*edge;
y1_np = rand*edge;
z1_np = rand*edge;

NP1 = [x1_np, y1_np, z1_np];

cube_offsets = dist/2 .* [-1 -1 -1;  % Compute the 8 vertices of the cube around NP1
                 -1 -1  1;  % Vertices are L/2 from NP1 in each direction
                 -1  1 -1;
                 -1  1  1;
                  1 -1 -1;
                  1 -1  1;
                  1  1 -1;
                  1  1  1];

cube_vertices = NP1 + cube_offsets; % Calculating the coordinates for NP2-8

plotted_1 = cat(1, NP1, cube_vertices);

plotted = plotted_1;

%% Make other clusters

for l = 2:27
    
    diff = plotted_1 - M(1,:);
    matrix = ones(particles_number,3).*M(l,:);
    plotted(particles_number*l-(particles_number-1):particles_number*l,1:3) = matrix+diff;

end 

np_correction = 1*edge;

plotted = plotted - np_correction;  % center the 27 cubes at origin

%% Plot Configuration

if h==max_h && s==max_s

figure;

plotcube([edge edge edge], [0 0 0],.1,[1 0 0]); 
hold on;

% Plot particles
plot3 (plotted(:,1), plotted(:,2), plotted(:,3), 'r.','MarkerSize',12);

xlabel('nanometers');
ylabel('nanometers');
zlabel('nanometers');
hold off;

caption = sprintf('In this configuration, NP radius = %f nm, liposome radius = %d nm, tau_c_p = %f s', rp, r_liposome, tau_cp);
title(caption,'FontSize', 14);

saveas(gcf, [strcat(path,'/np_dist.fig')]);
pause(5);
close;

end

%% Make Protons

for p=1:max_p
    
n_prot=1; %this is really always 1, not proton number
x_protons=rand(1)*edge;
y_protons=rand(1)*edge; %
z_protons=rand(1)*edge; %
protons = [x_protons y_protons z_protons];

%plot3 (x_protons,y_protons,z_protons,'k.') 

%% Protons Perform Random Walks

sigma = sqrt(6*D*delta_t);      %in the unit of m
walk = sigma/1e-9;       %in the unit of nm

rand_proton = 1;
load_proton = 0;

if rand_proton

% proton_direction_point = [rand(steps,1)*edge rand(steps,1)*edge rand(steps,1)*edge];  
proton_walks = ones(steps,3);
proton_walks(1,:) = protons(1,:);

alpha_angle = rand(steps,1)*2*pi;
beta_angle = rand(steps,1)*2*pi;

for j=2:steps
    proton_walks(j,1) = proton_walks(j-1,1)+ walk*cos(alpha_angle(j))*cos(beta_angle(j));
    proton_walks(j,2) = proton_walks(j-1,2)+ walk*cos(alpha_angle(j))*sin(beta_angle(j));
    proton_walks(j,3) = proton_walks(j-1,3)+ walk*sin(alpha_angle(j));
end


%% Periodic Boundary Conditions (PBC)

for k=1:3
parfor q=1:steps
    
   if proton_walks(q,k) < 0
       n_pbc = ceil(-proton_walks(q,k)/edge);
       proton_walks(q,k) = proton_walks(q,k) + n_pbc*edge;
   end
   
   if proton_walks(q,k) > edge
       n_pbc = floor(proton_walks(q,k)/edge);
       proton_walks(q,k) = proton_walks(q,k) - n_pbc*edge;
   end   
   
end
end

proton_walks = round(proton_walks, sig_digits, 'significant');


elseif load_proton
    
%     trajectory_num = randi([1 max_saved_p]);
    trajectory_num = p;
    cd(path_proton);
    proton_file = sprintf('proton_walks_%dthH_%dthP', h, trajectory_num);
    load(proton_file);
end

if h==max_h && s==max_s && p==max_p

dimension_1 = round((steps)^0.25); %avg a number of x/y/z elements where the number = dimension_1
dimension_2 = 1; %no cross-avg among x, y, and z 
proton_walks_avg = blockproc(proton_walks, [dimension_1 dimension_2], @(x) mean(x.data, 'all'));  %need image processing toolbox
figure;
plot3(proton_walks_avg(:,1),proton_walks_avg(:,2),proton_walks_avg(:,3));
xlabel('nanometers');
ylabel('nanometers');
zlabel('nanometers');
saveas(gcf, [strcat(path,'/proton_walk.fig')]);
close;

end

%% Magnetic Field

a = [0,0,0]; %edge of our cube
b = [0,0,edge];

[num c] = size(plotted);

Bz = zeros(n_prot*steps,1);

parfor j = 1:n_prot*steps   %parfor

% use matrix multiplication to calc particle->proton distance

vector_NP_to_P = repmat(proton_walks(j,:),num,1) - plotted;
vector_for_theta = vector_NP_to_P;

d_P = vecnorm(vector_NP_to_P,2,2)';

% for any d_P < NP radius, set that d_P = NP radius (aka impenetrable NPs)

d_P(d_P<rp) = rp;

% calc theta, which is the angle between the z-axis and where Bz is evaluated 

u = repmat(b,size(vector_for_theta,1),1);

denominator = vecnorm(b,2,2).*vecnorm(vector_for_theta,2,2);
denominator_invert = denominator.^(-1);

costheta = (dot(u,vector_for_theta,2).*(denominator_invert))';
% theta = real(acos(costheta));      

%calc Bz

Bz_P = sqrt(5/4)*(rp^3*wr)./(gamma*d_P.^3).*(3.*(costheta).^2-1); 

Bz(j,1) = sum(Bz_P);
end


psi = gamma*Bz*delta_t;  

%calculate the sum of psi between inversions

num_inversion = round(TE_max/tau_cp);

for j=1:num_inversion
    start = floor((steps/num_inversion)*(j-1)+1);
    final = floor((steps/num_inversion)*j);
    eval(sprintf('psi_%d = psi(%d:%d);',j,start,final));
    eval(sprintf('psisum_%d = sum(psi_%d);',j,j));
end

coefficient = ones(1,num_inversion);

for k=1:ceil(num_inversion/2)       % MSME sequence 
    coefficient(1,(2*k-1))=(-1)^(k-1);
    coefficient(1,(2*k))=(-1)^(k);
end

psiall = zeros(1,num_inversion);

for i=1:num_inversion
    eval(sprintf('psiall(1,%d) = (psisum_%d).*(coefficient(1,%d));',i,i,i));
end

phi = ones(1,TE_points);

for j=1:TE_points
    jj = 2*j;
    eval(sprintf('phi(1,%d) = sum(psiall(1,1:%d)) ;',j,jj));    
end

eval(sprintf('phi_all(%d,:) = phi(1,:);',p));

end    %end of protons loop

I_all = cos(phi_all);

I_avg = nanmean(I_all,1);  % average the signals from all protons

eval(sprintf('I_norm(%d,:) = I_avg(1,:);',s));

end    %end of clusters loop

%% Save data and fit
I_norm_final = nanmean(I_norm,1);       % average the signals from all clusters
minus_ln_signal = -log(I_norm_final);

TE = [(TE_max/TE_points):(TE_max/TE_points):TE_max].*1e3; %ms

W = polyfit(TE,minus_ln_signal,1);
t_plot = min(TE):(max(TE)-min(TE))/300:max(TE);
y = minus_ln_signal;
yfit = W(1)*TE+W(2);
yplot = W(1)*t_plot+W(2);

Rsq = 1 - sum(rmmissing((y - yfit).^2))/sum(rmmissing((y - mean(y)).^2));
R2 = W(1)*1e3;      % in the unit of s^-1

R2_matrix(h,1) = R2;

figure;
plot(TE,minus_ln_signal,'bo',t_plot,yplot,'r');
legend(' Results',' Fitting','FontSize', 16);
xlabel('TE time (ms)');
ylabel('R2 fitting (-ln[I(t)/Io)]'); 
caption = sprintf('In this simulation, there are %d protons, %d arenas, and %d steps. R_2 = %f s^-^1 and Rsq = %f', max_p, max_s, steps, R2, Rsq);
title(caption,'FontSize', 14);
fullscreen = gcf;
fullscreen.WindowState = 'maximized';
cd(path);
saveas(gcf, sprintf('R2_parameters_%dth_tauD.fig',h));
close;

save(strcat(path,'/R2_matrix'),'R2_matrix');

end 

R2_matrix_real = real(R2_matrix);
save(strcat(path,'/R2_matrix_real'),'R2_matrix_real');
save(strcat(path,'/tau_D_matrix'),'tau_D_matrix');

% figure;
% loglog(tau_D_matrix(:)*10^3, R2_matrix_real(:), 'ko')
% xlim([10^(-6) 1]);
% xlabel('Tau_D (ms)');
% ylabel('R2 (s^-1)'); 
% % set(gca, 'XScale', 'log10')
% % set(gca, 'YScale', 'log10')
% saveas(gcf, [strcat(path,'/trend_plot_log.fig')]);
% close;

figure;
plot(dp_input(:), R2_matrix_real(:), 'ko');
% xlim([0 11]);
xlabel('NP diameter (nm)');
ylabel('R2 (s^-^1)'); 
caption = sprintf('simulation result of %d arenas and %d spin-trajectories', max_s, max_p);
title(caption,'FontSize', 14);
% set(gca, 'XScale', 'log10')
% set(gca, 'YScale', 'log10')
saveas(gcf, [strcat(path,'/trend_plot.fig')]);
close;


%include fit curve
figure;
f = fit(log10(tau_D_matrix(:)*10^3), log10(R2_matrix_real(:)), 'poly2');
plot(f, log10(tau_D_matrix(:)*10^3), log10(R2_matrix_real(:)), 'ko')
% xlim([-6 -4.5]);
xlabel('log Tau_D (ms)','FontSize', 20);
ylabel('log R2 (s^-1)','FontSize', 20); 

% spin_traject = max_s*max_p;
caption = sprintf('Aggregated NPs: simulation result of %d arenas and %d spin-trajectories', max_s, max_p);
title(caption,'FontSize', 14);
set(gca, 'XScale')
set(gca, 'YScale')
saveas(gcf, [strcat(path,'/trend_plot_log_fitting.fig')]);