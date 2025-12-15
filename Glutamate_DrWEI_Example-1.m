% Glutamate MaCaReNa 2.0 Simulation Code
%
% Diya
% Last updated: 2025-12-08
%
% 100x smaller delta T
% 3 clusters
% 3 protons
% steps 1.5e5
% matrix multiplication
% changing rp by tau

function [delta_R2_ratio, Rsq_agg, Rsq_free, R2_matrix_real_agg, R2_matrix_real_free] = Glutamate_DrWEI_Example(a_i)

%MAY ADD THIS BACK LATER
%addpath(genpath('/Users/abigailfrey/Dropbox (MIT)/He-Abigail/matlab/Toolbox/'));
%addpath(genpath('B:\Dropbox\Data\Matlab\Map'));
%addpath(genpath('C:\Dropbox\Data\Matlab\Map'));
%addpath(genpath('D:\Dropbox\Data\Matlab\Map'));
% addpath 'C:\Users\polto\Desktop\MatLab';
%addpath(genpath('C:\Users\Study Room\Dropbox\HW-Lab-Share\Data\Matlab\Map'));  %for S1-257

dp_input = [a_i];     %diameter of NP in nanometers

rng('shuffle');       % Reshuffle the seed

%MAY ADD THIS BACK LATER
%path_proton = 'C:\Dropbox\Protocol_Data_Analysis\MRI modeling simulation\MaCaReNa simulation\free\20210918\proton_trajectories';

%% Tau D Loop

D = 2.5e-9; %(m^2/s)

Cl_coordinates = [-24.265907, 8.822892, 2.053962;      %binding site #1
-10.597951, 1.409974, -33.878379;      %binding site #2
-5.757979, 6.389922, 1.764968;         %binding site #3
6.465975, 11.915850, -14.044743;       %binding site #4
-3.453985, 25.533678, -30.517441;      %binding site #5
-6.245980, 24.110701, -6.942873;
1.176995, 18.230771, -14.854728;
2.691990, 19.332756, -21.732602;]./10;  %the coordinates of 8 chloride anions in nanometers. 
                                        %They are also the 8 binding sites to NPs. 

rp_input = dp_input./2; %radius of NP in nanometers

tau_D_matrix = ((rp_input./1e9).^2)./D;

max_h = length(dp_input);

TE_matrix = 10e-3.*ones(max_h,1);  % TE must be larger than T2 (1/R2)
R2_matrix = ones(max_h,1);
tau_D_matrix_used = ones(max_h,1);

for h = 1:max_h

    %% Variables

    tau_D = tau_D_matrix(h);
    tau_D_matrix_used (h,1) = tau_D;
    D = 2.5e-9; %(m^2/s)
    rp = sqrt(tau_D*D)*1e9; %nm
    PT_ZDS = 0.8;  %nm
    C2AB = 6.5;    %nm
    Ca = 0.231;
    r_liposome = 50; %nm
    dnp = 2*rp; %nm % interparticle distance from center of particle
    edge = 100*dnp;
    num_of_np_units = 1;
    particles_number = 8;
    sig_digits = 5;

    dist = C2AB + dnp; % IF C2ABs serve as bridges for NPs

    dist_1to2 = dnp + norm(Cl_coordinates(1,:) - Cl_coordinates(2,:));
    dist_2to3 = dnp + norm(Cl_coordinates(2,:) - Cl_coordinates(3,:));
    dist_3to4 = dnp + norm(Cl_coordinates(3,:) - Cl_coordinates(4,:));
    dist_4to5 = dnp + norm(Cl_coordinates(4,:) - Cl_coordinates(5,:));

    wr = 2.36e7; %rad/s 
    gamma = 2.6752e8; % rad s^-1 T^-1 % gyromagnetic ratio of a proton

    TE_max = TE_matrix(h); %s
    tau_cp0 = TE_max/100;
    tau_cp = max(tau_cp0,100*tau_D);   %tau_cp should satisfy the condition: tau_cp >> tau_D

    delta_t0 = (dnp*1e-9)^2/(6*D);
    delta_t = min(delta_t0, 0.2*tau_cp); %s

    max_s = 1;        % number of different cluster positions
    max_p = 15;       %number of protons
    steps = round(TE_max/delta_t);      %number of steps
    TE_points = floor(TE_max/tau_cp/2);
    I_norm = zeros(max_s,TE_points);
    phi_all = zeros(max_p,TE_points);


    %% Make random particles

    for s=1:max_s

        %% Make Liposome/Particle Coordinates

        M = [edge/2 edge/2 edge/2; edge/2 edge/2 edge*3/2; edge/2 edge*3/2 edge*3/2; edge/2 edge*3/2 edge/2; edge*3/2 edge/2 edge/2; edge*3/2 edge/2 edge*3/2; edge*3/2 edge*3/2 edge*3/2; 
            edge*3/2 edge*3/2 edge/2; edge*3/2 edge/2 edge*5/2; edge*3/2 edge*3/2 edge*5/2; edge/2 edge/2 edge*5/2; edge/2 edge*3/2 edge*5/2; edge/2 edge*5/2 edge*3/2; edge/2 edge*5/2 edge/2;
            edge*3/2 edge*5/2 edge*3/2; edge*3/2 edge*5/2 edge/2; edge*3/2 edge*5/2 edge*5/2; edge/2 edge*5/2 edge*5/2; edge*5/2 edge/2 edge/2; edge*5/2 edge/2 edge*3/2; edge*5/2 edge*3/2 edge*3/2;
            edge*5/2 edge*3/2 edge/2; edge*5/2 edge/2 edge*5/2; edge*5/2 edge*3/2 edge*5/2; edge*5/2 edge*5/2 edge*3/2; edge*5/2 edge*5/2 edge/2; edge*5/2 edge*5/2 edge*5/2];  % 27*3 matrix 

        %% Plot NPs  – cube geometry

        % Random centroid for NP cluster
        x1_np = rand*edge;
        y1_np = rand*edge;
        z1_np = rand*edge;
        centroid = [x1_np, y1_np, z1_np];

        % 8-point cube offsets
        cube_offsets = [-0.33795, -0.33795, -0.33795;
                        -0.33795, -0.33795,  0.33795;
                        -0.33795,  0.33795, -0.33795;
                        -0.33795,  0.33795,  0.33795;
                         0.33795, -0.33795, -0.33795;
                         0.33795, -0.33795,  0.33795;
                         0.33795,  0.33795, -0.33795;
                         0.33795,  0.33795,  0.33795];

        cube_points = centroid + cube_offsets;

        % Random rotation angles
        gamma_angle = rand * 2*pi;
        delta_angle = rand * 2*pi;

        % Rotation matrices
        RotationZ = [cos(gamma_angle), -sin(gamma_angle), 0;
                     sin(gamma_angle),  cos(gamma_angle), 0;
                     0,                 0,                1];

        RotationY = [cos(delta_angle),  0, sin(delta_angle);
                     0,                1, 0;
                    -sin(delta_angle), 0, cos(delta_angle)];

        % Apply rotation
        cube_final = (RotationY * RotationZ * cube_points')';

        % Final aggregated geometry (8 NPs)
        plotted_1 = cube_final;
        plotted   = plotted_1;

        %% Make other clusters

        for l = 2:27
            diff = plotted_1 - M(1,:);
            matrix = ones(particles_number,3).*M(l,:);
            plotted(particles_number*l-(particles_number-1):particles_number*l,1:3) = matrix+diff;
        end 

        np_correction = 1*edge;
        plotted = plotted - np_correction;  % center the 27 cubes at origin

        %% Make Protons

        for p=1:max_p

            n_prot=1;
            x_protons=rand(1)*edge;
            y_protons=rand(1)*edge;
            z_protons=rand(1)*edge;
            protons = [x_protons y_protons z_protons];

            %% Protons Perform Random Walks

            sigma = sqrt(6*D*delta_t);      %in the unit of m
            walk = sigma/1e-9;             %in the unit of nm

            rand_proton = 1;
            load_proton = 0;

            if rand_proton

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

                trajectory_num = p;
                cd(path_proton);
                proton_file = sprintf('proton_walks_%dthH_%dthP', h, trajectory_num);
                load(proton_file);

            end

            %% Magnetic Field

            a = [0,0,0]; %edge of our cube
            b = [0,0,edge];

            [num, c] = size(plotted); %#ok<ASGLU>

            Bz = zeros(n_prot*steps,1);

            parfor j = 1:n_prot*steps

                vector_NP_to_P = repmat(proton_walks(j,:),num,1) - plotted;
                vector_for_theta = vector_NP_to_P;

                d_P = vecnorm(vector_NP_to_P,2,2)';

                d_P(d_P<rp) = rp;

                u = repmat(b,size(vector_for_theta,1),1);

                denominator = vecnorm(b,2,2).*vecnorm(vector_for_theta,2,2);
                denominator_invert = denominator.^(-1);

                costheta = (dot(u,vector_for_theta,2).*(denominator_invert))';

                Bz_P = sqrt(5/4)*(rp^3*wr)./(gamma*d_P.^3).*(3.*(costheta).^2-1); 

                Bz(j,1) = sum(Bz_P);
            end

            psi = gamma*Bz*delta_t;  

            %% calculate the sum of psi between inversions

            num_inversion = round(TE_max/tau_cp);

            for j=1:num_inversion
                start = floor((steps/num_inversion)*(j-1)+1);
                final = floor((steps/num_inversion)*j);
                eval(sprintf('psi_%d = psi(%d:%d);',j,start,final));
                eval(sprintf('psisum_%d = sum(psi_%d);',j,j));
            end

            coefficient = ones(1,num_inversion);

            for kk=1:ceil(num_inversion/2)       % MSME sequence 
                coefficient(1,(2*kk-1))=(-1)^(kk-1);
                coefficient(1,(2*kk))=(-1)^(kk);
            end

            psiall = zeros(1,num_inversion);

            for i2=1:num_inversion
                eval(sprintf('psiall(1,%d) = (psisum_%d).*(coefficient(1,%d));',i2,i2,i2));
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

    %% Save data and fit (aggregated)

    I_norm_final = nanmean(I_norm,1);       % average the signals from all clusters
    minus_ln_signal = -log(I_norm_final);

    TE = [(TE_max/TE_points):(TE_max/TE_points):TE_max].*1e3; %ms

    W = polyfit(TE,minus_ln_signal,1);
    t_plot = min(TE):(max(TE)-min(TE))/300:max(TE);
    y = minus_ln_signal;
    yfit = W(1)*TE+W(2);
    yplot = W(1)*t_plot+W(2);

    Rsq_agg = 1 - sum(rmmissing((y - yfit).^2))/sum(rmmissing((y - mean(y)).^2));
    R2 = W(1)*1e3;      % in the unit of s^-1

    R2_matrix(h,1) = R2;

end 

R2_matrix_real_agg = real(R2_matrix);

%% free state

rng('shuffle');          % Reshuffle the seed

path_proton = 'C:\Dropbox\Protocol_Data_Analysis\MRI modeling simulation\MaCaReNa simulation\free\20210918\proton_trajectories';

D = 2.5e-9; %(m^2/s)

rp_input = dp_input./2; %radius of NP in nanometers
tau_D_matrix = ((rp_input./1e9).^2)./D;

max_h = length(dp_input);

TE_matrix = 10e-3.*ones(max_h,1);  % TE must be larger than T2 (1/R2)
R2_matrix = ones(max_h,1);
tau_D_matrix_used = ones(max_h,1);

for h = 1:max_h

    %% Variables

    tau_D = tau_D_matrix(h);
    tau_D_matrix_used (h,1) = tau_D;
    D = 2.5e-9; %(m^2/s)
    rp = sqrt(tau_D*D)*1e9; %nm
    PT_ZDS = 0.8;  %nm
    C2AB = 6.5;    %nm
    Ca = 0.231; 
    r_liposome = 50; %nm
    dnp = 2*rp; %nm 
    edge = 100*dnp;
    num_of_np_units = 1;
    particles_number = 5*num_of_np_units;
    sig_digits = 5;

    dist = C2AB + dnp; % C2ABs serve as bridges for NPs

    wr = 2.36e7; %rad/s 
    gamma = 2.6752e8; % rad s^-1 T^-1

    TE_max = TE_matrix(h); %s
    tau_cp0 = TE_max/100;
    tau_cp = max(tau_cp0,100*tau_D);   %tau_cp should satisfy the condition: tau_cp >> tau_D

    delta_t0 = (dnp*1e-9)^2/(6*D);
    delta_t = min(delta_t0, 0.2*tau_cp); %s

    max_s = 1;
    max_p = 15;
    steps = round(TE_max/delta_t);
    TE_points = floor(TE_max/tau_cp/2);
    I_norm = zeros(max_s,TE_points);
    phi_all = zeros(max_p,TE_points);

    %% Make random particles

    for s=1:max_s

        M = [edge/2 edge/2 edge/2; edge/2 edge/2 edge*3/2; edge/2 edge*3/2 edge*3/2; edge/2 edge*3/2 edge/2; edge*3/2 edge/2 edge/2; edge*3/2 edge/2 edge*3/2; edge*3/2 edge*3/2 edge*3/2; 
            edge*3/2 edge*3/2 edge/2; edge*3/2 edge/2 edge*5/2; edge*3/2 edge*3/2 edge*5/2; edge/2 edge/2 edge*5/2; edge/2 edge*3/2 edge*5/2; edge/2 edge*5/2 edge*3/2; edge/2 edge*5/2 edge/2;
            edge*3/2 edge*5/2 edge*3/2; edge*3/2 edge*5/2 edge/2; edge*3/2 edge*5/2 edge*5/2; edge/2 edge*5/2 edge*5/2; edge*5/2 edge/2 edge/2; edge*5/2 edge/2 edge*3/2; edge*5/2 edge*3/2 edge*3/2;
            edge*5/2 edge*3/2 edge/2; edge*5/2 edge/2 edge*5/2; edge*5/2 edge*3/2 edge*5/2; edge*5/2 edge*5/2 edge*3/2; edge*5/2 edge*5/2 edge/2; edge*5/2 edge*5/2 edge*5/2];  % 27*3 matrix 

        %% Free NPs

        x_particles=rand(particles_number,1)*edge;
        y_particles=rand(particles_number,1)*edge;
        z_particles=rand(particles_number,1)*edge;

        particles = [x_particles y_particles z_particles];
        plotted_1 = particles;
        plotted = plotted_1;

        for l = 2:27
            diff = plotted_1 - M(1,:);
            matrix = ones(particles_number,3).*M(l,:);
            plotted(particles_number*l-(particles_number-1):particles_number*l,1:3) = matrix+diff;
        end 

        np_correction = 1*edge;
        plotted = plotted - np_correction;

        %% Make Protons

        for p=1:max_p

            n_prot=1;
            x_protons=rand(1)*edge;
            y_protons=rand(1)*edge;
            z_protons=rand(1)*edge;
            protons = [x_protons y_protons z_protons];

            %% Protons Perform Random Walks

            sigma = sqrt(6*D*delta_t);      %in the unit of m
            walk = sigma/1e-9;             %in the unit of nm

            rand_proton = 1;
            load_proton = 0;

            if rand_proton

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

                trajectory_num = p;
                cd(path_proton);
                proton_file = sprintf('proton_walks_%dthH_%dthP', h, trajectory_num);
                load(proton_file);
            end

            %% Magnetic Field

            a = [0,0,0]; 
            b = [0,0,edge];

            [num, c] = size(plotted); %#ok<ASGLU>

            Bz = zeros(n_prot*steps,1);

            parfor j = 1:n_prot*steps

                vector_NP_to_P = repmat(proton_walks(j,:),num,1) - plotted;
                vector_for_theta = vector_NP_to_P;

                d_P = vecnorm(vector_NP_to_P,2,2)';

                d_P(d_P<rp) = rp;

                u = repmat(b,size(vector_for_theta,1),1);

                denominator = vecnorm(b,2,2).*vecnorm(vector_for_theta,2,2);
                denominator_invert = denominator.^(-1);

                costheta = (dot(u,vector_for_theta,2).*(denominator_invert))';

                Bz_P = sqrt(5/4)*(rp^3*wr)./(gamma*d_P.^3).*(3.*(costheta).^2-1); 

                Bz(j,1) = sum(Bz_P);
            end

            psi = gamma*Bz*delta_t;  

            num_inversion = round(TE_max/tau_cp);

            for j=1:num_inversion
                start = floor((steps/num_inversion)*(j-1)+1);
                final = floor((steps/num_inversion)*j);
                eval(sprintf('psi_%d = psi(%d:%d);',j,start,final));
                eval(sprintf('psisum_%d = sum(psi_%d);',j,j));
            end

            coefficient = ones(1,num_inversion);

            for kk=1:ceil(num_inversion/2)       % MSME sequence 
                coefficient(1,(2*kk-1))=(-1)^(kk-1);
                coefficient(1,(2*kk))=(-1)^(kk);
            end

            psiall = zeros(1,num_inversion);

            for i2=1:num_inversion
                eval(sprintf('psiall(1,%d) = (psisum_%d).*(coefficient(1,%d));',i2,i2,i2));
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

    %% Save data and fit (free)

    I_norm_final = nanmean(I_norm,1);       
    minus_ln_signal = -log(I_norm_final);

    TE = [(TE_max/TE_points):(TE_max/TE_points):TE_max].*1e3; %ms

    W = polyfit(TE,minus_ln_signal,1);
    t_plot = min(TE):(max(TE)-min(TE))/300:max(TE);
    y = minus_ln_signal;
    yfit = W(1)*TE+W(2);
    yplot = W(1)*t_plot+W(2);

    Rsq_free = 1 - sum(rmmissing((y - yfit).^2))/sum(rmmissing((y - mean(y)).^2));
    R2 = W(1)*1e3;      % in the unit of s^-1

    R2_matrix(h,1) = R2;

end 

R2_matrix_real_free = real(R2_matrix);

%% Calculate delta R2 ratio

delta_R2 = R2_matrix_real_agg - R2_matrix_real_free;
delta_R2_ratio = delta_R2./R2_matrix_real_free.*100;

end
