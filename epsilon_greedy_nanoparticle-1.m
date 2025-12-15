%---------------------------------------------------------------------------------------------------------------------------------------------------------------------
% This code uses the epsilon-greedy reinforcement learning strategy to
% optimize the design of nanometers.
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------

clear;
close all;

%%
% specify the input parameters as follows: 
% 1) epsilon: random threshold;
% 2) dp_input: a set of actions a_i (mapping to different diameter of NP in
%    nanometers);
% 3) Q_estmates[i]: initial value-estimates of the action a_i;
% 4) k[i]: the number of times the action a_i has been taken;
% 5) numTests: total number of tests.
% 6) resMatrix: 
%   (a) dimension: numTests * 6, 
%   (b) each row: [dnp, R2_free, Rsq_free, R2_agg, Rsq_agg, delta_R2_ratio]

folder = uigetdir;   
date = '/20241108_';
mkdir(folder, date);    % Change path 
path = strcat(folder, date);

epsilon    = 0.25;
dp_input   = (1/12:1/12:1) .* 60;    % candidate NP diameters (nm)
Q_estmates = zeros(size(dp_input));
k          = ones(size(dp_input));   % visit counts
numTests   = 1;                      % number of epsilon-greedy iterations

iniMatrix  = zeros(max(size(dp_input)), 6);
resMatrix  = zeros(numTests, 6);

%%
% To initialize Q_estmates[i], call each action once from the set of dp_input;
for i = 1:length(dp_input)
    [delta_R2_ratio, Rsq_agg, Rsq_free, R2_matrix_real_agg, R2_matrix_real_free] = ...
        Glutamate_DrWEI_Example-1(dp_input(i));

    % initial value estimate = observed delta_R2_ratio
    Q_estmates(i) = delta_R2_ratio;

    % record initial results
    iniMatrix(i, 1) = dp_input(i);            % NP diameter
    iniMatrix(i, 2) = R2_matrix_real_free;    % R2_free
    iniMatrix(i, 3) = Rsq_free;               % R^2 fit free
    iniMatrix(i, 4) = R2_matrix_real_agg;     % R2_agg
    iniMatrix(i, 5) = Rsq_agg;                % R^2 fit agg
    iniMatrix(i, 6) = delta_R2_ratio;         % delta R2 ratio
end

initial_results = real(iniMatrix);
save(strcat(path, '/initial_results'), 'initial_results');

%% 
% exploitation and exploration algorithm by epsilon-greedy.
for idx = 1:numTests
    % random number p in (0,1) for epsilon-greedy
    p = rand(1);

    % current best estimated value
    maxVal      = max(Q_estmates);
    multipleIdx = find(Q_estmates == maxVal);

    if length(multipleIdx) > 1
        randChooseIdx = randi(length(multipleIdx), 1);
        OptIndex      = multipleIdx(randChooseIdx);
    else 
        [~, OptIndex] = max(Q_estmates);
    end

    if p >= epsilon
        updateIndex = OptIndex;
    else
        tempIdxSet = [];
        for otherIdx = 1:length(dp_input)
            if Q_estmates(otherIdx) == maxVal
                continue;
            else
                tempIdxSet = [tempIdxSet; otherIdx];
            end
        end

        if isempty(tempIdxSet)
            updateIndex = randi(length(dp_input), 1);
        else
            randIndex   = randi(length(tempIdxSet), 1);
            updateIndex = tempIdxSet[randIndex];
        end
    end

    [delta_R2_ratio, Rsq_agg, Rsq_free, R2_matrix_real_agg, R2_matrix_real_free] = ...
        SW_MaCaReNa2_MSME_agg_free_function_v1b(dp_input(updateIndex));

    rewardVal = delta_R2_ratio;

    resMatrix(idx, 1) = dp_input(updateIndex);
    resMatrix(idx, 2) = R2_matrix_real_free;
    resMatrix(idx, 3) = Rsq_free;
    resMatrix(idx, 4) = R2_matrix_real_agg;
    resMatrix(idx, 5) = Rsq_agg;
    resMatrix(idx, 6) = delta_R2_ratio;

    Q_estmates(updateIndex) = ...
        (Q_estmates(updateIndex) .* k(updateIndex) + rewardVal) ./ (1 + k(updateIndex));
    k(updateIndex) = k(updateIndex) + 1;
end

[optVal, optI] = max(Q_estmates);

simulation_results = real(resMatrix);
save(strcat(path, '/simulation_results'), 'simulation_results');

figure;

plot(dp_input, Q_estmates, '--rs', 'LineWidth', 2, ...
                       'MarkerEdgeColor', 'k', ...
                       'MarkerFaceColor', 'g', ...
                       'MarkerSize', 10);

xlabel('NP diameter (nm)');
ylabel('deltaR2 ratio (a.u.)'); 
caption = sprintf('simulation result of %d tests and epsilon = %0.2f', numTests, epsilon);
title(caption, 'FontSize', 14);
saveas(gcf, strcat(path, '/deltaR2_ratio_result.fig'));
close;
