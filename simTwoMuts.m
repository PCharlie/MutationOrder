function [] = simTwoMuts(p, mutx, muty)
% Runs OrderMattersSim for all values of the two mutations mux and muty
% p is parameter vector with all unmutated values
% Saves all figures in:
savepath='C:/Users/Charlie/Documents/Thesis/Github/Simulation/';

% Which parameters are altered?
params={'EGF', 'ERK_{pp}-SOS1 FB', 'ERK_{pp}-MEK-RAF1 FB', 'EGFR_{tot}', 'RAS_{tot}','SOS_{tot}','RasGAP_{tot}', 'RAF_{tot}','MEK_{tot}','ERK_{tot}',...
'd_1', 'b_1', 'u_{1a}', 'u_{1b}', 'b_{2a}', 'u_{2a}', 'b_{2b}', 'u_{2b}', 'k_{2a}', 'k_{2b}', 'b_3', 'u_3', 'k_3', 'a_2', 'd_2',...
'p_1', 'q_1', 'p_2', 'p_3', 'p_4', 'p_6', 'q_2', 'q_3', 'q_4', 'q_5', 'q_6'};
pars = [mutx.params muty.params];
numpars = length(pars);
secstart = length(mutx.params)+1; % Params changed by second mutation starts at this index
% Starting values of these parameters
vals1 = p(pars);

% Gather all possible mutated parameter values for mutation X
simParamsx=allMutParam(p,mutx);
% and for mutation Y
simParamsy=allMutParam(p,muty);

% Combine, simulate, save
n=size(simParamsx,2);
m=size(simParamsy,2);

disp(['Starting to do ', num2str(n*m), ' simulations.'])

% Loop over sim. parameters for X (columns of simParamsx)
num = 1;    % Simulation number

% Save statistics of simulations in stats:
% p1 a1 a1 a1 a2 a2 a2 (mutx)
% p2 b1 b2 b3 b1 b2 b3 (muty)                                   (row p)
% ERK_mutx/ERK_SS
% ERK_muty/ERK_SS
% Added later: Behavior: 'on', 'off', 'oscill', 'bistability'   (row p+3)
stats = {};
behavior = zeros(1,6);  % on, off, osc, bistable, ND, other

for i=1:n
    % ...and Y
    disp(['Doing simulation ', num2str(num), '...'])
    for j=1:m
        % Combine two columns and simulate
        vals2 = [simParamsx(:,i)'    simParamsy(:,j)'];
        for k=1:numpars
            stats{k,num} = vals2(k);
        end

        % Create name of figure
        foldername = [savepath, mutx.name, '_', muty.name, '/'];
        warning('off', 'MATLAB:MKDIR:DirectoryExists');
        mkdir([savepath, mutx.name, '_', muty.name, '/']);
        mutName=[savepath, mutx.name, '_', muty.name, '/', mutx.name, '_', muty.name, '_', num2str(num)];
        [a b c] = OrderMattersSim(pars, vals1, vals2, secstart, 'none', mutName); % 'none', 0, 1 (create stats, jpg, pop-up figs)
        stats{ (numpars+1), num} = a;
        stats{ (numpars+2), num} = b;
        stats{ (numpars+3), num} = c;
        behavior = addToBehavior(behavior, c);
        num=num+1;
    end
end

% Get overall stats and save them
% behavior = array2table([mean(ismember('on',stats(numpars+3,:))), mean(ismember('off',stats(numpars+3,:))), mean(ismember('osc',stats(numpars+3,:))), ...
%                         mean(ismember('bistable',stats(numpars+3,:)))],  'VariableNames', {'on', 'off', 'osc', 'bistable'});
behavior = array2table(behavior, 'VariableNames', {'on', 'off', 'osc', 'order_matters', 'ND', 'other'});
save([savepath, mutx.name, '_', muty.name, '/_Behavior_stats_', mutx.name, '_', muty.name], 'behavior');

% Save stats
pars_cell = {params{pars}};
varnames = strcat('sim' , strtrim(cellstr(num2str((1:n*m)'))'));
stats = array2table(stats, 'RowNames', {pars_cell{:}, 'ERK_{pp} act mut_x', 'ERK_{pp} act mut_y', 'behavior'}, 'VariableNames', varnames);
writetable(stats, [savepath, mutx.name, '_', muty.name, '/_Statstable_', mutx.name, '_', muty.name], 'WriteRowNames', true);
end

function behavior = addToBehavior(behavior, c)
if strcmp('on', c)
    behavior(1)=behavior(1)+1;
elseif strcmp('off', c)
    behavior(2)=behavior(2)+1;
elseif strcmp('osc', c)
    behavior(3)=behavior(3)+1;
elseif strcmp('order matters', c)
    behavior(4)=behavior(4)+1;
elseif strcmp('ND', c)
    behavior(5)=behavior(5)+1;
elseif strcmp('other', c)
    behavior(6)=behavior(6)+1;
else
    error('Behavior undefined. Check function OrderMattersSim.');
end
end