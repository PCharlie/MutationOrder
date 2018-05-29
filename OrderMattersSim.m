function [ERK_mutx, ERK_muty, behav]=OrderMattersSim(pars,vals1,vals2,secstart, display, mutName)
% Outputs a simulation (t vs ERKpp) wtih two mutations in both orders
% pars are the indices of the parameters that is mutated; vector of length n
% vals1 is a 1xn vector with the starting values for all parameters that are mutated
% vals2 contains the mutated values
% secstart is the index of what is considered the second mutation
% pars = [pars1, pars2, pars3, pars4] and secstart=2, then the second mutation consists of three par changes
% display can be True (show Figure), False (save Figure) or 'None' (don't create plots)

% Default parameters
%       EGF      ERKpp_SOS1_FB,  ERKpp_MEK_RAF1_FB,  EGFR_tot, RAS_tot,SOS_tot,RasGAP_tot, RAF_tot,MEK_tot,ERK_tot
p_init=[5e-5     1               1                   3e5       6e4     1e5     6e3         5e5     2e5     3e6    ...
... d1,     b1      u1a,    u1b     b2a     u2a b2b     u2b     k2a     k2b     b3      u3,     k3,     a2      d2,     
    0.01    1e-5    0.01    100     1e-6    1   1e-7    1       1e-4    1e-5    1e-5    0.01    100    1e-7    0.01 ...
... p1      q1,     p2      p3      p4      p6      q2,     q3,     q4,     q5,     q6
    1e-7    0.01    3e-6    3e-9    6e-10   6e-10   0.01    3e-4    3e-4    100     3e-4 ...
];
p_init(1)=(1e-9)/(5e-5);
EGF_realistic=(700e-9)/(5e-5);
p_init(1)=EGF_realistic;
p_init(1)=(1000e-9)/(5e-5);

params={'EGF', 'ERK_{pp}-SOS1 FB', 'ERK_{pp}-MEK-RAF1 FB', 'EGFR_{tot}', 'RAS_{tot}','SOS_{tot}','RasGAP_{tot}', 'RAF_{tot}','MEK_{tot}','ERK_{tot}',...
'd_1', 'b_1', 'u_{1a}', 'u_{1b}', 'b_{2a}', 'u_{2a}', 'b_{2b}', 'u_{2b}', 'k_{2a}', 'k_{2b}', 'b_3', 'u_3', 'k_3', 'a_2', 'd_2',...
'p_1', 'q_1', 'p_2', 'p_3', 'p_4', 'p_6', 'q_2', 'q_3', 'q_4', 'q_5', 'q_6'};

species = {'EGFR_a', 'Sos_p', 'Sos_{pp}', 'Sos_{ppp}', 'Sos_{pppp}', 'EGFR_i-Sos_u', 'EGFR_a-Sos_u', 'Ras^{GTP}', 'EGFR_i-Sos_u-Ras^{GDP}', 'EGFR_i-Sos_u-Ras^{GTP}',...
'EGFR_a-Sos_u-Ras^{GDP}', 'EGFR_a-Sos_u-Ras^{GDP}', 'EGFR_a-Sos_u-Ras^{GTP}', 'RasGAP-Ras^{GDP}', 'RasGAP-Ras^{GTP}', 'RAF_a', 'RAF_p', 'MEK_(u,p}', 'MEK_{p,u}', 'MEK_{p,p}', 'MEK_{pp,u}',...
'MEK_{pp,p}', 'ERK_p', 'ERK_{pp}'};
c_init=zeros(1,23);

% Matrices that contain per column: param index, starting value, mutated value
mut1 = [pars(1:(secstart-1)); vals1(1:(secstart-1)); vals2(1:(secstart-1))]; % First mutation
mut2 = [pars(secstart:end); vals1(secstart:end); vals2(secstart:end)]; % Second mutation

%% Simulate
% SS with starting values
n = length(pars);
p0 = p_init;
for i=1:n       % Set all starting values
    index = pars(i);
    p0(index)=vals1(i);
end
[t1 x1] = modelsim(p0,c_init);
x0 = x1(end,:);

% First order of mutations: 1...(secstart-1)
p=p0;
for i=1:(secstart-1) % Mut 1
    index = pars(i);
    p(index)=vals2(i);
end
[t2a x2a] = modelsim(p, x0);
c_temp = x2a(end,:);

% Calculate ERK fold increase; given no oscillations at start
% Check for oscillations after mutx
x_sechalf     = x2a( int32((size(x2a,1))*3/4):end, 23);
% If ERK is within a 10% range of the final value, then there is no osc.
osc     = ~all(0.90*x2a(end,23) < x_sechalf) && ~all(1.10*x2a(end,23) > x_sechalf);
if osc
    ERK_mutx = NaN;
else
    ERK_mutx = x2a(end,23)/x0(23);
end

for i=secstart:n    % Mut 2
    index = pars(i);
    p(index)=vals2(i);
end
[t3a x3a] = modelsim(p, c_temp);

% Second order of mutations
p=p0;
for i=secstart:n    % Mut 2
    index = pars(i);
    p(index)=vals2(i);
end
[t2b x2b] = modelsim(p, x0);
c_temp=x2b(end,:);

% Calculate ERK fold increase
x_sechalf     = x2b( int32((size(x2b,1))*3/4):end, 23);
osc     = ~all(0.90*x2b(end,23) < x_sechalf) && ~all(1.10*x2b(end,23) > x_sechalf);
if osc
    ERK_muty = NaN;
else
    ERK_muty = x2b(end,23)/x0(23);
end

for i=1:(secstart-1) % Mut 1
    index = pars(i);
    p(index)=vals2(i);
end
[t3b x3b] = modelsim(p, c_temp);

behav = behavior(t3a, t3b, x3a, x3b, p_init);

%%
if display~='none'  % Plot

% Don't display if plot option if false
if ~display
    figure('Visible', 'off');
end

% Plot SS
plot(t1, x1(:,23), 'color', 'black')
t_high = t1(end);
hold on;

% Plot first order of mutations
t2a = t2a+t_high;
h1 = plot(t2a, x2a(:,23), 'color', 'blue');
t3a = t3a+2*t_high;
plot(t3a, x3a(:,23), 'color', 'blue')

% Plot second order of mutations
t2b = t2b+t_high;
h2 = plot(t2b, x2b(:,23), 'color', 'red');
t3b = t3b+2*t_high;
plot(t3b, x3b(:,23), 'color', 'red')
xlabel('Time'); ylabel('ERK_{pp}');
line([t_high t_high], [0 p0(10)], 'color', 'black', 'linestyle', '--');
line([2*t_high 2*t_high], [0 p0(10)], 'color', 'black', 'linestyle', '--');
% Legend
mut1_text = []; mut2_text = [];
for i=1:(secstart-1)
    mut1_text = [mut1_text, params{pars(i)}, ': ', num2str(vals1(i)), ' to ', num2str(vals2(i))];
    if i<(secstart-1)
        mut1_text = [mut1_text ' | ']; end
end
for i=secstart:n
    mut2_text = [mut2_text, params{pars(i)}, ': ', num2str(vals1(i)), ' to ', num2str(vals2(i))];
    if i<n
        mut2_text = [mut2_text ' | ']; end
end
legend([h1, h2],{[mut1_text 10 mut2_text], [mut2_text 10 mut1_text]});

% Save if display is false
if ~display
    saveas(gcf, mutName, 'jpg');
end
end
end

function [t, x]=modelsim(p,c)

model = MAPK_short_fbs_2in1;
tspan = linspace(0,5e4,2001);
[t, x] = ode23s(model{2}, tspan, c, [], p(1), p(2), p(3), ...
    p(4), p(5), p(6), p(7), p(8), p(9), p(10), p(11), p(12), p(13),  ...
    p(14), p(15), p(16), p(17), p(18), p(19), p(20), p(21), p(22), ... 
    p(23), p(24), p(25), p(26), p(27), p(28), p(29), p(30), p(31), ...
    p(32), p(33), p(34), p(35), p(36));
end

function [behav]=behavior(t3a, t3b, x3a, x3b, p_init)

% Find behavior after two mutations
% Check osc
x_sechalf = x3a( int32((size(x3a,1))*3/4):end, 23);
%osc1      = ~( issorted(x_sechalf) || issorted( fliplr(x_sechalf) ) );
osc1     = ~all(0.90*x3a(end,23) < x_sechalf) || ~all(1.10*x3a(end,23) > x_sechalf);
x_sechalf = x3b( int32((size(x3b,1))*3/4):end, 23);
%osc2      = ~( issorted(x_sechalf) || issorted( fliplr(x_sechalf) ) );
osc2     = ~all(0.90*x3b(end,23) < x_sechalf) || ~all(1.10*x3b(end,23) > x_sechalf);
if osc1 && osc2
    behav = 'osc';
elseif osc1 || osc2
    behav = 'ND';   % Mutations are not independent
else
    % Check bistability: The difference between the end states is more than 1% ERK activity
    if (abs(x3a(end,23)-x3b(end,23)) > p_init(10)*0.01)
        behav = 'order matters';
    else
        % Check monostable 'on': More than 10% of ERK is active
        if (x3a(end,23) > p_init(10)*0.25) || (x3b(end,23) > p_init(10)*0.10)
            behav = 'on';
        elseif (x3a(end,23) < p_init(10)*0.25) || (x3b(end,23) < p_init(10)*0.10)
            behav = 'off';
        else
            behav = 'other';    % Should not happen
        end
    end
end
end