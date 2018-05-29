% Script to use OrderMattersSim

%       EGF      ERKpp_SOS1_FB,  ERKpp_MEK_RAF1_FB,  EGFR_tot, RAS_tot,SOS_tot,RasGAP_tot, RAF_tot,MEK_tot,ERK_tot %10-
p_init=[5e-5     1               1                   3e5       6e4     1e5     6e3         5e5     2e5     3e6    ...
... d1,     b1      u1a,    u1b     b2a     u2a b2b     u2b     k2a     k2b     b3      u3,     k3,     a2      d2, %25-    
    0.01    1e-5    0.01    100     1e-6    1   1e-7    1       1e-4    1e-5    1e-5    0.01    100    1e-7    0.01 ...
... p1      q1,     p2      p3      p4      p6      q2,     q3,     q4,     q5,     q6
    1e-7    0.01    3e-6    3e-9    6e-10   6e-10   0.01    3e-4    3e-4    100     3e-4 ...
];
p_init(1)=(1e-9)/(5e-5);
EGF_realistic=(700e-9)/(5e-5);
p_init(1)=EGF_realistic;
p_init(1)=(1000e-9)/(5e-5);

%% Set mutation effects

%mut=struct;         % The struct that contains all mutations

% Mut1
mut1=struct;
mut1.name='EGFR_T790M';
mut1.params=[1];
mut1.direct={'+'};    % Direction on mutation (increasing, decreasing, both)
allMutParam(p_init, mut1)

% Mut2
mut2=struct;
mut2.name='SOS_active_unbound';
mut2.params=[12, 19];
mut2.direct={'+-', '+'};
allMutParam(p_init, mut2)

% Mut3
mut3=struct;
mut3.name='Less_SOS-Ras_binding';
mut3.params=[15 17 19 20];
mut3.direct={'-', '-', '-', '-'};
allMutParam(p_init, mut3)

% Mut4
mut4=struct;
mut4.name='EGFR_tot';
mut4.params=[4];
mut4.direct={'specific values'};    % Not a declaration of direction, but of parameter value
mut4.values=[1e3    1e4    5e4    8e4   3e5]; mut4.values=linspace(1e3, 3e5, 10);
allMutParam(p_init, mut4)

% Mut5
mut5=struct;
mut5.name='SOS1_trunc';
mut5.params=[13 14 33];
mut5.direct={'-', '-', '0'};
allMutParam(p_init, mut5)

% Mut6
mut6=struct;
mut6.name='RasGAP';
mut6.params=[21, 22, 23];
mut6.direct={'-', '-', '+'};
allMutParam(p_init, mut6)

% Mut9
mut9=struct;
mut9.name='Ras_G12V';
mut9.params=[23];
mut9.direct={'0'};
allMutParam(p_init, mut9)

% Mut10
mut10=struct;
mut10.name='Ras_F156L';
mut10.params=[19 23 24];
mut10.direct={'+', '-', '-'};
allMutParam(p_init, mut10)

% Mut13
mut13=struct;
mut13.name='ERK';
mut13.params=[29 31];
mut13.direct={'+-', '+-'};
allMutParam(p_init, mut13)

% Mut14
mut14=struct;
mut14.name='BRAF_V600E';
mut14.params=[8 24 26 31];
mut14.direct={'specific values', '+-', '+-', '+-'};
mut14.values=[5e4];
allMutParam(p_init, mut14)

% Mut15
mut15=struct;
mut15.name='SOS_Q977R';
mut15.params=[19, 20];
mut15.direct={'-', '-'};
allMutParam(p_init, mut15)

% MutSyn
mutSyn=struct;
mutSyn.name='EGF_decrease';
mutSyn.params=[1];
mutSyn.direct={'specific values'};
mutSyn.values=linspace((1000e-9)/(5e-5), (0.1e-9)/(5e-5), 10);
allMutParam(p_init, mutSyn)


simTwoMuts(p_init, mut4, mutSyn) 


%% Single simulation
% pars =  [1              2];
% vals1=  [(1000e-9)/(5e-5)   p_init(2)];
% vals2=  [(5e-9)/(5e-5)      0];
% secstart=2; display=1; mutName=' ';
% [a b c]=OrderMattersSim(pars,vals1,vals2,secstart, display, mutName)
