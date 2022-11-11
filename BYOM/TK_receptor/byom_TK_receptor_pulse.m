%% BYOM, byom_calanus_2016_onecomp.m 
%
% * Author: Tjalling Jager 
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* An organism is exposed to a chemical in its surrounding
% medium. The animal accumulates the chemical according to standard
% one-compartment first-order kinetics.  
%
% *This script:* One-compartment TK of C2-naphthalene in _Calanus finmarchicus_.
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 
%
%
% September 2022 
% Modifications by AMD to include (irrevirsible) receptor binding of an 
% antagonist. Here with the example of the antagonist thiacloprid (THI) 
% binding to the nicotinic-acetylcholine receptor (nAChR) as observed 
% in Gammarus pulex. in Raths et al. 202X


%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
diary off           % turn of the diary function (if it is accidentaly on)
set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
% set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(0) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

% Internal concentrations of THI in Gammarus pulex in [unit] 

DATA{3} = [ 0.5	 1	     1	    1	    2	    2	    2
            0.0	 0.0001	 0.0000	0.0000	0.0000	0.0000	0.0000
            1.0	 0.2148	 0.1780	0.1789	0.7390	0.7474	0.7777
            2.0	 0.3031	 0.2461	0.2827	0.7887	0.6594	0.8110
            3.0	 0.2037	 0.2047	0.2073	0.2643	0.3133	0.2680
            5.0	 0.1582	 0.1726	0.2434	0.2153	0.2072	0.2006
            6.0	 0.2465	 0.2685	0.3154	0.7725	0.7111	0.7260
            7.0	 0.2697	 0.2941	0.2847	0.7511	0.7651	0.7614
            8.0	 0.1662	 0.1879	0.2158	0.2253	0.2285	0.2332
           10.1	 0.2327	 0.2098	0.2296	0.2172	0.2323	0.1935
           11.1	 0.2859	 0.3029	0.3853	0.7698	0.7397	0.7423
           12.1	 0.3125	 0.3555	0.3843	0.7809	0.7601	0.8217
           13.1	 0.2049	 0.2347	0.1889	0.2339	0.2339	0.2465
           15.1	 0.1691	0.1914	0.2380	0.2100	0.1716	0.2111];

% If needed, weights for individual measurements can be defined
% For this, uncommend the following line and specify your weights

% W{1} = 21 * ones(size(DATA{1})-1); % each point is a pooled sample of 21 animals

% In this data set, exposure was time-varying and reported as a series of
% concentrations over time. Here, the scenario is used as a linear forcing
% series (which has an analytical solution, and is thus much faster than
% the ODE version). Double time entries are used, which is more efficient,
% and probably more accurate.
Cw1 = [ 0	1	    2       
        0.0	0.0227	0.2193  % µmol/L
        1.0	0.0227	0.2193
        2.0	0	    0
        3.0	0	    0
        5.0	0.0227	0.2193
        6.0	0.0227	0.2193
        7.0	0	    0
        8.0	0	    0
       10.1	0.0227	0.2193
       11.1	0.0227	0.2193
       12.1	0	    0
       13.1	0	    0
       15.1	0	    0 ];

make_scen(2,Cw1); % prepare as linear-forcing function interpolation (can use glo.use_ode = 0)  

% Create a table with nicer labels for the legends
Scenario = [1;2]; 
Label = {'Pulse 1'; 'Pulse 2'}; 
glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = Scenario; % scenarios (concentrations or identifiers)
X0mat(2,:) = 0;      % initial values state 1 (structure internal concentrations)
X0mat(3,:) = 0;      % initial values state 2 (receptor-antagonist complex concentration)
X0mat(4,:) = 0;      % initial values state 3 (total internal concentrations)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% syntax: par.name = [startvalue fit(0/1) minval maxval];
par.ke    = [3.121  1 0.01 100 1];  % elimination rate constant, d-1
par.ku    = [7.547   1 0.01 1e6 1];  % uptake rate constant, L/kg/d
par.kon   = [1       1 0.01 100 1];  % association of ligand-receptor complex
%par.koff  = [0.3116  1 0.01 100 1];  % dissociation of ligand-receptor complex
par.B_MAX = [0.29    1 0    100 1];  % maximal binding capacity, µmol/kg
par.Kd    = [0.6     1 0    100 1];  % equilibrium dissociation constant, nmol

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = ['internal concentration (',char(181),'mol/kg)'];
glo.ylab{2} = ['receptor-antagonist complex concentration (',char(181),'mol/kg)'];
glo.ylab{3} = ['total internal concentration (',char(181),'mol/kg)'];
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'M']; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 

glo.R_mod = 1; % choose kinetics for receptor model, (1) Michaelis-Menten Kinetics, or (2) second order kinetics

opt_optim.fit = 1; % fit the parameters (1), or don't (0)
opt_optim.it  = 0; % show iterations of the simplex optimisation (1, default) or not (0)
opt_plot.bw   = 0; % plot in black and white
opt_plot.cn   = 0; % if set to 1, connect model line to points (only for bw=1)
glo.useode    = 1; % use the analytical solution in simplefun.m (0) or the ODE solution in derivatives (1)

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

% % % Profiling the likelihood
% % By profiling you make robust confidence intervals for one or more of your
% % parameters. Use the names of the parameters as they occurs in your
% % parameter structure _par_ above. This can be a single string (e.g.,
% % 'kd'), a cell array of strings (e.g., {'kd','ke'}), or 'all' to profile
% % all fitted parameters. 
% % 
% % Options for profiling can be set using opt_prof (see prelim_checks.m).
% 
% opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
% opt_prof.subopt   = 10; % number of sub-optimisations to perform to increase robustness
% 
% % UNCOMMENT LINE(S) TO CALCULATE
% par_better = calc_proflik(par_out,{'all'},opt_prof,opt_optim);  % calculate a profile
% if ~isempty(par_better)                 % if the profiling found a better optimum ...
%     calc_and_plot(par_better,opt_plot); % calculate model lines and plot them
% end
