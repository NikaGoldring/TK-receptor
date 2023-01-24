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

% Internal concentrations of THI in Gammarus pulex in [µmol/kg] 
DATA{3} = [0.50	1	1	1	2	2	2	3	3	3	4	4	5	5	5
0.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.0000	0.0000	NaN	NaN	NaN
0.12	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.2573	0.1895	NaN	NaN	NaN
0.25	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.3612	0.4987	NaN	NaN	NaN
0.42	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.5398	0.5323	NaN	NaN	NaN
1.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.7147	0.7789	NaN	NaN	NaN
1.42	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.8121	0.7184	NaN	NaN	NaN
2.00	NaN	NaN	NaN	NaN	NaN	NaN	0.2674	0.2467	0.2269	0.7640	0.9192	10.8444	9.7258	10.4270
2.13	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.6125	0.5869	NaN	NaN	NaN
2.25	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.5525	0.5186	NaN	NaN	NaN
2.42	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.4769	0.4182	NaN	NaN	NaN
3.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.3570	0.3204	NaN	NaN	NaN
4.00	NaN	NaN	NaN	NaN	NaN	NaN	0.2215	0.2506	0.2312	0.2622	0.3317	0.2627	0.2896	0.3041
6.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.2856	0.2446	NaN	NaN	NaN
8.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.2471	0.1902	NaN	NaN	NaN
10.00	0.0332	0.0310	0.0463	0.2118	0.2007	0.1516	NaN	NaN	NaN	0.2462	0.2213	NaN	NaN	NaN
15.00	0.0435	0.0704	0.0475	0.2185	0.2094	0.1611	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
20.00	0.0895	0.0626	0.0783	0.2307	0.1944	0.2522	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
25.00	0.0899	0.0864	0.0995	0.2247	0.1847	0.1377	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
];

% If needed, weights for individual measurements can be defined
% For this, uncommend the following line and specify your weights

% W{1} = 21 * ones(size(DATA{1})-1); % each point is a pooled sample of 21 animals

% In this data set, exposure was time-varying and reported as a series of
% concentrations over time. Here, the scenario is used as a linear forcing
% series (which has an analytical solution, and is thus much faster than
% the ODE version). Double time entries are used, which is more efficient,
% and probably more accurate.
% [µmol/L]
Cw1 = [ 0.00	1	    2	3	    4	    5
    0.00	0.0002	0.0021	0.0219	0.2206	5.0173
    0.12	0.0002	0.0021	0.0219	0.2206	5.0173
    0.25	0.0002	0.0021	0.0219	0.2206	5.0173
    0.42	0.0002	0.0021	0.0219	0.2206	5.0173
    1.00	0.0002	0.0021	0.0219	0.2206	5.0173
    1.42	0.0002	0.0021	0.0219	0.2206	5.0173
    2.00	0.0002	0.0021	0.0219	0.2206	5.0173
    2.13	0.0002	0.0021	0.0000	0.0000	0.0000
    2.25	0.0002	0.0021	0.0000	0.0000	0.0000
    2.42	0.0002	0.0021	0.0000	0.0000	0.0000
    3.00	0.0002	0.0021	0.0000	0.0000	0.0000
    4.00	0.0002	0.0021	0.0000	0.0000	0.0000
    6.00	0.0002	0.0021	0.0000	0.0000	0.0000
    8.00	0.0002	0.0021	0.0000	0.0000	0.0000
    10.00	0.0002	0.0021	0.0000	0.0000	0.0000
    15.00	0.0002	0.0021	0.0000	0.0000	0.0000
    20.00	0.0000	0.0000	0.0000	0.0000	0.0000
    25.00	0.0000	0.0000	0.0000	0.0000	0.0000];

make_scen(2,Cw1); % prepare as linear-forcing function interpolation (can use glo.use_ode = 0)  

% Create a table with nicer labels for the legends
Scenario = [1;2;3;4;5]; 
Label = {'0.05 ug/L';'0.5 ug/L';'5 ug/L';'50 ug/L (TK-main)';'1500 ug/L'};
glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = Scenario; % scenarios (concentrations or identifiers)
X0mat(2,:) = 0;      % initial values state 1 (structure internal concentrations)
X0mat(3,:) = 0;      % initial values state 2 (receptor-antagonist complex concentration)
X0mat(4,:) = 0;      % initial values state 3 (total internal concentrations)

glo.R_mod = 1 % choose kinetics for receptor model, (1) Michaelis-Menten Kinetics, or (2) second order kinetics
%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

glo.FMS = 0.01 ; % normalize for membrane protein content assuming a density of 1 (VMP/VS)

% syntax: par.name = [startvalue fit(0/1) minval maxval scale];
par.ke    = [2.596  1    1 10  1];  % elimination rate constant, d-1
par.ku    = [5.289  1    1 20  1];  % uptake rate constant, L/kg/d
par.kon   = [83.94  1   10 200 1];  % association of ligand-receptor complex
%par.koff  = [0.00001  1     0 10 1];  % dissociation of ligand-receptor complex
par.B_MAX = [25   1 0    100   1];  % maximal binding capacity, µmol/kg
par.Kd    = [0.6    1 0.001 1   0];  % equilibrium dissociation constant, nmol

switch glo.R_mod % make sure that right parameters are fitted
    case 1 % (1) Michaelis-Menten Kinetics, or 
        par.kon(2) = 0; % Do not fit kon, the association of ligand-receptor complex
    case 2 % (2) second order kinetics
        par.Kd(2)  = 0; % Do not fit Kd, equilibrium dissociation constant
end
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
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % call the plotting routine again to plot fits with CIs

% %% Calculations and plotting
% % Here, the function is called that will do the calculation and the plotting.
% % Options for the plotting can be set using opt_plot (see prelim_checks.m).
% % Options for the optimsation routine can be set using opt_optim. Options
% % for the ODE solver are part of the global glo. 
% 
% opt_optim.type = 1; % optimisation method 1) simplex, 4) parameter-space explorer
% opt_optim.fit  = 1; % fit the parameters (1), or don't (0)
% opt_optim.it   = 0; % show iterations of the simplex optimisation (1, default) or not (0)
% opt_plot.bw    = 0; % plot in black and white
% opt_plot.cn    = 0; % if set to 1, connect model line to points (only for bw=1)
% opt_plot.annot = 1; % annotations in sub-plot: text box with parameter estimates or overall legend
% glo.useode     = 1; % use the analytical solution in simplefun.m (0) or the ODE solution in derivatives (1)
% glo.stiff      = 0; % use ode45 with very strict tolerances
% 
% opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots to monitor progress of parameter-space explorer
% opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
% opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
% opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
% 
% % optimise and plot (fitted parameters in par_out)
% par_out = calc_optim(par,opt_optim); % start the optimisation
% % no plotting here; we'll immediately plot with CIs below
% 
% %% Plot results with confidence intervals
% % The following code can be used to make a standard plot (the same as for
% % the fits), but with confidence intervals. Options for confidence bounds
% % on model curves can be set using opt_conf (see prelim_checks).
% % 
% % Use opt_conf.type to tell calc_conf which sample to use: 
% % -1) Skips CIs (zero does the same, and an empty opt_conf as well).
% % 1) Bayesian MCMC sample (default); CI on predictions are 95% ranges on 
% % the model curves from the sample 
% % 2) parameter sets from a joint likelihood region using the shooting 
% % method (limited sets can be used), which will yield (asymptotically) 95% 
% % CIs on predictions
% % 3) as option 2, but using the parameter-space explorer
% 
% opt_conf.type    = 3; % make intervals from 1) slice sampler, 2) likelihood region shooting, 3) parspace explorer
% opt_conf.lim_set = 0; % use limited set of n_lim points (1) or outer hull (2, likelihood methods only) to create CIs
% opt_conf.sens    = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time
% 
% out_conf = calc_conf(par_out,opt_conf);   % calculate confidence intervals on model curves
% calc_and_plot(par_out,opt_plot,out_conf); % call the plotting routine again to plot fits with CIs
