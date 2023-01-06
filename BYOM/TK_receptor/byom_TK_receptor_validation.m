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

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% Validation is done for pulsed and constant exposures seperatly,
% determinded by the switch case here:

glo.V_mod = 2 ; % Validate for (1) pulse or (2) constant exposures

switch glo.V_mod % make sure that right parameters are fitted
    case 1 % (1) Pulse exposure
        % Internal concentrations of THI in Gammarus pulex in [µg/kg] 
        % Pulse exposure
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

        % In this data set, exposure was time-varying and reported as a series of
        % concentrations over time. Here, the scenario is used as a linear forcing
        % series (which has an analytical solution, and is thus much faster than
        % the ODE version). Double time entries are used, which is more efficient,
        % and probably more accurate.
        % Pulse exposure    
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

        % Create a table with nicer labels for the legends
        % Pulse exposure
        Scenario = [1;2]; 
        Label = {'Pulse 1'; 'Pulse 2'}; 

    case 2 % (2) Constant exposure
        % Constant exposures
        DATA{3} = [ 0.50	3	3	3	4	4	4	5	5	5	6	6	6
                    0.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
                    2.00	NaN	NaN	NaN	0.7314	0.7308	0.6980	4.4546	3.8172	3.9223	35.5530	43.7306	36.7030
                    4.00	0.259215324	0.264833362	0.235430288	0.2830	0.2491	0.2046	0.2527	0.2420	0.2577	0.4499	0.4797	0.2505
                    6.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
                    8.00	0.212699237	0.148274607	0.205478371	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN];

        % Constant exposures
        Cw1 = [0.00 	3	    4	    5	    6
               0.00	0.0186	0.2102	1.8002	20.2265
               2.00	0.0186	0.2102	1.8002	20.2265
               4.00	0.0186	0.0000	0.0000	0.0000
               4.01	0.0000	0.0000	0.0000	0.0000
               6.00	0.0000	0.0000	0.0000	0.0000
               8.00	0.0000	0.0000	0.0000	0.0000];
        % Create a table with nicer labels for the legends
        % Constant exposures
        Scenario = [3;4;5;6]; 
        Label = {'5 ug/L'; '50 ug/L'; '500 ug/L'; '5000 ug/L'}; 
end
        

make_scen(2,Cw1); % prepare as linear-forcing function interpolation (can use glo.use_ode = 0)  

glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = Scenario; % scenarios (concentrations or identifiers)
X0mat(2,:) = 0;      % initial values state 1 (structure internal concentrations)
X0mat(3,:) = 0;      % initial values state 2 (receptor-antagonist complex concentration)
X0mat(4,:) = 0;      % initial values state 3 (total internal concentrations)

glo.R_mod = 2; % choose kinetics for receptor model, (1) Michaelis-Menten Kinetics, or (2) second order kinetics
%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

switch glo.R_mod % make sure that right parameters are fitted
    case 1 % (1) Michaelis-Menten Kinetics
        % syntax: par.name = [startvalue fit(0/1) minval maxval];
        par.ke    = [5.363  1 0.01 100 1];  % elimination rate constant, d-1
        par.ku    = [10.95   1 0.01 1e6 1];  % uptake rate constant, L/kg/d
        par.kon   = [0      0 0 100 1];  % needs to be defined but is not used
        %par.koff  = [0.3116  1 0.01 100 1];  % dissociation of ligand-receptor complex
        par.B_MAX = [0.1487    1 0    100 1];  % maximal binding capacity, µmol/kg
        par.Kd    = [0.01959     1 0    100 1];  % equilibrium dissociation constant, nmol
    case 2 % (2) second order kinetics
        % syntax: par.name = [startvalue fit(0/1) minval maxval];
        par.ke    = [2.595  1 0.01 100 1];  % elimination rate constant, d-1
        par.ku    = [5.291   1 0.01 1e6 1];  % uptake rate constant, L/kg/d
        par.kon   = [83.94       1 0.01 100 1];  % association of ligand-receptor complex
        %par.koff  = [0.3116  1 0.01 100 1];  % dissociation of ligand-receptor complex
        par.B_MAX = [0.2326    1 0    100 1];  % maximal binding capacity, µmol/kg
        par.Kd    = [0     0 0    100 1] % needs to be defined but is not used
end
%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

switch glo.V_mod % time vector according to validation modus
    case 1 % (1) pulse
        glo.t = linspace(0,20,500); % time vector for the model curves in days
    case 2 % (2) constant
        glo.t = linspace(0,10,500); % time vector for the model curves in days
end

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
% Also try plottype 5 for this example! See the SIMbyom package for more
% detailed information on the simulator options.
 
opt_sim.plottype = 3; % select a plotting option
% 1) 3d plot, 
% 2) 2d plot, 
% 3) states in subplots, 
% 4) scenarios in subplots,
% 5) dx/dt versus x
% 6) 2d plot for each scenario separate

opt_sim.plot_int = 1;      % interval for plotting (how many time points to do in one go) 

sim_and_plot(par,opt_sim); % call the script which calculates and plots (simulation only)
calc_and_plot(par,opt_plot); % call the plotting routine again to plot raw data
