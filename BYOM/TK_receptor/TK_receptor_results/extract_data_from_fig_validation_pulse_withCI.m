clear all;
close all;
clc;

%% Specify file for data extraction
open('20230106_validation_pulse_2nd.fig') ; % choose corresponding file name of MATLAB figure
h = findobj(gcf, 'Type', 'line');

% NOTE: As there were no input data for the metabolite provided at 07
% degrees, the locations of the data for the plotted lines within the
% matlab figure are different. Make sure to use the correct locations by
% identifying them first (i.e., open xdata and ydata, first lines in
% each section)

%% Extract xdata
xdata= get(h,'XData'); % open variables stored in XData for identification
xdata1= get(h(7,1),'XData'); % X data for all lines

%% Extract ydata
ydata= get(h,'YData'); % open variables stored in YData for identification
% % Total internal conentrations (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata1= get(h(8,1),'YData'); % line pulse1 µg
ydata2= get(h(7,1),'YData'); % line pulse2 µg

% % Internal conentrations in receptor-complex/membrane (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata3= get(h(10,1),'YData'); % line pulse 1 µg
ydata4= get(h(9,1),'YData'); % line pulse 2 µg

% % Internal conentrations supertenant+debris (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata5= get(h(12,1),'YData'); % line pulse 1 µg
ydata6= get(h(11,1),'YData'); % line pulse 2 µg


%% Export data to txt
fig01 = []; %create empty table
fig01(:,1) = xdata1 ; % xdata for all lines
fig01(:,2) = ydata1 ; 
fig01(:,3) = ydata2 ;
fig01(:,4) = ydata3 ;
fig01(:,5) = ydata4 ; 
fig01(:,6) = ydata5 ;
fig01(:,7) = ydata6 ;

dlmwrite('20230113_validation_pulse_2nd.txt', fig01, ','); % write dataframe in txt
