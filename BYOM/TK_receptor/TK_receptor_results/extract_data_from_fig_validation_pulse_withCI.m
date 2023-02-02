clear all;
close all;
clc;

%% Specify file for data extraction
open('20230130_validation_pulse_2nd.fig') ; % choose corresponding file name of MATLAB figure
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
ydata1= get(h(11,1),'YData'); % line pulse 1
ydata2= get(h(7,1),'YData'); % upper pulse 1
ydata3= get(h(8,1),'YData'); % lower pulse 1

ydata4= get(h(12,1),'YData'); % line pulse 2
ydata5= get(h(9,1),'YData'); % upper pulse 2
ydata6= get(h(10,1),'YData'); % lower pulse 2

% % Internal conentrations in receptor-complex/membrane (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata7= get(h(17,1),'YData'); % line pulse 1
ydata8= get(h(13,1),'YData'); % upper pulse 1
ydata9= get(h(14,1),'YData'); % lower pulse 1

ydata10= get(h(18,1),'YData'); % line pulse 2
ydata11= get(h(15,1),'YData'); % upper pulse 2
ydata12= get(h(16,1),'YData'); % lower pulse 2

% % Internal conentrations supertenant+debris (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata13= get(h(23,1),'YData'); % line pulse 1
ydata14= get(h(19,1),'YData'); % upper pulse 1
ydata15= get(h(20,1),'YData'); % lower pulse 1

ydata16= get(h(24,1),'YData'); % line pulse 2
ydata17= get(h(21,1),'YData'); % upper pulse 2
ydata18= get(h(22,1),'YData'); % lower pulse 2

%% Export data to txt
fig01 = []; %create empty table
fig01(:,1) = xdata1 ; % xdata for all lines
fig01(:,2) = ydata1 ; 
fig01(:,3) = ydata2 ;
fig01(:,4) = ydata3 ;
fig01(:,5) = ydata4 ; 
fig01(:,6) = ydata5 ;
fig01(:,7) = ydata6 ;
fig01(:,8) = ydata7 ; 
fig01(:,9) = ydata8 ;
fig01(:,10) = ydata9 ;
fig01(:,11) = ydata10 ; 
fig01(:,12) = ydata11 ;
fig01(:,13) = ydata12 ;
fig01(:,14) = ydata13 ;
fig01(:,15) = ydata14 ; 
fig01(:,16) = ydata15 ;
fig01(:,17) = ydata16 ;
fig01(:,18) = ydata17 ;
fig01(:,19) = ydata18 ; 

dlmwrite('20230131_validation_pulse_2nd.txt', fig01, ','); % write dataframe in txt
