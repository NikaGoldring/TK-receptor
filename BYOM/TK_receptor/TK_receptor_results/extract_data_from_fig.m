clear all;
close all;
clc;

%% Specify file for data extraction
open('20230130_calibration_fit.fig') ; % choose corresponding file name of MATLAB figure
h = findobj(gcf, 'Type', 'line');

% NOTE: As there were no input data for the metabolite provided at 07
% degrees, the locations of the data for the plotted lines within the
% matlab figure are different. Make sure to use the correct locations by
% identifying them first (i.e., open xdata and ydata, first lines in
% each section)

%% Extract xdata
xdata= get(h,'XData'); % open variables stored in XData for identification
xdata1= get(h(15,1),'XData'); % X data for all lines

%% Extract ydata
ydata= get(h,'YData'); % open variables stored in YData for identification
% % Total internal conentrations (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata1= get(h(25,1),'YData'); % line 1500 µg
ydata2= get(h(15,1),'YData'); % upper 1500 µg
ydata3= get(h(16,1),'YData'); % lower 1500 µg

ydata4= get(h(26,1),'YData'); % line 50 µg
ydata5= get(h(17,1),'YData'); % upper 50 µg
ydata6= get(h(18,1),'YData'); % lower 50 µg

ydata7= get(h(27,1),'YData'); % line 5 µg
ydata8= get(h(19,1),'YData'); % upper 5 µg
ydata9= get(h(20,1),'YData'); % lower 5 µg

ydata10= get(h(28,1),'YData'); % line 0.5 µg
ydata11= get(h(21,1),'YData'); % upper 0.5 µg
ydata12= get(h(22,1),'YData'); % lower 0.5 µg

ydata13= get(h(29,1),'YData'); % line 0.05 µg
ydata14= get(h(23,1),'YData'); % upper 0.05 µg
ydata15= get(h(24,1),'YData'); % lower 0.05 µg

% % Internal conentrations in receptor-complex/membrane (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata16= get(h(40,1),'YData'); % line 1500 µg
ydata17= get(h(30,1),'YData'); % upper 1500 µg
ydata18= get(h(31,1),'YData'); % lower 1500 µg

ydata19= get(h(41,1),'YData'); % line 50 µg
ydata20= get(h(32,1),'YData'); % upper 50 µg
ydata21= get(h(33,1),'YData'); % lower 50 µg

ydata22= get(h(42,1),'YData'); % line 5 µg
ydata23= get(h(34,1),'YData'); % upper 5 µg
ydata24= get(h(35,1),'YData'); % lower 5 µg

ydata25= get(h(43,1),'YData'); % line 0.5 µg
ydata26= get(h(36,1),'YData'); % upper 0.5 µg
ydata27= get(h(37,1),'YData'); % lower 0.5 µg

ydata28= get(h(44,1),'YData'); % line 0.05 µg
ydata29= get(h(38,1),'YData'); % upper 0.05 µg
ydata30= get(h(39,1),'YData'); % lower 0.05 µg

% % Internal conentrations supertenant+debris (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata31= get(h(55,1),'YData'); % line 1500 µg
ydata32= get(h(45,1),'YData'); % upper 1500 µg
ydata33= get(h(46,1),'YData'); % lower 1500 µg

ydata34= get(h(56,1),'YData'); % line 50 µg
ydata35= get(h(47,1),'YData'); % upper 50 µg
ydata36= get(h(48,1),'YData'); % lower 50 µg

ydata37= get(h(57,1),'YData'); % line 5 µg
ydata38= get(h(49,1),'YData'); % upper 5 µg
ydata39= get(h(50,1),'YData'); % lower 5 µg

ydata40= get(h(58,1),'YData'); % line 0.5 µg
ydata41= get(h(51,1),'YData'); % upper 0.5 µg
ydata42= get(h(52,1),'YData'); % lower 0.5 µg

ydata43= get(h(59,1),'YData'); % line 0.05 µg
ydata44= get(h(53,1),'YData'); % upper 0.05 µg
ydata45= get(h(54,1),'YData'); % lower 0.05 µg


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
fig01(:,20) = ydata19 ;
fig01(:,21) = ydata20 ;
fig01(:,22) = ydata21 ;
fig01(:,23) = ydata22 ; 
fig01(:,24) = ydata23 ;
fig01(:,25) = ydata24 ;
fig01(:,26) = ydata25 ; 
fig01(:,27) = ydata26 ;
fig01(:,28) = ydata27 ;
fig01(:,29) = ydata28 ;
fig01(:,30) = ydata29 ; 
fig01(:,31) = ydata30 ;
fig01(:,32) = ydata31 ;
fig01(:,33) = ydata32 ;
fig01(:,34) = ydata33 ; 
fig01(:,35) = ydata34 ;
fig01(:,36) = ydata35 ;
fig01(:,37) = ydata36 ;
fig01(:,38) = ydata37 ; 
fig01(:,39) = ydata38 ;
fig01(:,40) = ydata39 ;
fig01(:,41) = ydata40 ;
fig01(:,42) = ydata41 ;
fig01(:,43) = ydata42 ;
fig01(:,44) = ydata43 ; 
fig01(:,45) = ydata44 ;
fig01(:,46) = ydata45 ;

dlmwrite('20230131_calibration_fit_2nd.txt', fig01, ','); % write dataframe in txt
