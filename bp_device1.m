myTable = readtable('/Users/jennymao/Downloads/bioz_ecg_4_19.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic Data Table Stats  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%summary(myTable)
%myTable.Properties.VariableNames
% Print out some sample values of the x_Time2_s_ column
%disp(myTable.x_Time2_s_(1:10));

time_intervals_ecg = myTable.x_Time2_s_;
time_intervals_bioz = myTable.x_Time_s_;
bioz = myTable.BioZ_Transformed_ohm_;
ecg = myTable.ECG_mV_;


% Convert the time intervals to seconds
time_ecg = seconds(time_intervals_ecg); % assuming time intervals are in minutes
time_bioz = seconds(time_intervals_bioz);

y_new = interp1(time_ecg, myTable.('ECG_mV_'), time_bioz);

y_new

plot(time_bioz, bioz)
xlim([seconds(5), seconds(15)]); % Change the values as needed
ylim([-5, 5])


% Plot both sets of y values on the same graph
plot(time_bioz, bioz, 'b', time_bioz, y_new*0.9+0.2, 'r')
legend('BioZ', 'ECG')
xlabel('x')
ylabel('y')
title('Plot of BioZ and ECG Data')
xlim([seconds(8), seconds(11)]); % Change the values as needed
ylim([-2, 5])

%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleaning BioZ Data 
%%%%%%%%%%%%%%%%%%%%%%%%%