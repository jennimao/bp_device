addpath('toolbox/signal')


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
bioz_plot = myTable.BioZ_Transformed_ohm_;
ecg_plot = myTable.ECG_mV_;


% Convert the time intervals to seconds
time_ecg = seconds(time_intervals_ecg); % assuming time intervals are in minutes
time_bioz = seconds(time_intervals_bioz);

y_new = interp1(time_ecg, ecg_plot, time_bioz);

y_new

%plot(time_bioz, bioz)
%xlim([seconds(5), seconds(15)]); % Change the values as needed
%ylim([-5, 5])


% Plot both sets of y values on the same graph
plot(time_bioz, bioz_plot, 'b', time_bioz, y_new*0.9+0.2, 'r')
legend('BioZ', 'ECG')
xlabel('x')
ylabel('y')
title('Plot of BioZ and ECG Data')
xlim([seconds(8), seconds(11)]); % Change the values as needed
ylim([-2, 5])

%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleaning BioZ Data 
%%%%%%%%%%%%%%%%%%%%%%%%%

bioz = myTable.BioZ_ohm_;
nan_indices = isnan(bioz);
% Remove the NaN values
bioz = bioz(~nan_indices);
nan_indices


mean_val = mean(bioz); % Mean value
std_val = std(bioz); % Standard deviation
max_val = max(bioz); % Maximum value

bioz

mean_val, std_val, max_val

% Extract features
pp_amp = peak2peak(bioz); % Peak-to-peak amplitude
rms_val = rms(bioz); % RMS value



%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleaning ECG Data 
%%%%%%%%%%%%%%%%%%%%%%%%%

ecg = myTable.ECG_mV_;

%{

% bandpass filter? 
[b,a] = butter(2, [0.05 30]/(fs/2), 'bandpass');
ecg_filtered = filtfilt(b, a, ecg_data);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating PAT 
%%%%%%%%%%%%%%%%%%%%%%%%%

% Resample both signals -- COMMENT OUT FOR NOW 
%{
fs_resample = 250; % desired sampling rate
ecg_resampled = resample(ecg, fs_resample, fs);
bioz_resampled = resample(bioz, fs_resample, fs_bioimp);
%}


% find delay between R peak and bioimpedance upstroke 
[acor,lag] = xcorr(ecg_resampled, bioz_resampled);
[~,I] = max(abs(acor));
delay = lag(I);
bioz_shifted = circshift(bioz_resampled, [delay 0]);

% calculate pulse transit time (PTT) and pulse arrival time (PAT) 
[~,r_peak_locs] = findpeaks(ecg_resampled, 'MinPeakHeight', 0.6, 'MinPeakDistance', fs_resample*0.6);
[~,bioz_upstroke_locs] = findpeaks(bioz_shifted, 'MinPeakHeight', 100);
ptt = (bioz_upstroke_locs(1:length(r_peak_locs))-r_peak_locs)/fs_resample;
pat = ptt - 1/fs_resample*r_peak_locs';


%}


