
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BioZ + ECG Overlay Plot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert the time intervals to seconds
time_ecg = seconds(time_intervals_ecg); % assuming time intervals are in minutes
time_bioz = seconds(time_intervals_bioz);
%{
y_new = interp1(time_ecg, ecg_plot, time_bioz);

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
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleaning BioZ Data and Feature Extraction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bioz = myTable.BioZ_ohm_;
nan_indices = isnan(bioz);
bioz = bioz(~nan_indices); % Remove the NaN values


mean_val = mean(bioz); % Mean value
std_val = std(bioz); % Standard deviation
max_val = max(bioz); % Maximum value
%mean_val, std_val, max_val

% Extract features
pp_amp = peak2peak(bioz); % Peak-to-peak amplitude
rms_val = rms(bioz); % RMS value
%pp_amp, rms_val



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleaning ECG Data and Feature Extraction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ecg = myTable.ECG_mV_;
nan_indices = isnan(ecg);
ecg = ecg(~nan_indices); % Remove the NaN values

mean_val_ecg = mean(ecg); % Mean value
std_val_ecg = std(ecg); % Standard deviation
max_val_ecg = max(ecg); % Maximum value

%mean_val_ecg, std_val_ecg, max_val_ecg

% bandpass filter? 
%[b,a] = butter(2, [0.05 30]/(fs/2), 'bandpass');
%ecg_filtered = filtfilt(b, a, ecg_data);



%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating PAT 
%%%%%%%%%%%%%%%%%%%%%%%%%

% find initial sampling rate of ecg 
dt = diff(time_ecg);  % Get time differences
fs_ecg = 1 / mean(seconds(dt));    % Calculate average sampling frequency

% find initial sampling rate of bioz 
nan_indices = isnan(time_bioz); % Remove the NaN values
time_bioz = time_bioz(~nan_indices);
dt = diff(time_bioz);  % Get time differences
fs_bioz = 1 / mean(seconds(dt));    % Calculate average sampling frequency

% original sample rates
round(fs_ecg)  % f = 256
round(fs_bioz) % f = 64


% Resample both signals
fs_resample = 64; % desired sampling rate = bioz sampling rate
ecg_resampled = resample(ecg, fs_resample, round(fs_ecg));
bioz_resampled = resample(bioz, fs_resample, round(fs_bioz));

%{
bioz_resampled = bioz_resampled(x+320:end);
ecg_resampled = ecg_resampled(x+320:end);
%}

% Smooth BioZ with moving average filter
window_size = 20;

% Create the moving average filter
moving_avg_filter = ones(1, window_size) / window_size;

% Apply the filter to the bioimpedance signal
bioz_smoothed = conv(bioz_resampled, moving_avg_filter, 'same');

% find delay between R peak and bioimpedance upstroke 
[acor,lag] = xcorr(ecg_resampled, bioz_smoothed);
[~,I] = max(abs(acor));
delay = lag(I);
bioz_shifted = circshift(bioz_smoothed, [delay 0]);

% calculate pulse transit time (PTT) and pulse arrival time (PAT) 
[pks,r_peak_locs] = findpeaks(ecg_resampled, 'MinPeakHeight', 0.6, 'MinPeakDistance', fs_resample*0.6);

% Find bioimpedance upstroke locations
bioz_diff = diff(bioz_shifted);
bioz_diff_smoothed = smoothdata(bioz_diff, 'gaussian', 0.1 * fs_resample);
bioz_upstroke_locs = find(bioz_diff_smoothed > 0);

% Calculate PTT
ptt = (bioz_upstroke_locs(1:length(r_peak_locs))-r_peak_locs)/fs_resample;

% Calculate PAT
bioz_upstroke_times = bioz_upstroke_locs / fs_resample; % Convert to time
ecg_r_times = r_peak_locs / fs_resample; % Convert to time
pat = ecg_r_times - bioz_upstroke_times(1:length(ecg_r_times));
pat_real = diff(pat);
disp(pat_real);


% Plot results
subplot(3,1,1);
plot(ecg_resampled); hold on;
plot(r_peak_locs, ecg_resampled(r_peak_locs), 'bo', 'MarkerSize', 3);
legend('ECG', 'R peak');
title('ECG Data');
xlabel('Sample Number');
ylabel('Amplitude');
xlim([600, 800])
ylim([-5,5])

subplot(3,1,2);
plot(bioz_resampled); hold on;
plot(bioz_upstroke_locs, bioz_resampled(bioz_upstroke_locs), 'ro', 'MarkerSize', 3);
legend('Bioimpedance', 'Upstroke');
title('Bioimpedance Data');
xlabel('Sample Number');
ylabel('Amplitude');
xlim([600, 800])
ylim([33.2,34.2])

subplot(3,1,3);
plot(pat_real, 'bo-', 'LineWidth', 2);
hold on;
%plot(ptt, 'ro-', 'LineWidth', 2); --> figure out why pat and ppt r same
xlim([1, length(ecg_r_times)]);
legend('PAT');
title('Pulse Arrival Times');
xlabel('Beat Number');
ylabel('Time (s)');








%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimentation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
t = (0:length(bioz_resampled)-1) / fs_resample; % Time vector
start_idx = find(t >= 8, 1); % Index of the first sample after 10 seconds
bioz_cut = bioz_resampled(start_idx:end); % Cut out the portion of the signal



% Find the minimum peak in the bioimpedance signal
[peaks, locs] = findpeaks(bioz_resampled);
max_peak_loc = locs(1); % Location of the minimum peak
max_peak_val = peaks(1); % Value of the minimum peak (flip the sign back to positive)

% Plot the signal and the minimum peak
figure;
plot(bioz_resampled);
hold on;
plot(t(max_peak_loc), max_peak_val, 'ro');
xlabel('Time (s)');
ylabel('Bioimpedance');
legend('Signal', 'Minimum Peak');

%}