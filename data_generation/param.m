
rng = 32; %set random seed

seq_length = 40; % number of waves in PRI modulation sequence - 10

% Database parameters
n_samples_PRI = 80; % number of samples of each PRI modulation class - 80
n_samples_PM = 10; % number of samples of each pulse modulation class - 10

% PRI characteristics
PRI_min = 10; % Minimum PRI (ms)
PRI_max = 30; % Maximum PRI (ms)

% Waveform characteristics
n_pulses_min = 1; % Minimum number of pulses
n_pulses_max = 3; % Maximum number of pulses - was 6
Tw_min = 10; % Mininum Pulse Width - Units: ms
Tw_max = 16; % Maximum Pulse Width - Units: ms
Td_min = 1; % Minimum Pulse Time Delay (Frequency Offset) - Units: ms
Td_max = 30; % Maximum Pulse Time Delay (Frequency Offset) - Units: ms


% Vector Size
log_limit = 3000; % vector limit for waveform - TBC

%% PRI MODULATION %%

% Constant PRI %

% Jittered PRI %
min_dev = 3;
max_dev = 30;

% Staggered PRI %
min_len = 1;
max_len = 10; %seq_length

%% PM %%

numPhases = 4; %Polyphase only
Code = 'Barker'; %'Frank','P1','P2','Barker'
NumChips = 13; %Up to 13 for Barker Code / 16 for Frank Code

%% White Noise
noise_min = 0; %dB
noise_max = 20; %dB 20

%% Create Modulation parameter database

%Note - Tw set to always be less than allowable PRI
for i = 1 : n_samples_PM
    n_pulses = randi([n_pulses_min n_pulses_max]); % randomise number of pulses
    Tw = randi([Tw_min Tw_max])*1e-6; % Pulse Width
    Td = randi([Td_min Td_max])*1e-6; % Pulse Time Delay (Frequency Offset)

    LFM_params(i,:) = [n_pulses, Tw, Td]; % Create LFM database 
end