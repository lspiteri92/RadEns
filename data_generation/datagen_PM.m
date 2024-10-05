
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constant PRI

pri_class = 1; % constant

wav_ct = zeros([log_limit, seq_length]); %create log
wav_concat_ct = zeros([log_limit*seq_length, 1]); %create log
labels_int = zeros([seq_length, 7]); %create log

% Create Pulse Chain
for k = 1 : n_samples_PM % feed in iterated PM parameters
    for j = 1 : n_samples_PRI % iterate over number of PRI samples %%%%%%%%%%%%%%%% can now increase this to as randomising PM within PRI sequence
        
        % Noise
        noisedB = randi([noise_min noise_max]);
        
        for i = 1 : seq_length % iterate over PRI in single modulation sequence

            pri = constPRI(i,j);
            SampleRate = 1000*(1/pri); 
            NumPulses = LFM_params(k,1);
            Tw = LFM_params(k,2);
            Td = LFM_params(k,3);

            random_wave = randi(3,1); %4 is polyphase
            %random_wave = 4;
            if random_wave == 1
                
                pm_class = 1; % constant

                wav = phased.RectangularWaveform('OutputFormat','Pulses','SampleRate',SampleRate, ...
                'PulseWidth',Tw,'PRF',1/(pri),'NumPulses',NumPulses,'DurationSpecification','Pulse width'); %'FrequencyOffset',1/Td,

                seq_size = size(step(wav),1);
                wav_ct(1:seq_size,i) = step(wav);
                wav_ct(:,i) = awgn(wav_ct(:,i), noisedB); %add noise

            elseif random_wave == 2

                    pm_class = 2; % LFM
                    
                    wav = phased.LinearFMWaveform('OutputFormat','Pulses','SampleRate',SampleRate, ...
                        'PulseWidth',Tw,'PRF',1/(pri),'NumPulses',NumPulses, ...
                        'DurationSpecification','Pulse width'); %create modulated pulse 'FrequencyOffset',1/Td,

                    seq_size = size(step(wav),1);
                    wav_ct(1:seq_size,i) = step(wav);
                    wav_ct(:,i) = awgn(wav_ct(:,i), noisedB); %add noise

            elseif random_wave == 3

                    pm_class = 3; % Phase Code
                    
                    Cw = Tw/NumChips; 

                    wav = phased.PhaseCodedWaveform('SampleRate',SampleRate,'Code',Code,'ChipWidth',round(Cw*1000*(1/pri))/(1000*(1/pri)), ...
                    'NumChips',NumChips, 'OutputFormat','Pulses','PRF',(1/pri),'NumPulses',NumPulses); %'FrequencyOffset',1/Td

                    seq_size = size(step(wav),1);
                    wav_ct(1:seq_size,i) = step(wav);
                    wav_ct(:,i) = awgn(wav_ct(:,i), noisedB); %add noise

            % elseif random_wave == 4
            % 
            %         pm_class = 4; % Poly Phase Code
            % 
            %         Cw = Tw/NumChips; 
            % 
            %         SampPerChip = round(Cw * SampleRate);
            % 
            %         wav = phased.PhaseCodedWaveform('SampleRate',SampleRate,'Code',Code,'ChipWidth',round(Cw*1000*(1/pri))/(1000*(1/pri)), ...
            %         'NumChips',NumChips, 'OutputFormat','Pulses','PRF',(1/pri),'NumPulses',NumPulses); %'FrequencyOffset',1/Td
            % 
            %         sing_wav = step(wav);
            %         seq_size = size(sing_wav,1);
            % 
            %         % Define phase shifts for polyphase components
            %         phaseShifts = linspace(0, 2*pi, numPhases+1);
            %         phaseShifts = phaseShifts(1:end-1); % Exclude the last phase shift to avoid duplication
            % 
            %         % Initialize matrix to hold polyphase components
            %         waveformMatrix = zeros(seq_size, numPhases);
            % 
            %         % Generate each polyphase component
            %         for phaseIdx = 1:numPhases
            %             % Apply phase shift to Barker code
            %             phaseCode = sing_wav .* exp(1j * phaseShifts(phaseIdx));
            % 
            %             % Create waveform for the current phase
            %             phaseWaveform = zeros(seq_size, 1);
            %             for chipIdx = 1:NumChips
            %                 startIdx = (chipIdx - 1) * SampPerChip + 1;
            %                 endIdx = chipIdx * SampPerChip;
            %                 phaseWaveform(startIdx:endIdx) = phas eCode(chipIdx);
            %             end
            % 
            %             % Store the phase waveform in matrix
            %             waveformMatrix(:, phaseIdx) = phaseWaveform;
            %         end
            % 
            %         % Combine the polyphase components
            %         wav = sum(waveformMatrix, 2);
            % 
            %         % Normalize the waveform
            %         wav = wav / max(abs(wav));
            % 
            %         wav_ct(1:seq_size,i) = wav(:,:);
            
            else disp('error in random wave generation')
            end

            labels_int(i,:) = [pm_class, pri, NumPulses, Tw, Td, SampleRate, noisedB]; 
        end

        % Time Delay - Shift Pulse
        shiftWav = zeros(log_limit,seq_length);
        shiftIdx = ceil(Td/(1/SampleRate));
        shiftWav(1:shiftIdx-1,:) = 0;
        shiftWav(shiftIdx:end,:) = wav_ct(1:(end-shiftIdx+1),:);
        
        % Append waveforms together into single column
        wav_db_ct(:,(j-1)*seq_length+1:(j-1)*seq_length+seq_length) = shiftWav(:,:); %create database of waveforms with PM

        pri_db_ct(:,j) = labels_int(:,2); %extract pri only into database
        labels_pm_int((j-1)*seq_length+1:(j-1)*seq_length+seq_length,:) = labels_int(:,:);
        labels_pri_int(j,:) = [pri_class];
      
        wav_ct = zeros([log_limit, seq_length]); %clear log
        wav_concat_ct = zeros([log_limit*seq_length, 1]); %clear log

    end

    wav_db_ct_var(:,((k-1)*n_samples_PRI*seq_length+1:(((k-1)*seq_length*n_samples_PRI+1)+ n_samples_PRI*seq_length -1))) = wav_db_ct(:,:); % 2d array, PRI modulation and then PM stacked column-wise
    labels_ct_pm((k-1)*n_samples_PRI*seq_length+1:((k-1)*n_samples_PRI*seq_length+n_samples_PRI*seq_length),:) = labels_pm_int(:,:); % compile labels for database
    pri_db_ct_var(:,(k-1)*n_samples_PRI+1:((k-1)*n_samples_PRI+n_samples_PRI)) = pri_db_ct(:,:); % compile pri for database
    labels_ct_pri((k-1)*n_samples_PRI+1:((k-1)*n_samples_PRI+n_samples_PRI),:) = labels_pri_int(:,:);

end

%% Jittered PRI

pri_class = 2; % jittered

% confirm if this is the right way to do it - creating a zero value array
% to make sure all samples are of fixed length

wav_jt = zeros([log_limit, seq_length]); %create log
wav_concat_jt = zeros([log_limit*seq_length, 1]); %create log 
labels_int = zeros([seq_length, 7]); %create log

% Create LFM PM Pulses
for k = 1 : n_samples_PM % feed in iterated PM parameters
    for j = 1 : n_samples_PRI % iterate over number of PRI samples
            
        % Noise
        noisedB = randi([noise_min noise_max]);
        
        for i = 1 : seq_length % iterate over PRI in single modulation sequence

            pri = jitteredPRI(i,j);
            SampleRate = 1000*(1/pri); 
            NumPulses = LFM_params(k,1);
            Tw = LFM_params(k,2);
            Td = LFM_params(k,3);
            
            random_wave = randi(3,1);
            if random_wave == 1
                
                pm_class = 1; % constant

                wav = phased.RectangularWaveform('OutputFormat','Pulses','SampleRate',SampleRate, ...
                'PulseWidth',Tw,'PRF',1/(pri),'NumPulses',NumPulses,'DurationSpecification','Pulse width');

            elseif random_wave == 2
                
                pm_class = 2; % LFM
                    
                wav = phased.LinearFMWaveform('OutputFormat','Pulses','SampleRate',SampleRate, ...
                        'PulseWidth',Tw,'PRF',1/(pri),'NumPulses',NumPulses, ...
                        'DurationSpecification','Pulse width'); %create modulated pulse

            elseif random_wave == 3

                pm_class = 3; % Phase Code
                    
                Cw = Tw/NumChips; 

                wav = phased.PhaseCodedWaveform('SampleRate',SampleRate,'Code',Code,'ChipWidth',round(Cw*1000*(1/pri))/(1000*(1/pri)), ...
                    'NumChips',NumChips, 'OutputFormat','Pulses','PRF',(1/pri),'NumPulses',NumPulses); 

            else disp('error in random wave generation')
            end

            SampleRate = wav.SampleRate;
            seq_size = size(step(wav),1);
            wav_jt(1:seq_size,i) = step(wav);

            wav_jt(:,i) = awgn(wav_jt(:,i), noisedB);


            labels_int(i,:) = [pm_class, pri, NumPulses, Tw, Td, SampleRate, noisedB]; 
        end

        % Time Delay - Shift Pulse
        shiftWav = zeros(log_limit,seq_length);
        shiftIdx = ceil(Td/(1/SampleRate));
        shiftWav(1:shiftIdx-1,:) = 0;
        shiftWav(shiftIdx:end,:) = wav_jt(1:(end-shiftIdx+1),:);
        
        % Append waveforms together into single column
        %wav_concat_jt = wav_jt(:); %concat into single vector
        wav_db_jt(:,(j-1)*seq_length+1:(j-1)*seq_length+seq_length) = shiftWav(:,:); %create database of waveforms with PM
        %wav_db_jt(:,j) = wav_concat_jt(:,1); %create database of waveforms with PRI modulation
        pri_db_jt(:,j) = labels_int(:,2); %extract pri only into database
        labels_pm_int((j-1)*seq_length+1:(j-1)*seq_length+seq_length,:) = labels_int(:,:);
        labels_pri_int(j,:) = [pri_class];
        
        wav_jt = zeros([log_limit, seq_length]); %clear log
        wav_concat_jt = zeros([log_limit*seq_length, 1]); %clear log

    end

    wav_db_jt_var(:,((k-1)*n_samples_PRI*seq_length+1:(((k-1)*seq_length*n_samples_PRI+1)+ n_samples_PRI*seq_length -1))) = wav_db_jt(:,:); % 2d array, PRI modulation and then PM stacked column-wise
    labels_jt_pm((k-1)*n_samples_PRI*seq_length+1:((k-1)*n_samples_PRI*seq_length+n_samples_PRI*seq_length),:) = labels_pm_int(:,:); % compile labels for database
    pri_db_jt_var(:,(k-1)*n_samples_PRI+1:((k-1)*n_samples_PRI+n_samples_PRI)) = pri_db_jt(:,:); % compile pri for database
    labels_jt_pri((k-1)*n_samples_PRI+1:((k-1)*n_samples_PRI+n_samples_PRI),:) = labels_pri_int(:,:);
end
   
 
%% Staggered PRI

pri_class = 3; % staggered

wav_st = zeros([log_limit, seq_length]); %create log
wav_concat_st = zeros([log_limit*seq_length, 1]); %create log
labels_int = zeros([seq_length, 7]); %create log

% Create LFM PM Pulses
for k = 1 : n_samples_PM % feed in iterated PM parameters
    for j = 1 : n_samples_PRI % iterate over number of PRI samples
       
        % Noise
        noisedB = randi([noise_min noise_max]);

        for i = 1 : seq_length % iterate over PRI in single modulation sequence

            pri = staggeredPRI(i,j);
            SampleRate = 1000*(1/pri); 
            NumPulses = LFM_params(k,1);
            Tw = LFM_params(k,2);
            Td = LFM_params(k,3);

            random_wave = randi(3,1);

            if random_wave == 1
                
                pm_class = 1; % constant

                wav = phased.RectangularWaveform('OutputFormat','Pulses','SampleRate',SampleRate, ...
                'PulseWidth',Tw,'PRF',1/(pri),'NumPulses',NumPulses,'DurationSpecification','Pulse width');

            elseif random_wave == 2

                    pm_class = 2; % LFM
                    
                    wav = phased.LinearFMWaveform('OutputFormat','Pulses','SampleRate',SampleRate, ...
                        'PulseWidth',Tw,'PRF',1/(pri),'NumPulses',NumPulses, ...
                        'DurationSpecification','Pulse width'); %create modulated pulse

            elseif random_wave == 3

                    pm_class = 3; % Phase Code
                    
                    Cw = Tw/NumChips; 

                    wav = phased.PhaseCodedWaveform('SampleRate',SampleRate,'Code',Code,'ChipWidth',round(Cw*1000*(1/pri))/(1000*(1/pri)), ...
                    'NumChips',NumChips, 'OutputFormat','Pulses','PRF',(1/pri),'NumPulses',NumPulses); 

            else disp('error in random wave generation')
            end

            SampleRate = wav.SampleRate;
            seq_size = size(step(wav),1);
            wav_st(1:seq_size,i) = step(wav);
            wav_st(:,i) = awgn(wav_st(:,i), noisedB);

            labels_int(i,:) = [pm_class, pri, NumPulses, Tw, Td, SampleRate, noisedB]; 
        end

        % Time Delay - Shift Pulse
        shiftWav = zeros(log_limit,seq_length);
        shiftIdx = ceil(Td/(1/SampleRate));
        shiftWav(1:shiftIdx-1,:) = 0;
        shiftWav(shiftIdx:end,:) = wav_st(1:(end-shiftIdx+1),:);
        
        % Append waveforms together into single column
        %wav_concat_st = wav_st(:); %concat into single vector
        wav_db_st(:,(j-1)*seq_length+1:(j-1)*seq_length+seq_length) = shiftWav(:,:); %create database of waveforms with PM
        %wav_db_st(:,j) = wav_concat_st(:,1); %create database of waveforms with PRI modulation
        labels_pm_int((j-1)*seq_length+1:(j-1)*seq_length+seq_length,:) = labels_int(:,:);
        pri_db_st(:,j) = labels_int(:,2); %extract pri only into database
        labels_pri_int(j,:) = [pri_class];

        wav_st = zeros([log_limit, seq_length]); %clear log
        wav_concat_st = zeros([log_limit*seq_length, 1]); %clear log

    end

    wav_db_st_var(:,((k-1)*n_samples_PRI*seq_length+1:(((k-1)*seq_length*n_samples_PRI+1)+ n_samples_PRI*seq_length -1))) = wav_db_st(:,:); % 2d array, PRI modulation and then PM stacked column-wise
    labels_st_pm((k-1)*n_samples_PRI*seq_length+1:((k-1)*n_samples_PRI*seq_length+n_samples_PRI*seq_length),:) = labels_pm_int(:,:); % compile labels for database
    pri_db_st_var(:,(k-1)*n_samples_PRI+1:((k-1)*n_samples_PRI+n_samples_PRI)) = pri_db_st(:,:); % compile pri for database
    labels_st_pri((k-1)*n_samples_PRI+1:((k-1)*n_samples_PRI+n_samples_PRI),:) = labels_pri_int(:,:);
end

