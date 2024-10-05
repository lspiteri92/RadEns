%%

% Combine all labels

% class 1 - ct_UM
% class 2 - ct_LFM
% class 3 - ct_PC
% class 4 - jt_UM
% class 5 - jt_LFM
% class 6 - jt_PC
% class 7 - st_UM
% class 8 - st_LFM
% class 9 - st_PC

%labels_comb = [labels_ct_UM; labels_ct_LFM; labels_ct_PC; labels_jt_UM; labels_jt_LFM; labels_jt_PC; labels_st_UM; labels_st_LFM; labels_st_PC];
%labels_comb = [labels_ct_var; labels_jt_var; labels_st_var];
comb_labels_pri = [labels_ct_pri; labels_jt_pri; labels_st_pri];
comb_labels_pm = [labels_ct_pm; labels_jt_pm; labels_st_pm];

% %% UID Counter
% 
% %pm
% increment = seq_length;
% numElements = size(labels_comb_pm,1)+(seq_length-2);
% uID_pm = zeros(numElements,1);
% 
% currentVal = 0;
% 
% for i = 1:numElements+1
%     if mod(i, increment) == 0
%         currentVal = currentVal+1;
%     end
%     uID_pm(i) = currentVal;
% end
% 
% uID_pm = uID_pm(seq_length:end);
% 
% labels_comb_pm = [uID_pm, labels_comb_pm];
% 
% 
% %pri
% increment = 1;
% numElements = size(labels_comb_pri,1);
% uID_pri = zeros(numElements,1);
% currentVal = 0;
% 
% for i = 1:numElements
%     if mod(i, increment) == 0
%         currentVal = currentVal+1;
%     end
%     uID_pri(i) = currentVal;
% end
% 
% uID_pri = uID_pri(1:end);
% 
% labels_comb_pri = [uID_pri, labels_comb_pri];


%%

% Combine all wave data

%wav_comb = [wav_db_ct_UM, wav_db_ct_LFM, wav_db_ct_PC, wav_db_jt_UM, wav_db_jt_LFM, wav_db_jt_PC, wav_db_st_UM, wav_db_st_LFM, wav_db_st_PC];
comb_wav = [wav_db_ct_var, wav_db_jt_var, wav_db_st_var];
comb_pri = [pri_db_ct_var, pri_db_jt_var, pri_db_st_var];

%% 

% Write to file - matlab requires you to split real and imaginary parts
% wav_r_ct = real(wav_db_ct_var);
% wav_i_ct = imag(wav_db_ct_var);
% wav_r_jt = real(wav_db_jt_var);
% wav_i_jt = imag(wav_db_jt_var);
% wav_r_st = real(wav_db_st_var);
% wav_i_st = imag(wav_db_st_var);
comb_wav_r = real(comb_wav);
comb_wav_i = imag(comb_wav);
% 
% h5create('dataset_test.h5', '/wav_r_ct',size(wav_r_ct));
% h5write('dataset_test.h5','/wav_r_ct',wav_r_ct);
% h5create('dataset_test.h5', '/wav_i_ct',size(wav_i_ct));
% h5write('dataset_test.h5','/wav_i_ct',wav_i_ct);
% 
% h5create('dataset_test.h5', '/wav_r_jt',size(wav_r_jt));
% h5write('dataset_test.h5','/wav_r_jt',wav_r_jt);
% h5create('dataset_test.h5', '/wav_i_jt',size(wav_i_jt));
% h5write('dataset_test.h5','/wav_i_jt',wav_i_jt);
% 
% h5create('dataset_test.h5', '/wav_r_st',size(wav_r_st));
% h5write('dataset_test.h5','/wav_r_st',wav_r_st);
% h5create('dataset_test.h5', '/wav_i_st',size(wav_i_st));
% h5write('dataset_test.h5','/wav_i_st',wav_i_st);

h5create('dataset_2.h5', '/comb_wav_r',size(comb_wav_r));
h5write('dataset_2.h5','/comb_wav_r',comb_wav_r);
h5create('dataset_2.h5', '/comb_wav_i',size(comb_wav_i));
h5write('dataset_2.h5','/comb_wav_i',comb_wav_i);

h5create('dataset_2.h5', '/comb_pri',size(comb_pri));
h5write('dataset_2.h5','/comb_pri',comb_pri);

h5create('dataset_2.h5', '/labels_pri',size(comb_labels_pri));
h5write('dataset_2.h5','/labels_pri',comb_labels_pri);
h5create('dataset_2.h5', '/labels_pm',size(comb_labels_pm));
h5write('dataset_2.h5','/labels_pm',comb_labels_pm);