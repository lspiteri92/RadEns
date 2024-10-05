%% Constant PRI %%

for i = 1 : n_samples_PRI
    PRI_const = randi([PRI_min PRI_max]);
    for j = 1 : seq_length % build out pulse train
        constPRI(j,i) = [PRI_const]; %in rows rather than columns
    end
end

constPRI = constPRI.*10e-6;


%% Jittered PRI %%

for i = 1 : n_samples_PRI
    PRI_jit_0 = randi([PRI_min PRI_max]);
    nRows = 1; % only sample 1 row of data
    for j = 1 : seq_length % build out pulse train
        dev = randi([min_dev max_dev])/100;
        PRI_jit_min = PRI_jit_0 * (1-dev);
        PRI_jit_max = PRI_jit_0 * (1+dev);
        data = [PRI_jit_min; PRI_jit_max];

        rowsToExtract = randperm(size(data, 1), nRows);
        randomlySelectedRows = data(rowsToExtract, :); %sample from jitter
        jitteredPRI(j,i) = [randomlySelectedRows];
    end
end

jitteredPRI = jitteredPRI.*10e-6;


%% Staggered PRI %%

for i = 1 : n_samples_PRI
    M = randi([min_len max_len]); %length of sequence
    num_seq = ceil(seq_length / M);
    for j = 1 : M
        PRI = randi([PRI_min PRI_max]);
        PRI_st(j,1) = [PRI]; % stagger PRI sequence
    end
    PRI_st_sq = repmat(PRI_st, [num_seq 1]);
    PRI_st_cut = PRI_st_sq(1:seq_length);
    staggeredPRI(:,i) = PRI_st_cut(:,1);
end

staggeredPRI = staggeredPRI.*10e-6;
