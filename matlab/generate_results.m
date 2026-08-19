function generate_results(varargin)
% =============================================================================
% Flexible driver for channel discrimination hierarchy SDPs.
% Supports unitary + non-unitary channels.
% =============================================================================

% Print which physical copy of this file is actually running. If you ever
% see a "not a recognized parameter" error again, check this line first --
% it almost always means a stale/older generate_results.m earlier on the
% MATLAB path is shadowing this one (run `which -all generate_results` to
% list every copy MATLAB can see).
fprintf('generate_results running from: %s\n', mfilename('fullpath'));

% =========================================================================
% PARSE INPUTS
% =========================================================================
p = inputParser();
p.addParameter('protocols',    [1 2 9],   @isnumeric);
p.addParameter('N_start',      2,         @isnumeric);
p.addParameter('N_max',        20,        @isnumeric);
p.addParameter('n_samples',    1,         @isnumeric);
p.addParameter('seed',         42);
p.addParameter('prior',        'uniform', @ischar);
p.addParameter('channel_type', 'unitary', @ischar);
p.addParameter('channels',     'diff',    @ischar);  % Renamed from unitary_mode
p.addParameter('unitary_mode', '',        @ischar);  % Backward-compat alias; see below
p.addParameter('commuting',    0,         @(x) isnumeric(x) && ismember(x, [0 2 3]));
p.addParameter('savedir',      fullfile(getenv('HOME'), 'MATLAB', 'results'), @ischar);
p.addParameter('run_id',       '',        @ischar);  % caller-supplied run folder name
p.parse(varargin{:});
cfg = p.Results;

% Backward compatibility: if an older caller still passes 'unitary_mode'
% instead of 'channels', honor it instead of silently ignoring it.
if ~isempty(cfg.unitary_mode)
    cfg.channels = cfg.unitary_mode;
end

fprintf('cfg.commuting = %d\n', cfg.commuting);

if ~ismember(cfg.channel_type, {'unitary','channels','pauli','superoperator'})
    error('channel_type must be ''unitary'', ''channels'', or ''pauli''.');
end

if cfg.commuting > 0 && ~strcmpi(cfg.channel_type, 'unitary')
    warning('commuting=%d is only supported for channel_type=''unitary''. Ignoring.', cfg.commuting);
end

% =========================================================================
% SETUP
% =========================================================================
home = getenv('HOME');
addpath(genpath(fullfile(home, 'MATLAB', 'Library')));
addpath(genpath(fullfile(home, 'MATLAB', 'SDP')));
addpath(genpath(pwd));

cvx_quiet true;

% Mosek license
mac_lic = fullfile(home, 'mosek', 'mosek.lic');
lin_lic = fullfile(home, 'MATLAB', 'Library', 'mosek', 'mosek.lic');
if exist(mac_lic, 'file')
    setenv('MOSEKLM_LICENSE_FILE', mac_lic);
elseif exist(lin_lic, 'file')
    setenv('MOSEKLM_LICENSE_FILE', lin_lic);
end

% Protocol names
proto_names = containers.Map( ...
    {1,2,3,4,5,6,7,8,9}, ...
    {'PAR','SEQ','QC-convFO','QC-NICC','pQC-CC','QC-SupFO','QC-NIQC','QC-QC','GEN'});

d = 2; k = 3;
dim = d^(2*k);
dIn = d; dOut = d;

protocols  = cfg.protocols;
N_start    = cfg.N_start;
N_max      = cfg.N_max;
n_samples  = cfg.n_samples;

n_proto = numel(protocols);
n_N     = N_max - N_start + 1;
N_list  = (N_start:N_max)';

% Labels
proto_labels = cell(1, n_proto);
for j = 1:n_proto
    proto_labels{j} = proto_names(protocols(j));
end

% RNG
if ischar(cfg.seed) && strcmpi(cfg.seed,'random')
    rng('shuffle');
    s_ = rng; seed_used = s_.Seed; seed_str = 'random';
else
    seed_used = cfg.seed; rng(seed_used);
    seed_str = sprintf('seed%d', seed_used);
end

% -------------------------------------------------------------------------
% Output folder
% -------------------------------------------------------------------------
tag     = strjoin(proto_labels,'-');
basedir = cfg.savedir;
if basedir(end) == '/'; basedir = basedir(1:end-1); end

if ~isempty(cfg.run_id)
    run_id = cfg.run_id;
else
    ts     = datestr(now,'yyyymmdd_HHMMSS');
    run_id = [ts '_' tag '_' cfg.channel_type '_' cfg.channels '_' num2str(n_samples) 'samples'];
end

rundir = [basedir '/' run_id];
system(['mkdir -p "' rundir '"']);

% Modes to loop over
if strcmpi(cfg.channels,'both')
    modes = {'diff','same'};
else
    modes = {cfg.channels};
end

% =========================================================================
% PRINT CONFIG
% =========================================================================
fprintf('\n%s\n', repmat('=',1,72));
fprintf('generate_results\n');
fprintf('%s\n', repmat('=',1,72));
fprintf('Protocols    : %s\n', strjoin(proto_labels, ', '));
fprintf('N range      : %d .. %d\n', N_start, N_max);
fprintf('n_samples    : %d\n', n_samples);
fprintf('channel_type : %s\n', cfg.channel_type);
fprintf('channels     : %s\n', cfg.channels);
fprintf('run_id       : %s\n', run_id);
fprintf('save to      : %s\n', rundir);
fprintf('%s\n\n', repmat('=', 1, 72));

write_meta_json(rundir, cfg, run_id, seed_used, seed_str, proto_labels);

% =========================================================================
% MAIN LOOP
% =========================================================================
for m = 1:numel(modes)

    mode_str = modes{m};
    fprintf('--- %s | %s ---\n\n', cfg.channel_type, mode_str);

    results  = zeros(n_N, n_proto, n_samples);
    times    = zeros(n_N, n_proto, n_samples);
    statuses = cell(n_N, n_proto, n_samples);

    % Reproducibility store
    channels_used = cell(n_N, n_samples);

    for idx = 1:n_N
        N = N_start + idx - 1;

        for s = 1:n_samples

            % Prior
            if strcmpi(cfg.prior, 'uniform')
                p_i = ones(1, N) / N;
            else
                e   = -log(rand(1, N));
                p_i = e / sum(e);
            end

            % -----------------------------------------------------------------
            % Build C3  (dim x dim x N)
            % -----------------------------------------------------------------
            C3  = zeros(dim, dim, N);
            Ch1_store = zeros(d*d, d*d, N);
            Ch2_store = zeros(d*d, d*d, N);
            Ch3_store = zeros(d*d, d*d, N);

            if strcmpi(mode_str, 'identical')
                for i = 1:N
                    Ch = draw_one_channel(cfg.channel_type, d);
                    C3(:,:,i) = kron(kron(Ch, Ch), Ch);
                    Ch1_store(:,:,i) = Ch;
                    Ch2_store(:,:,i) = Ch;
                    Ch3_store(:,:,i) = Ch;
                end
            else
                for i = 1:N
                    if strcmpi(mode_str, 'same')
                        Ch  = draw_one_channel(cfg.channel_type, d);
                        Ch1 = Ch; Ch2 = Ch; Ch3_slot = Ch;
                    else
                        Ch1      = draw_one_channel(cfg.channel_type, d);
                        Ch2      = draw_one_channel(cfg.channel_type, d);
                        Ch3_slot = draw_one_channel(cfg.channel_type, d);
                    end
                    C3(:,:,i) = kron(kron(Ch1, Ch2), Ch3_slot);
                    Ch1_store(:,:,i) = Ch1;
                    Ch2_store(:,:,i) = Ch2;
                    Ch3_store(:,:,i) = Ch3_slot;
                end
            end

            channels_used{idx, s} = struct( ...
                'N', N, 'sample', s, ...
                'Ch1', Ch1_store, 'Ch2', Ch2_store, 'Ch3', Ch3_store);

            if N >= 2
                fprintf('C3(1) vs C3(2) diff: %.6e\n', norm(C3(:,:,1) - C3(:,:,2),'fro'));
            end
            fprintf('commuting = %d\n', cfg.commuting);

            proto_results = zeros(1,n_proto);
            proto_times   = zeros(1,n_proto);
            proto_status  = cell(1,n_proto);

            for j = 1:n_proto
                t_sdp = tic;
                [pS,~,~,st] = channel_discrimination_3copies_primal_9classes( ...
                    C3, protocols(j), [dIn dOut], p_i);
                proto_times(j)   = toc(t_sdp);
                proto_results(j) = real(pS);
                proto_status{j}  = st;
            end

            cvx_clear;  

            results(idx,:,s)  = proto_results;
            times(idx,:,s)    = proto_times;
            statuses(idx,:,s) = proto_status;

            % Print
            fprintf('N=%2d s=%d/%d [%s|%s]\n', ...
                N, s, n_samples, cfg.channel_type, mode_str);

            for j = 1:n_proto
                fprintf('  %-10s p=%.6f t=%.1fs [%s]\n', ...
                    proto_labels{j}, proto_results(j), ...
                    proto_times(j), proto_status{j});
            end
            fprintf('\n');

        end
    end

    % =========================================================================
    % SAVE (MAT files)
    % =========================================================================
    avg_results = mean(results,3);
    avg_times   = mean(times,3);

    save([rundir '/data_' mode_str '.mat'], ...
        'results','avg_results','times','avg_times','statuses', ...
        'N_list','proto_labels','protocols','cfg','run_id','-v7.3');

    save([rundir '/channels_used_' mode_str '.mat'], ...
        'channels_used','N_list','n_samples','d','k','run_id','mode_str','-v7.3');

    % =========================================================================
    % SAVE (TXT files)
    % =========================================================================
    fid_res = fopen([rundir '/results_' mode_str '.txt'], 'w');
    fprintf(fid_res, '============================================================\n');
    fprintf(fid_res, 'Simulation Results Summary\n');
    fprintf(fid_res, 'Run ID      : %s\n', run_id);
    fprintf(fid_res, 'Channel Type: %s\n', cfg.channel_type);
    fprintf(fid_res, 'Mode        : %s\n', mode_str);
    fprintf(fid_res, '============================================================\n\n');
    
    fprintf(fid_res, '--- AVERAGE SUCCESS PROBABILITIES ---\n');
    fprintf(fid_res, '%-5s', 'N');
    for j = 1:n_proto
        fprintf(fid_res, '\t%-12s', proto_labels{j});
    end
    fprintf(fid_res, '\n');
    for idx = 1:n_N
        N_val = N_list(idx);
        fprintf(fid_res, '%-5d', N_val);
        for j = 1:n_proto
            fprintf(fid_res, '\t%12.8f', avg_results(idx, j));
        end
        fprintf(fid_res, '\n');
    end
    
    fprintf(fid_res, '\n--- AVERAGE COMPUTATION TIMES (seconds) ---\n');
    fprintf(fid_res, '%-5s', 'N');
    for j = 1:n_proto
        fprintf(fid_res, '\t%-12s', proto_labels{j});
    end
    fprintf(fid_res, '\n');
    for idx = 1:n_N
        N_val = N_list(idx);
        fprintf(fid_res, '%-5d', N_val);
        for j = 1:n_proto
            fprintf(fid_res, '\t%12.4f', avg_times(idx, j));
        end
        fprintf(fid_res, '\n');
    end
    fclose(fid_res);

    fid_ch = fopen([rundir '/channels_used_' mode_str '.txt'], 'w');
    fprintf(fid_ch, '============================================================\n');
    fprintf(fid_ch, 'Choi Matrices for Run ID: %s\n', run_id);
    fprintf(fid_ch, 'Mode                    : %s\n', mode_str);
    fprintf(fid_ch, '============================================================\n');
    fprintf(fid_ch, 'Format per matrix: %dx%d complex values (Real + Imag*i)\n\n', d*d, d*d);
    
    for idx = 1:n_N
        N_val = N_list(idx);
        for s = 1:n_samples
            struct_data = channels_used{idx, s};
            
            fprintf(fid_ch, '=========================================\n');
            fprintf(fid_ch, 'N = %d | Sample = %d\n', N_val, s);
            fprintf(fid_ch, '=========================================\n\n');
            
            for h = 1:N_val
                fprintf(fid_ch, '--- Hypothesis %d ---\n', h);
                
                fprintf(fid_ch, '  -- Ch1 --\n');
                print_matrix_to_file(fid_ch, struct_data.Ch1(:,:,h));
                
                fprintf(fid_ch, '  -- Ch2 --\n');
                print_matrix_to_file(fid_ch, struct_data.Ch2(:,:,h));
                
                fprintf(fid_ch, '  -- Ch3 --\n');
                print_matrix_to_file(fid_ch, struct_data.Ch3(:,:,h));
                fprintf(fid_ch, '\n');
            end
        end
    end
    fclose(fid_ch);

    fprintf('Saved results + channels for mode: %s (.mat and .txt formats)\n\n', mode_str);

end

fprintf('\n=== DONE ===\n');
fprintf('RUN_ID:%s\n', run_id);   
end

% =========================================================================
% LOCAL HELPERS
% =========================================================================
function C = unitary_to_choi(U, d)
    Id    = eye(d);
    UId   = kron(U, Id);
    omega = reshape(Id, [d*d, 1]);
    C     = UId * (omega * omega') * UId';
end

function print_matrix_to_file(fid, M)
    [rows, cols] = size(M);
    for r = 1:rows
        fprintf(fid, '    ');
        for c = 1:cols
            val = M(r,c);
            if imag(val) >= 0
                fprintf(fid, '%11.8f + %11.8fi\t', real(val), imag(val));
            else
                fprintf(fid, '%11.8f - %11.8fi\t', real(val), abs(imag(val)));
            end
        end
        fprintf(fid, '\n');
    end
end

function write_meta_json(rundir, cfg, run_id, seed_used, seed_str, proto_labels)
    fid = fopen([rundir '/meta.json'], 'w');
    fprintf(fid, '{\n');
    fprintf(fid, '  "run_id": "%s",\n', run_id);
    fprintf(fid, '  "job_id": "%s",\n', getenv('SLURM_JOB_ID'));
    fprintf(fid, '  "node": "%s",\n', getenv('HOSTNAME'));
    fprintf(fid, '  "protocols": [%s],\n', strjoin(cellstr(num2str(cfg.protocols(:))), ', '));
    fprintf(fid, '  "proto_labels": ["%s"],\n', strjoin(proto_labels, '", "'));
    fprintf(fid, '  "N_start": %d,\n', cfg.N_start);
    fprintf(fid, '  "N_max": %d,\n', cfg.N_max);
    fprintf(fid, '  "n_samples": %d,\n', cfg.n_samples);
    fprintf(fid, '  "channel_type": "%s",\n', cfg.channel_type);
    fprintf(fid, '  "channels": "%s",\n', cfg.channels);
    fprintf(fid, '  "commuting": %d,\n', cfg.commuting);
    fprintf(fid, '  "prior": "%s",\n', cfg.prior);
    fprintf(fid, '  "seed": "%s",\n', seed_str);
    fprintf(fid, '  "seed_used": %s,\n', mat2str(seed_used));
    fprintf(fid, '  "source": "generate_results.m"\n');
    fprintf(fid, '}\n');
    fclose(fid);
end