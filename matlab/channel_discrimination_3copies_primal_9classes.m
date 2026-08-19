function [pS, T, W, status] = channel_discrimination_3copies_primal_9classes(C, protocol, varargin)
% =========================================================================
% Primal SDP for k=3 slot channel discrimination.
%
% PROTOCOLS (full hierarchy, low -> high):
%   1 = PAR       Parallel
%   2 = SEQ       Sequential
%   3 = QC-convFO Convex mixture of fixed causal orders  (Eq. 5.15)
%   4 = QC-NICC   Non-influenceable classical control    (Eq. 5.17)
%   5 = pQC-CC    Probabilistic classical control
%   6 = QC-SupFO  Superposition of fixed causal orders   (Eq. 5.20)
%   7 = QC-NIQC   Non-influenceable quantum control      (Eq. 5.22)
%   8 = QC-QC     Quantum control of causal order
%   9 = GEN       General / unconstrained
%
% Reference: Bavaresco, Murao & Quintino, arXiv:2105.13369
% Solver: CVX + Mosek
%
% INPUT:
%   C           -- d^6 x d^6 x N  full 3-slot Choi matrix (built by caller)
%                  OR d^2 x d^2 x N single-slot Choi (Tensor called internally)
%   protocol    -- integer 1..9
%   varargin{1} -- [dIn dOut]  (default: d=2)
%   varargin{2} -- prior p_i   (default: uniform)
%
% OUTPUT:
%   pS     -- optimal success probability
%   T      -- d^6 x d^6 x N  per-hypothesis testers
%   W      -- d^6 x d^6      total process matrix
%   status -- CVX status string
% =========================================================================

N = size(C, 3);
k = 3;

switch length(varargin)
    case 0
        d    = 2;
        dIn  = d;  dOut = d;
        p_i  = ones(1,N) / N;
    case 1
        dIn  = varargin{1}(1);  dOut = varargin{1}(2);
        d    = dIn;
        p_i  = ones(1,N) / N;
    case 2
        dIn  = varargin{1}(1);  dOut = varargin{1}(2);
        d    = dIn;
        p_i  = varargin{2};
    otherwise
        error('Too many input arguments.');
end

dim  = dIn^(2*k);   % 64 for d=2, k=3
DIM  = [dIn dOut dIn dOut dIn dOut];

% Accept 64x64xN (pre-built) or 4x4xN (single-slot, Tensor called here)
C3 = zeros(dim, dim, N);
if size(C,1) == dim
    C3 = C;
else
    for i = 1:N
        C3(:,:,i) = Tensor(C(:,:,i), k);
    end
end

dim3  = [d d d];
dim4  = [d d d d];
dim5  = [d d d d d];
dim6  = [d d d d d d];
dim4r = [d d d d];

N_loc   = double(N);
d_loc   = double(d);
dim_loc = double(dim);

% =========================================================================
switch protocol

% =========================================================================
% PROTOCOL 1: PAR
% Source: channel_discrimination_3copies_primal.m, case 1
% =========================================================================
case 1

% CVX_CLEAR   Clears all active CVX data.
%    CVX_CLEAR clears the current CVX model in progress. This is useful if, for
%    example, you have made an error typing in your model and wish to start 
%    over. Typing this before entering another CVX_BEGIN again avoids the 
%    warning message that occurs if CVX_BEGIN detects a model in progress.
cvx_clear

cvx_begin SDP
cvx_solver mosek

    variable T(dim, dim, N) complex semidefinite
    expression pS; pS = 0;
    expression W;  W  = 0;
    for i = 1:N
        % CVX/Mosek works in floaring point airthmetic. The imaginary part
        % of the result is never exactly 0 (1 + 3.2e-16i). Both T_i and C_i
        % are Hermitian Positive Semidefinite. CVX's "maximise" requires a
        % real-valued expression
        pS = pS + real(trace(p_i(i) * T(:,:,i) * C3(:,:,i))); % is the probability of success
        W  = W  + T(:,:,i); 
    end
    % DIM  = [dIn dOut dIn dOut dIn dOut];
    % Pos  = [1   2    3   4    5   6]
    W == TR(W, [2 4 6], DIM); % W must be independent of all output systems
    trace(W) == dOut^k;
    maximise real(pS);
cvx_end
status = cvx_status;

% =========================================================================
% PROTOCOL 2: SEQ
% Source: channel_discrimination_3copies_primal.m, case 2
% =========================================================================
case 2
% =========================================================================
% PROTOCOL 2: SEQ — optimal over ALL 6 sequential orderings
%
% A sequential (comb) strategy fixes one ordering sigma in S3 and applies
% the 3 comb conditions in that ordering's natural slot order.
% The true SEQ optimal is:
%   p^S_SEQ = max_{sigma in S3} p^S_{SEQ,sigma}
%
% The 6 orderings and their comb conditions (slots always refer to
% sigma-natural order [I_{k1},O_{k1},I_{k2},O_{k2},I_{k3},O_{k3}] = 1..6):
%   Cond 1: W = TR(W, 6, DIM_sigma)          [last output trivial]
%   Cond 2: PT(W,[5,6]) = kron(PT(W,[4,5,6]), eye(d)/d)
%   Cond 3: PT(W,[3,4,5,6]) = kron(PT(W,[2,3,4,5,6]), eye(d)/d)
%
% For each sigma, the SDP is solved with C3 permuted to sigma-natural order,
% then the tester is permuted back to canonical order for output.
% =========================================================================

% All 6 orderings of {1,2,3} and their permutation vectors
% perm_vec maps sigma-natural -> canonical I1O1I2O2I3O3
orderings = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];
perm_vecs = [1 2 3 4 5 6;   % (1,2,3): identity
             1 2 5 6 3 4;   % (1,3,2): I1O1 I3O3 I2O2
             3 4 1 2 5 6;   % (2,1,3): I2O2 I1O1 I3O3
             5 6 1 2 3 4;   % (2,3,1): I2O2 I3O3 I1O1 -- wait, 231 means k1=2,k2=3,k3=1
             3 4 5 6 1 2;   % (3,1,2): I3O3 I1O1 I2O2
             5 6 3 4 1 2];  % (3,2,1): I3O3 I2O2 I1O1

% Permutation vectors: sigma-natural order -> canonical I1O1I2O2I3O3
% For sigma=(k1,k2,k3), sigma-natural slot j holds I_{kj} or O_{kj}.
% To go from sigma-natural to canonical, canonical slot 2*ki-1 = sigma slot 2*pos(ki)-1
% These are precomputed correctly:
%   sigma=(1,2,3): canonical = sigma-natural -> [1 2 3 4 5 6]
%   sigma=(1,3,2): sigma-natural=[I1,O1,I3,O3,I2,O2], canonical needs
%                  I1->pos1, O1->pos2, I2->pos5, O2->pos6, I3->pos3, O3->pos4
%                  so PermuteSystems with [1 2 5 6 3 4] maps natural->canonical
%   sigma=(2,1,3): sigma-natural=[I2,O2,I1,O1,I3,O3]
%                  I1->slot3, O1->slot4, I2->slot1, O2->slot2 -> perm [3 4 1 2 5 6]
%   sigma=(2,3,1): sigma-natural=[I2,O2,I3,O3,I1,O1]
%                  I1->slot5, O1->slot6, I2->slot1, O2->slot2 -> perm [5 6 1 2 3 4]
%   sigma=(3,1,2): sigma-natural=[I3,O3,I1,O1,I2,O2]
%                  I1->slot3, O1->slot4, I2->slot5, O2->slot6, I3->slot1 -> perm [3 4 5 6 1 2]
%   sigma=(3,2,1): sigma-natural=[I3,O3,I2,O2,I1,O1]
%                  I1->slot5, O1->slot6, I2->slot3, O2->slot4, I3->slot1 -> perm [5 6 3 4 1 2]

best_pS  = -Inf;
best_T   = zeros(dim, dim, N);
best_W   = zeros(dim, dim);
best_st  = '';

for ord_idx = 1:6

    sigma    = orderings(ord_idx, :);   % e.g. [2 3 1]
    pvec     = perm_vecs(ord_idx, :);   % permutation vector to canonical

    % Permute C3 hypotheses from canonical to sigma-natural order
    % (apply inverse permutation: canonical -> sigma-natural)
    inv_pvec = zeros(1,6);
    inv_pvec(pvec) = 1:6;              % inverse of pvec
    C3_sigma = zeros(dim, dim, N);
    for i = 1:N
        C3_sigma(:,:,i) = PermuteSystems(C3(:,:,i), inv_pvec, dim6);
    end

    cvx_clear
    cvx_begin SDP quiet
    cvx_solver mosek

        variable T_sigma(dim, dim, N) complex semidefinite
        expression pS_sigma; pS_sigma = 0;
        expression W_sigma;  W_sigma  = 0;

        for i = 1:N
            pS_sigma = pS_sigma + real(trace(p_i(i) * T_sigma(:,:,i) * C3_sigma(:,:,i)));
            W_sigma  = W_sigma  + T_sigma(:,:,i);
        end

        % SEQ comb conditions in sigma-natural order (slots 1..6)
        % Cond 1: last output A^O_{k3} (slot 6) is trivial
        W_sigma == TR(W_sigma, [6], DIM);

        % Cond 2: A^O_{k2} (slot 4) cannot be influenced by slot 3 onwards
        PartialTrace(W_sigma, [6 5],     dim6) == kron(PartialTrace(W_sigma, [6 5 4],     dim6), eye(d)/d);

        % Cond 3: A^O_{k1} (slot 2) cannot be influenced by slots 3 onwards
        PartialTrace(W_sigma, [6 5 4 3], dim6) == kron(PartialTrace(W_sigma, [6 5 4 3 2], dim6), eye(d)/d);

        trace(W_sigma) == dOut^k;
        maximise real(pS_sigma);

    cvx_end

    val = real(pS_sigma);

    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        if val > best_pS
            best_pS = val;
            best_st = cvx_status;

            % Map testers back to canonical order
            for i = 1:N
                best_T(:,:,i) = PermuteSystems(T_sigma(:,:,i), pvec, dim6);
            end
            best_W = zeros(dim, dim);
            for i = 1:N
                best_W = best_W + best_T(:,:,i);
            end
        end
    end

end  % orderings

pS     = best_pS;
T      = best_T;
W      = best_W;
status = best_st;
% =========================================================================
% PROTOCOL 3: QC-convFO  (Eq. 5.15)
% Source: channel_discrimination_qcconvfo_qcsupfo.m, case 7
%
% W is a convex mixture of 6 independent SEQ processes, one per sigma in S3.
% Each W^sigma (64x64, sigma-natural order) satisfies the same 3 SEQ comb
% conditions. Permutations: sigma-natural -> canonical I1O1I2O2I3O3.
% =========================================================================
case 3

cvx_clear
cvx_begin SDP
cvx_solver mosek

    % Defining W_sigma >= 0 
    % 64 x 64 PSD matrices, one per ordering
    variable W_123(dim,dim) hermitian semidefinite
    variable W_132(dim,dim) hermitian semidefinite
    variable W_213(dim,dim) hermitian semidefinite
    variable W_231(dim,dim) hermitian semidefinite
    variable W_312(dim,dim) hermitian semidefinite
    variable W_321(dim,dim) hermitian semidefinite

    % For each hypothesis n, a 64 x 64 PSD matrix T_123{n}. 
    % There are 6 x N matrices in total.
    T_123 = cell(N,1); for n=1:N; T_123{n} = semidefinite(dim,true); end
    T_132 = cell(N,1); for n=1:N; T_132{n} = semidefinite(dim,true); end
    T_213 = cell(N,1); for n=1:N; T_213{n} = semidefinite(dim,true); end
    T_231 = cell(N,1); for n=1:N; T_231{n} = semidefinite(dim,true); end
    T_312 = cell(N,1); for n=1:N; T_312{n} = semidefinite(dim,true); end
    T_321 = cell(N,1); for n=1:N; T_321{n} = semidefinite(dim,true); end

    expression S_123(dim,dim); S_123 = zeros(dim,dim);
    expression S_132(dim,dim); S_132 = zeros(dim,dim);
    expression S_213(dim,dim); S_213 = zeros(dim,dim);
    expression S_231(dim,dim); S_231 = zeros(dim,dim);
    expression S_312(dim,dim); S_312 = zeros(dim,dim);
    expression S_321(dim,dim); S_321 = zeros(dim,dim);
    expression pS_expr; pS_expr = 0;

    
    for n = 1:N
        % Sum of testers over hypotheses for each ordering
        S_123 = S_123 + T_123{n};
        S_132 = S_132 + T_132{n};
        S_213 = S_213 + T_213{n};
        S_231 = S_231 + T_231{n};
        S_312 = S_312 + T_312{n};
        S_321 = S_321 + T_321{n};
        expression Tn(dim,dim);
        Tn = T_123{n} ...
           + PermuteSystems(T_132{n}, [1 2 5 6 3 4], dim6) ...
           + PermuteSystems(T_213{n}, [3 4 1 2 5 6], dim6) ...
           + PermuteSystems(T_231{n}, [5 6 1 2 3 4], dim6) ...
           + PermuteSystems(T_312{n}, [3 4 5 6 1 2], dim6) ...
           + PermuteSystems(T_321{n}, [5 6 3 4 1 2], dim6);
        pS_expr = pS_expr + real(trace(p_i(n) * Tn * C3(:,:,n)));
    end

    % sum_n T_sigma{n} == W_sigma  (sigma-natural order)
    S_123 == W_123;  S_132 == W_132;  S_213 == W_213;
    S_231 == W_231;  S_312 == W_312;  S_321 == W_321;

    % SEQ comb conditions — identical for all 6 orderings
    % slots 1..6 refer to sigma-natural order of each W^sigma
    W_123 == kron(PartialTrace(W_123, 6, dim6), eye(d)/d);
    PartialTrace(W_123, [5,6],     dim6) == kron(PartialTrace(W_123, [4,5,6],     dim6), eye(d)/d);
    PartialTrace(W_123, [3,4,5,6], dim6) == kron(PartialTrace(W_123, [2,3,4,5,6], dim6), eye(d)/d);

    W_132 == kron(PartialTrace(W_132, 6, dim6), eye(d)/d);
    PartialTrace(W_132, [5,6],     dim6) == kron(PartialTrace(W_132, [4,5,6],     dim6), eye(d)/d);
    PartialTrace(W_132, [3,4,5,6], dim6) == kron(PartialTrace(W_132, [2,3,4,5,6], dim6), eye(d)/d);

    W_213 == kron(PartialTrace(W_213, 6, dim6), eye(d)/d);
    PartialTrace(W_213, [5,6],     dim6) == kron(PartialTrace(W_213, [4,5,6],     dim6), eye(d)/d);
    PartialTrace(W_213, [3,4,5,6], dim6) == kron(PartialTrace(W_213, [2,3,4,5,6], dim6), eye(d)/d);

    W_231 == kron(PartialTrace(W_231, 6, dim6), eye(d)/d);
    PartialTrace(W_231, [5,6],     dim6) == kron(PartialTrace(W_231, [4,5,6],     dim6), eye(d)/d);
    PartialTrace(W_231, [3,4,5,6], dim6) == kron(PartialTrace(W_231, [2,3,4,5,6], dim6), eye(d)/d);

    W_312 == kron(PartialTrace(W_312, 6, dim6), eye(d)/d);
    PartialTrace(W_312, [5,6],     dim6) == kron(PartialTrace(W_312, [4,5,6],     dim6), eye(d)/d);
    PartialTrace(W_312, [3,4,5,6], dim6) == kron(PartialTrace(W_312, [2,3,4,5,6], dim6), eye(d)/d);

    W_321 == kron(PartialTrace(W_321, 6, dim6), eye(d)/d);
    PartialTrace(W_321, [5,6],     dim6) == kron(PartialTrace(W_321, [4,5,6],     dim6), eye(d)/d);
    PartialTrace(W_321, [3,4,5,6], dim6) == kron(PartialTrace(W_321, [2,3,4,5,6], dim6), eye(d)/d);

    trace(W_123)+trace(W_132)+trace(W_213)+trace(W_231)+trace(W_312)+trace(W_321) == dOut^k;

    maximise real(pS_expr);
cvx_end
status = cvx_status;

T = zeros(dim,dim,N);
W = zeros(dim,dim);
for n = 1:N
    T(:,:,n) = T_123{n} ...
             + PermuteSystems(T_132{n}, [1 2 5 6 3 4], dim6) ...
             + PermuteSystems(T_213{n}, [3 4 1 2 5 6], dim6) ...
             + PermuteSystems(T_231{n}, [5 6 1 2 3 4], dim6) ...
             + PermuteSystems(T_312{n}, [3 4 5 6 1 2], dim6) ...
             + PermuteSystems(T_321{n}, [5 6 3 4 1 2], dim6);
    W = W + T(:,:,n);
end
pS = real(pS_expr);

% =========================================================================
% PROTOCOL 4: QC-NICC  (Eq. 5.17)
% Source: channel_discrimination_qcnicc_qcniqc.m, case 7
%
% = pQC-CC routing (C1)-(C5) + NI conditions on all routing variables.
% =========================================================================
case 4

w123_cell = cell(N_loc,1); w132_cell = cell(N_loc,1);
w213_cell = cell(N_loc,1); w231_cell = cell(N_loc,1);
w312_cell = cell(N_loc,1); w321_cell = cell(N_loc,1);

cvx_clear
cvx_begin SDP
cvx_solver mosek

    variable w1(d_loc,d_loc) hermitian semidefinite
    variable w2(d_loc,d_loc) hermitian semidefinite
    variable w3(d_loc,d_loc) hermitian semidefinite

    variable w12(d_loc^3,d_loc^3) hermitian semidefinite
    variable w13(d_loc^3,d_loc^3) hermitian semidefinite
    variable w21(d_loc^3,d_loc^3) hermitian semidefinite
    variable w23(d_loc^3,d_loc^3) hermitian semidefinite
    variable w31(d_loc^3,d_loc^3) hermitian semidefinite
    variable w32(d_loc^3,d_loc^3) hermitian semidefinite

    variable w123(d_loc^5,d_loc^5) hermitian semidefinite
    variable w132(d_loc^5,d_loc^5) hermitian semidefinite
    variable w213(d_loc^5,d_loc^5) hermitian semidefinite
    variable w231(d_loc^5,d_loc^5) hermitian semidefinite
    variable w312(d_loc^5,d_loc^5) hermitian semidefinite
    variable w321(d_loc^5,d_loc^5) hermitian semidefinite

    for ii = 1:N_loc
        w123_cell{ii} = semidefinite(d_loc^6,true);
        w132_cell{ii} = semidefinite(d_loc^6,true);
        w213_cell{ii} = semidefinite(d_loc^6,true);
        w231_cell{ii} = semidefinite(d_loc^6,true);
        w312_cell{ii} = semidefinite(d_loc^6,true);
        w321_cell{ii} = semidefinite(d_loc^6,true);
    end

    expression w123F(d_loc^6,d_loc^6); w123F = zeros(d_loc^6,d_loc^6);
    expression w132F(d_loc^6,d_loc^6); w132F = zeros(d_loc^6,d_loc^6);
    expression w213F(d_loc^6,d_loc^6); w213F = zeros(d_loc^6,d_loc^6);
    expression w231F(d_loc^6,d_loc^6); w231F = zeros(d_loc^6,d_loc^6);
    expression w312F(d_loc^6,d_loc^6); w312F = zeros(d_loc^6,d_loc^6);
    expression w321F(d_loc^6,d_loc^6); w321F = zeros(d_loc^6,d_loc^6);

    expression W_cvx(dim_loc,dim_loc);  W_cvx = zeros(dim_loc,dim_loc);

    expression pS; pS = 0;

    for ii = 1:N_loc
        w123F = w123F + w123_cell{ii};
        w132F = w132F + w132_cell{ii};
        w213F = w213F + w213_cell{ii};
        w231F = w231F + w231_cell{ii};
        w312F = w312F + w312_cell{ii};
        w321F = w321F + w321_cell{ii};

        expression Ti_slice(dim_loc,dim_loc);

        Ti_slice = w123_cell{ii}+w132_cell{ii}+w213_cell{ii}+w231_cell{ii}+w312_cell{ii}+w321_cell{ii};

        W_cvx = W_cvx + Ti_slice;
        
        pS = pS + real(trace(p_i(ii)*Ti_slice*C3(:,:,ii)));
    end

    % pQC-CC constraints as it is
    trace(w1)+trace(w2)+trace(w3) == 1;

    PartialTrace(w12,3,dim3)+PartialTrace(w13,3,dim3) == kron(w1,eye(d_loc));
    PartialTrace(w21,3,dim3)+PartialTrace(w23,3,dim3) == kron(w2,eye(d_loc));
    PartialTrace(w31,3,dim3)+PartialTrace(w32,3,dim3) == kron(w3,eye(d_loc));

    PartialTrace(w123,5,dim5) == kron(w12,eye(d_loc));
    PartialTrace(w132,5,dim5) == kron(w13,eye(d_loc));
    PartialTrace(w213,5,dim5) == kron(w21,eye(d_loc));
    PartialTrace(w231,5,dim5) == kron(w23,eye(d_loc));
    PartialTrace(w312,5,dim5) == kron(w31,eye(d_loc));
    PartialTrace(w321,5,dim5) == kron(w32,eye(d_loc));

    kron(w123,eye(d_loc))                                    == w123F;
    PermuteSystems(kron(w132,eye(d_loc)),[1 2 5 6 3 4],dim6) == w132F;
    PermuteSystems(kron(w213,eye(d_loc)),[3 4 1 2 5 6],dim6) == w213F;
    PermuteSystems(kron(w231,eye(d_loc)),[5 6 1 2 3 4],dim6) == w231F;
    PermuteSystems(kron(w312,eye(d_loc)),[3 4 5 6 1 2],dim6) == w312F;
    PermuteSystems(kron(w321,eye(d_loc)),[5 6 3 4 1 2],dim6) == w321F;

    real(trace(W_cvx)) == dOut^k;

    % NI-L1: [1-A^O_k1] Tr_{A^I_k2} W_(k1,k2) = 0
    PartialTrace(w12,3,dim3) == TR(PartialTrace(w12,3,dim3),2,[d_loc d_loc]);
    PartialTrace(w13,3,dim3) == TR(PartialTrace(w13,3,dim3),2,[d_loc d_loc]);
    PartialTrace(w21,3,dim3) == TR(PartialTrace(w21,3,dim3),2,[d_loc d_loc]);
    PartialTrace(w23,3,dim3) == TR(PartialTrace(w23,3,dim3),2,[d_loc d_loc]);
    PartialTrace(w31,3,dim3) == TR(PartialTrace(w31,3,dim3),2,[d_loc d_loc]);
    PartialTrace(w32,3,dim3) == TR(PartialTrace(w32,3,dim3),2,[d_loc d_loc]);

    % NI-L2: three conditions per level-2 variable
    % (a) K={k1}: trace pos{3,4,5}, apply [1-A^O_k1] at pos 2
    % (b) K={k2}: trace pos{1,2,5}, apply [1-A^O_k2] at pos 2
    % (c) K={k1,k2}: inclusion-exclusion at pos{2,4} after tracing pos 5
    PartialTrace(PartialTrace(PartialTrace(w123,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w123,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(PartialTrace(PartialTrace(w123,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w123,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(w123,5,dim5) - TR(PartialTrace(w123,5,dim5),2,dim4r) ...
        - TR(PartialTrace(w123,5,dim5),4,dim4r) ...
        + TR(TR(PartialTrace(w123,5,dim5),2,dim4r),4,dim4r) == zeros(d_loc^4,d_loc^4);

    PartialTrace(PartialTrace(PartialTrace(w132,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w132,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(PartialTrace(PartialTrace(w132,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w132,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(w132,5,dim5) - TR(PartialTrace(w132,5,dim5),2,dim4r) ...
        - TR(PartialTrace(w132,5,dim5),4,dim4r) ...
        + TR(TR(PartialTrace(w132,5,dim5),2,dim4r),4,dim4r) == zeros(d_loc^4,d_loc^4);

    PartialTrace(PartialTrace(PartialTrace(w213,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w213,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(PartialTrace(PartialTrace(w213,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w213,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(w213,5,dim5) - TR(PartialTrace(w213,5,dim5),2,dim4r) ...
        - TR(PartialTrace(w213,5,dim5),4,dim4r) ...
        + TR(TR(PartialTrace(w213,5,dim5),2,dim4r),4,dim4r) == zeros(d_loc^4,d_loc^4);

    PartialTrace(PartialTrace(PartialTrace(w231,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w231,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(PartialTrace(PartialTrace(w231,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w231,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(w231,5,dim5) - TR(PartialTrace(w231,5,dim5),2,dim4r) ...
        - TR(PartialTrace(w231,5,dim5),4,dim4r) ...
        + TR(TR(PartialTrace(w231,5,dim5),2,dim4r),4,dim4r) == zeros(d_loc^4,d_loc^4);

    PartialTrace(PartialTrace(PartialTrace(w312,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w312,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(PartialTrace(PartialTrace(w312,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w312,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(w312,5,dim5) - TR(PartialTrace(w312,5,dim5),2,dim4r) ...
        - TR(PartialTrace(w312,5,dim5),4,dim4r) ...
        + TR(TR(PartialTrace(w312,5,dim5),2,dim4r),4,dim4r) == zeros(d_loc^4,d_loc^4);

    PartialTrace(PartialTrace(PartialTrace(w321,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w321,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(PartialTrace(PartialTrace(w321,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(w321,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(w321,5,dim5) - TR(PartialTrace(w321,5,dim5),2,dim4r) ...
        - TR(PartialTrace(w321,5,dim5),4,dim4r) ...
        + TR(TR(PartialTrace(w321,5,dim5),2,dim4r),4,dim4r) == zeros(d_loc^4,d_loc^4);

    maximise real(pS);
cvx_end
status = cvx_status;

T = zeros(dim_loc,dim_loc,N_loc);
W = zeros(dim_loc,dim_loc);
for ii = 1:N_loc
    T(:,:,ii) = w123_cell{ii}+w132_cell{ii}+w213_cell{ii}+w231_cell{ii}+w312_cell{ii}+w321_cell{ii};
    W = W + T(:,:,ii);
end
pS = real(pS);

% =========================================================================
% PROTOCOL 5: pQC-CC
% Source: channel_discrimination_3copies_qccc_qcqc.m, case 5
% =========================================================================
case 5

w123_cell = cell(N_loc,1); w132_cell = cell(N_loc,1);
w213_cell = cell(N_loc,1); w231_cell = cell(N_loc,1);
w312_cell = cell(N_loc,1); w321_cell = cell(N_loc,1);

cvx_clear
cvx_begin SDP
cvx_solver mosek

    variable w1(d_loc,d_loc) hermitian semidefinite
    variable w2(d_loc,d_loc) hermitian semidefinite
    variable w3(d_loc,d_loc) hermitian semidefinite

    variable w12(d_loc^3,d_loc^3) hermitian semidefinite
    variable w13(d_loc^3,d_loc^3) hermitian semidefinite
    variable w21(d_loc^3,d_loc^3) hermitian semidefinite
    variable w23(d_loc^3,d_loc^3) hermitian semidefinite
    variable w31(d_loc^3,d_loc^3) hermitian semidefinite
    variable w32(d_loc^3,d_loc^3) hermitian semidefinite

    variable w123(d_loc^5,d_loc^5) hermitian semidefinite
    variable w132(d_loc^5,d_loc^5) hermitian semidefinite
    variable w213(d_loc^5,d_loc^5) hermitian semidefinite
    variable w231(d_loc^5,d_loc^5) hermitian semidefinite
    variable w312(d_loc^5,d_loc^5) hermitian semidefinite
    variable w321(d_loc^5,d_loc^5) hermitian semidefinite

    for ii = 1:N_loc
        w123_cell{ii} = semidefinite(d_loc^6,true);
        w132_cell{ii} = semidefinite(d_loc^6,true);
        w213_cell{ii} = semidefinite(d_loc^6,true);
        w231_cell{ii} = semidefinite(d_loc^6,true);
        w312_cell{ii} = semidefinite(d_loc^6,true);
        w321_cell{ii} = semidefinite(d_loc^6,true);
    end

    expression w123F(d_loc^6,d_loc^6); w123F = zeros(d_loc^6,d_loc^6);
    expression w132F(d_loc^6,d_loc^6); w132F = zeros(d_loc^6,d_loc^6);
    expression w213F(d_loc^6,d_loc^6); w213F = zeros(d_loc^6,d_loc^6);
    expression w231F(d_loc^6,d_loc^6); w231F = zeros(d_loc^6,d_loc^6);
    expression w312F(d_loc^6,d_loc^6); w312F = zeros(d_loc^6,d_loc^6);
    expression w321F(d_loc^6,d_loc^6); w321F = zeros(d_loc^6,d_loc^6);

    expression W_cvx(dim_loc,dim_loc);  W_cvx = zeros(dim_loc,dim_loc);

    expression pS; pS = 0;

    for ii = 1:N_loc
        w123F = w123F + w123_cell{ii};
        w132F = w132F + w132_cell{ii};
        w213F = w213F + w213_cell{ii};
        w231F = w231F + w231_cell{ii};
        w312F = w312F + w312_cell{ii};
        w321F = w321F + w321_cell{ii};

        expression Ti_slice(dim_loc,dim_loc);

        Ti_slice = w123_cell{ii}+w132_cell{ii}+w213_cell{ii}+w231_cell{ii}+w312_cell{ii}+w321_cell{ii};
        
        W_cvx = W_cvx + Ti_slice;

        pS = pS + real(trace(p_i(ii)*Ti_slice*C3(:,:,ii)));
    end

    trace(w1)+trace(w2)+trace(w3) == 1;

    PartialTrace(w12,3,dim3)+PartialTrace(w13,3,dim3) == kron(w1,eye(d_loc));
    PartialTrace(w21,3,dim3)+PartialTrace(w23,3,dim3) == kron(w2,eye(d_loc));
    PartialTrace(w31,3,dim3)+PartialTrace(w32,3,dim3) == kron(w3,eye(d_loc));

    PartialTrace(w123,5,dim5) == kron(w12,eye(d_loc));
    PartialTrace(w132,5,dim5) == kron(w13,eye(d_loc));
    PartialTrace(w213,5,dim5) == kron(w21,eye(d_loc));
    PartialTrace(w231,5,dim5) == kron(w23,eye(d_loc));
    PartialTrace(w312,5,dim5) == kron(w31,eye(d_loc));
    PartialTrace(w321,5,dim5) == kron(w32,eye(d_loc));

    kron(w123,eye(d_loc))                                    == w123F;
    PermuteSystems(kron(w132,eye(d_loc)),[1 2 5 6 3 4],dim6) == w132F;
    PermuteSystems(kron(w213,eye(d_loc)),[3 4 1 2 5 6],dim6) == w213F;
    PermuteSystems(kron(w231,eye(d_loc)),[5 6 1 2 3 4],dim6) == w231F;
    PermuteSystems(kron(w312,eye(d_loc)),[3 4 5 6 1 2],dim6) == w312F;
    PermuteSystems(kron(w321,eye(d_loc)),[5 6 3 4 1 2],dim6) == w321F;

    real(trace(W_cvx)) == dOut^k;

    maximise real(pS);
cvx_end
status = cvx_status;

T = zeros(dim_loc,dim_loc,N_loc);
W = zeros(dim_loc,dim_loc);
for ii = 1:N_loc
    T(:,:,ii) = w123_cell{ii}+w132_cell{ii}+w213_cell{ii}+w231_cell{ii}+w312_cell{ii}+w321_cell{ii};
    W = W + T(:,:,ii);
end
pS = real(pS);

% =========================================================================
% PROTOCOL 6: QC-SupFO  (Eq. 5.20)
%
% TESTER STRUCTURE: identical to QC-QC (protocol 8).
% A single variable T(dim,dim,N) complex semidefinite is declared.
% W = sum_n T(:,:,n) is then constrained to lie in the QC-SupFO set
% via the group routing variables. No per-hypothesis group variables.
%
% This is strictly more general than declaring separate per-hypothesis
% group testers, and guarantees QC-SupFO = QC-NIQC for N <= 3.
% =========================================================================
case 6
 
cvx_clear
cvx_begin SDP
cvx_solver mosek
 
    % ---- Testers: same as QC-QC ----
    variable T(dim,dim,N) complex semidefinite
 
    % ---- Level-0: 2x2 per first-party choice ----
    variable W23_1(d,d) hermitian semidefinite
    variable W32_1(d,d) hermitian semidefinite
    variable W13_2(d,d) hermitian semidefinite
    variable W31_2(d,d) hermitian semidefinite
    variable W12_3(d,d) hermitian semidefinite
    variable W21_3(d,d) hermitian semidefinite
 
    % ---- Level-1: 8x8 per ordered pair ----
    variable W3_12(d^3,d^3) hermitian semidefinite
    variable W2_13(d^3,d^3) hermitian semidefinite
    variable W3_21(d^3,d^3) hermitian semidefinite
    variable W1_23(d^3,d^3) hermitian semidefinite
    variable W2_31(d^3,d^3) hermitian semidefinite
    variable W1_32(d^3,d^3) hermitian semidefinite
 
    % ---- Level-2: 32x32 per GROUP (3, not 6) ----
    variable W0g_123(d^5,d^5) hermitian semidefinite  % group {1,2}, last=3
    variable W0g_132(d^5,d^5) hermitian semidefinite  % group {1,3}, last=2
    variable W0g_231(d^5,d^5) hermitian semidefinite  % group {2,3}, last=1
 
    % ---- Build W and objective from T directly ----
    expression pS_expr; pS_expr = 0;
    expression W_total(dim,dim); W_total = zeros(dim,dim);
    for n = 1:N
        W_total  = W_total + T(:,:,n);
        pS_expr  = pS_expr + real(trace(p_i(n) * T(:,:,n) * C3(:,:,n)));
    end
 
    % ---- Process matrix decomposition (Eq. 5.20, line 1) ----
    W_total == kron(W0g_123,eye(d)) ...
             + PermuteSystems(kron(W0g_132,eye(d)),[1 2 5 6 3 4],dim6) ...
             + PermuteSystems(kron(W0g_231,eye(d)),[5 6 1 2 3 4],dim6);
 
    % ---- Level-2 comb: SUM over both orderings in each group ----
    PartialTrace(W0g_123,5,dim5) == PermuteSystems(kron(W3_21,eye(d)),[3 4 1 2],dim4) ...
                                  + kron(W3_12,eye(d));
    PartialTrace(W0g_132,5,dim5) == PermuteSystems(kron(W2_31,eye(d)),[3 4 1 2],dim4) ...
                                  + kron(W2_13,eye(d));
    PartialTrace(W0g_231,5,dim5) == PermuteSystems(kron(W1_32,eye(d)),[3 4 1 2],dim4) ...
                                  + kron(W1_23,eye(d));
 
    % ---- Level-1 comb: individual per ordered pair ----
    PartialTrace(W3_12,3,dim3) == kron(W23_1,eye(d));
    PartialTrace(W2_13,3,dim3) == kron(W32_1,eye(d));
    PartialTrace(W3_21,3,dim3) == kron(W13_2,eye(d));
    PartialTrace(W1_23,3,dim3) == kron(W31_2,eye(d));
    PartialTrace(W2_31,3,dim3) == kron(W12_3,eye(d));
    PartialTrace(W1_32,3,dim3) == kron(W21_3,eye(d));
 
    % ---- Normalisation ----
    trace(W23_1)+trace(W32_1)+trace(W13_2)+trace(W31_2)+trace(W12_3)+trace(W21_3) == 1;
    real(trace(W_total)) == dOut^k;
 
    maximise real(pS_expr);
cvx_end
status = cvx_status;
 
W = zeros(dim,dim);
for n = 1:N
    W = W + T(:,:,n);
end
pS = real(pS_expr);
status = cvx_status;

% =========================================================================
% PROTOCOL 7: QC-NIQC  (Eq. 5.22)
% Source: channel_discrimination_qcnicc_qcniqc.m, case 8
%
% = QC-QC routing constraints + NI conditions on all routing variables.
% =========================================================================
case 7

cvx_clear
cvx_begin SDP
cvx_solver mosek

    variable wq1(d_loc,d_loc) hermitian semidefinite
    variable wq2(d_loc,d_loc) hermitian semidefinite
    variable wq3(d_loc,d_loc) hermitian semidefinite

    variable wq12(d_loc^3,d_loc^3) hermitian semidefinite
    variable wq13(d_loc^3,d_loc^3) hermitian semidefinite
    variable wq21(d_loc^3,d_loc^3) hermitian semidefinite
    variable wq23(d_loc^3,d_loc^3) hermitian semidefinite
    variable wq31(d_loc^3,d_loc^3) hermitian semidefinite
    variable wq32(d_loc^3,d_loc^3) hermitian semidefinite

    variable wq123(d_loc^5,d_loc^5) hermitian semidefinite
    variable wq132(d_loc^5,d_loc^5) hermitian semidefinite
    variable wq231(d_loc^5,d_loc^5) hermitian semidefinite

    variable T(dim_loc,dim_loc,N_loc) complex semidefinite

    expression pS; pS = 0;
    expression W_cvx(dim_loc,dim_loc); W_cvx = zeros(dim_loc,dim_loc);
    for ii = 1:N_loc
        pS    = pS    + real(trace(p_i(ii)*T(:,:,ii)*C3(:,:,ii)));
        W_cvx = W_cvx + T(:,:,ii);
    end

    % QC-QC conditions as it is
    trace(wq1)+trace(wq2)+trace(wq3) == 1;

    PartialTrace(wq12,3,dim3)+PartialTrace(wq13,3,dim3) == kron(wq1,eye(d_loc));
    PartialTrace(wq21,3,dim3)+PartialTrace(wq23,3,dim3) == kron(wq2,eye(d_loc));
    PartialTrace(wq31,3,dim3)+PartialTrace(wq32,3,dim3) == kron(wq3,eye(d_loc));

    PartialTrace(wq123,5,dim5) == PermuteSystems(kron(wq21,eye(d_loc)),[3 4 1 2],dim4)+kron(wq12,eye(d_loc));
    PartialTrace(wq132,5,dim5) == PermuteSystems(kron(wq31,eye(d_loc)),[3 4 1 2],dim4)+kron(wq13,eye(d_loc));
    PartialTrace(wq231,5,dim5) == PermuteSystems(kron(wq32,eye(d_loc)),[3 4 1 2],dim4)+kron(wq23,eye(d_loc));
    
    W_cvx == kron(wq123,eye(d_loc)) ...
           + PermuteSystems(kron(wq132,eye(d_loc)),[1 2 5 6 3 4],dim6) ...
           + PermuteSystems(kron(wq231,eye(d_loc)),[5 6 1 2 3 4],dim6);

    real(trace(W_cvx)) == dOut^k;
    
    % ∏​[1−AiO​]\ TrAN∖XIO ​​M = 0 [Validity condition for intermediate
    % matrices]

    % NI-L1
    % [1−Ak1​O​]\ TrAk2​I ​​wq({k1​},k2​) ​= 0 (X = {k1})
    % \TrAk2​I​​ wq({k1​},k2​) ​= TR( \TrAk2​I ​​wq({k1​},k2​) ​,2,[d,d])
    PartialTrace(wq12,3,dim3) == TR(PartialTrace(wq12,3,dim3),2,[d_loc d_loc]);
    PartialTrace(wq13,3,dim3) == TR(PartialTrace(wq13,3,dim3),2,[d_loc d_loc]);
    PartialTrace(wq21,3,dim3) == TR(PartialTrace(wq21,3,dim3),2,[d_loc d_loc]);
    PartialTrace(wq23,3,dim3) == TR(PartialTrace(wq23,3,dim3),2,[d_loc d_loc]);
    PartialTrace(wq31,3,dim3) == TR(PartialTrace(wq31,3,dim3),2,[d_loc d_loc]);
    PartialTrace(wq32,3,dim3) == TR(PartialTrace(wq32,3,dim3),2,[d_loc d_loc]);

    % NI-L2 (3 group variables)
    PartialTrace(PartialTrace(PartialTrace(wq123,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(wq123,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(PartialTrace(PartialTrace(wq123,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(wq123,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(wq123,5,dim5) - TR(PartialTrace(wq123,5,dim5),2,dim4r) ...
        - TR(PartialTrace(wq123,5,dim5),4,dim4r) ...
        + TR(TR(PartialTrace(wq123,5,dim5),2,dim4r),4,dim4r) == zeros(d_loc^4,d_loc^4);

    PartialTrace(PartialTrace(PartialTrace(wq132,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(wq132,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(PartialTrace(PartialTrace(wq132,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(wq132,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(wq132,5,dim5) - TR(PartialTrace(wq132,5,dim5),2,dim4r) ...
        - TR(PartialTrace(wq132,5,dim5),4,dim4r) ...
        + TR(TR(PartialTrace(wq132,5,dim5),2,dim4r),4,dim4r) == zeros(d_loc^4,d_loc^4);

    PartialTrace(PartialTrace(PartialTrace(wq231,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(wq231,5,dim5),4,dim4r),3,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(PartialTrace(PartialTrace(wq231,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]) == ...
        TR(PartialTrace(PartialTrace(PartialTrace(wq231,5,dim5),2,dim4r),1,[d_loc d_loc d_loc]),2,[d_loc d_loc]);
    PartialTrace(wq231,5,dim5) - TR(PartialTrace(wq231,5,dim5),2,dim4r) ...
        - TR(PartialTrace(wq231,5,dim5),4,dim4r) ...
        + TR(TR(PartialTrace(wq231,5,dim5),2,dim4r),4,dim4r) == zeros(d_loc^4,d_loc^4);

    maximise real(pS);
cvx_end
status = cvx_status;


W = zeros(dim_loc,dim_loc);
for ii = 1:N_loc; W = W + T(:,:,ii); end
pS = real(pS);

% =========================================================================
% PROTOCOL 8: QC-QC
% Source: channel_discrimination_3copies_qccc_qcqc.m, case 6
% =========================================================================
case 8

cvx_clear
cvx_begin SDP
cvx_solver mosek

    variable T(dim,dim,N) complex semidefinite

    variable wq1(d,d) hermitian semidefinite
    variable wq2(d,d) hermitian semidefinite
    variable wq3(d,d) hermitian semidefinite

    variable wq12(d^3,d^3) hermitian semidefinite
    variable wq13(d^3,d^3) hermitian semidefinite
    variable wq21(d^3,d^3) hermitian semidefinite
    variable wq23(d^3,d^3) hermitian semidefinite
    variable wq31(d^3,d^3) hermitian semidefinite
    variable wq32(d^3,d^3) hermitian semidefinite

    variable wq123(d^5,d^5) hermitian semidefinite
    variable wq132(d^5,d^5) hermitian semidefinite
    variable wq231(d^5,d^5) hermitian semidefinite

    expression pS; pS = 0;
    expression W;  W  = 0;
    for i = 1:N
        pS = pS + real(trace(p_i(i)*T(:,:,i)*C3(:,:,i)));
        W  = W  + T(:,:,i);
    end
    
    % wq1 lives on subsystem A^I_k1
    % ∑k1​​\ TrAk1​I​​ W_(∅,k1​)​ = 1
    trace(wq1)+trace(wq2)+trace(wq3) == 1;
    
    % ∑k2​​\ TrAk2​I​ ​W({1},k2​) ​= W(∅,1)​ ⊗ \id^A1O​
    PartialTrace(wq12,3,dim3)+PartialTrace(wq13,3,dim3) == kron(wq1,eye(d));
    PartialTrace(wq21,3,dim3)+PartialTrace(wq23,3,dim3) == kron(wq2,eye(d));
    PartialTrace(wq31,3,dim3)+PartialTrace(wq32,3,dim3) == kron(wq3,eye(d));
    
    % \TrAk3​I ​​W({k1​,k2​},k3​) ​= W({k2​},k1​) ​⊗ \id + W({k1​},k2​)​ ⊗ \id
    PartialTrace(wq123,5,dim5) == PermuteSystems(kron(wq21,eye(d)),[3 4 1 2],dim4)+kron(wq12,eye(d));
    PartialTrace(wq132,5,dim5) == PermuteSystems(kron(wq31,eye(d)),[3 4 1 2],dim4)+kron(wq13,eye(d));
    PartialTrace(wq231,5,dim5) == PermuteSystems(kron(wq32,eye(d)),[3 4 1 2],dim4)+kron(wq23,eye(d));
    
    % W=∑kN ​​W(N∖kN​,kN​)​ ⊗ \idAkN​O​
    W == PermuteSystems(kron(wq231,eye(d)),[5 6 1 2 3 4],dim6) ...
       + PermuteSystems(kron(wq132,eye(d)),[1 2 5 6 3 4],dim6) ...
       + kron(wq123,eye(d));

    real(trace(W)) == dOut^k;

    maximise real(pS);
cvx_end
status = cvx_status;

% =========================================================================
% PROTOCOL 9: GEN
% Source: channel_discrimination_3copies_primal.m, case 4
% =========================================================================
case 9

cvx_clear
cvx_begin SDP
cvx_solver mosek
    variable T(dim,dim,N) complex semidefinite
    expression pS; pS = 0;
    expression W;  W  = 0;
    for i = 1:N
        pS = pS + real(trace(p_i(i)*T(:,:,i)*C3(:,:,i)));
        W  = W  + T(:,:,i);
    end
    % (Araújo) Witness Causal Nonseparability [Section II]
    % ∏​[1−AkO​]\TrA{1,2,3}∖XIO ​​W = 0 ∀∅ != X⊆{1,2,3}
    
    % [1−A3O​]\ TrA1IO​A2IO​​ W = 0 (X = {3})
    TR(W,[1 2 3 4],DIM) == TR(W,[1 2 3 4 6],DIM);

    % [1−A1O​]\ TrA2IO​A3IO ​​W = 0 (X = {1})
    TR(W,[3 4 5 6],DIM) == TR(W,[2 3 4 5 6],DIM);

    % [1−A2O​]\ TrA1IO​A3IO​​ W = 0 (X = {2})
    TR(W,[1 2 5 6],DIM) == TR(W,[1 2 4 5 6],DIM);
    
    % [1−A2O​][1−A3O​]\ TrA1IO​​ W = 0 (X = {2,3})
    TR(W,[1 2],DIM)+TR(W,[1 2 4 6],DIM) == TR(W,[1 2 4],DIM)+TR(W,[1 2 6],DIM);

    % [1−A1O​][1−A3O​]\ TrA2IO ​​W = 0 (X = {1,3})
    TR(W,[3 4],DIM)+TR(W,[2 3 4 6],DIM) == TR(W,[2 3 4],DIM)+TR(W,[3 4 6],DIM);

    % [1−A1O​][1−A2O​]\ TrA3IO ​​W = 0 (X = {1,2})
    TR(W,[5 6],DIM)+TR(W,[2 4 5 6],DIM) == TR(W,[4 5 6],DIM)+TR(W,[2 5 6],DIM);
    
    % [1−A1O​][1−A2O​][1−A3O​] W = 0 (X = {1, 2, 3})
    W == TR(W,[2 4 6],DIM)+TR(W,2,DIM)+TR(W,4,DIM)+TR(W,6,DIM) ...
       - TR(W,[2 4],DIM)-TR(W,[2 6],DIM)-TR(W,[4 6],DIM);

    trace(W) == dOut^k;

    maximise real(pS);
cvx_end
status = cvx_status;

% =========================================================================
otherwise
    error('Protocol must be 1..9. See function header for the full list.');
end

end