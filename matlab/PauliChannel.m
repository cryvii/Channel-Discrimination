function C = PauliChannel(varargin)
% PauliChannel  Choi matrix of a Pauli channel on a single qubit.
%
% USAGE:
%   C = PauliChannel()                        % fully random
%   C = PauliChannel('seed', 42)              % seeded random
%   C = PauliChannel('theta', 0.8, 'q', [0.5 0.3 0.2])  % explicit
%   C = PauliChannel('seed', 'random')        % shuffle RNG

% --- Parse inputs ---
p = inputParser();
p.addParameter('theta', [],  @(x) isempty(x) || (isscalar(x) && x>=0 && x<=1));
p.addParameter('q',     [],  @(x) isempty(x) || (numel(x)==3 && abs(sum(x)-1)<1e-10));
p.addParameter('seed',  []);
p.parse(varargin{:});
cfg = p.Results;

% --- Set RNG ---
if ~isempty(cfg.seed)
    if ischar(cfg.seed) && strcmpi(cfg.seed, 'random')
        rng('shuffle');
    else
        rng(cfg.seed);
    end
end

% --- Draw parameters ---
theta = cfg.theta;
if isempty(theta)
    theta = rand();
end

q = cfg.q;
if isempty(q)
    e = -log(rand(1,3));
    q = e / sum(e);
else
    q = q(:)';
end

% --- Pauli matrices ---
I = eye(2);
X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1  0; 0 -1];
Paulis = {X, Y, Z};

% --- Kraus operators ---
K = cell(1,4);
K{1} = sqrt(theta) * I;
for i = 1:3
    K{i+1} = sqrt((1-theta) * q(i)) * Paulis{i};
end

% --- Choi matrix ---
C = ChoiMatrix(K(:));

% --- Print summary ---
TP_check = zeros(2,2);
for i = 1:4
    TP_check = TP_check + K{i}' * K{i};
end
fprintf('\n=== PauliChannel ===\n');
fprintf('  theta        : %.6f\n', theta);
fprintf('  q (X,Y,Z)   : [%.4f  %.4f  %.4f]  (sum=%.6f)\n', q(1),q(2),q(3),sum(q));
fprintf('  TP check     : max|sum K_i^dag K_i - I| = %.2e\n', max(max(abs(TP_check - I))));
fprintf('  Choi trace   : %.6f  (expected 2)\n', real(trace(C)));
fprintf('  Choi min eig : %.6f  (should be >= 0)\n', min(real(eig(C))));
fprintf('====================\n\n');
end
