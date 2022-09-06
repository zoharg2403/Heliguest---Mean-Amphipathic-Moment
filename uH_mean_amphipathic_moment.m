function [uH] = uH_mean_amphipathic_moment (seq)
%% Calculation of the mean amphipathic moment <uH>:
% Calculation of the mean amphipathic moment according to the formula:
% uH = 1/N * SQRT((sum(Hn*sin(nd)))^2 + (sum(Hn*cos(nd)))^2)
%       * N = The sequence length
%       * Hn = The hydrophobicity of the n-th amino acid in the sequence
%       * d = delta =  the angle separating side chains along the backbone 
%                      d = 100 for an alpha helix

% The µH value ranges from to 0 to 3.26

% Hn_values:
%       loaded from 'Hn_values.mat':       
%           Ala(A):0.310	Arg(R):-1.010	Asn(N):-0.600   Asp(D):-0.770	
%           Cys(C):1.540    Gln(Q):-0.220   Glu(E):-0.640   Gly(G):0.000	
%           His(H):0.130    Ile(I):1.800    Leu(L):1.700	Lys(K):-0.990  
%           Met(M):1.230	Phe(F):1.790	Pro(P):0.720    Ser(S):-0.040
%           Thr(T):0.260	Trp(W):2.250    Tyr(Y):0.960    Val(V):1.220

% Input:
%       * seq = Amino acid sequence (protein) -> string
%               The amino acids in the sequence shuld be represented by the
%               Code (A, R, N...) ('aminolookup' matlab function can be used for conversio)
% Output:
%       * uH = Mean amphipathic moment

%% 

d = 100 * pi / 180; % delta = 100 degree for an alpha helix (*pi/180 -> convert to rad)
Hn_values = load('Hn_values_2.mat'); % load Hn_values
Hn_values = Hn_values.Hn_values;

AAlist = fieldnames(Hn_values); % list of AA in 'AAcount' struct 

uH_sin = 0; % init sin part of the formula
uH_cos = 0; % init cos part of the formula

for a = 1:length(AAlist) % for each AA in AAlist
    curAA = AAlist{a};
    curAA_indices = find(seq == curAA); % get list of indices that 'curAA' appears
    uH_sin = uH_sin + Hn_values.(curAA) * sum(sin(curAA_indices .* d)); % sin part of the formula
    uH_cos = uH_cos + Hn_values.(curAA) * sum(cos(curAA_indices .* d)); % cos part of the formula
end

uH = sqrt(uH_sin^2 + uH_cos^2) / length(seq);

end

