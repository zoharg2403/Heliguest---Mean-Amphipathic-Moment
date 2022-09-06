function [H] = H_mean_hydrophobicity(seq)
%% Calculation of the mean hydrophobicity <H>:
% Calculation of the mean hydrophobicity according to the formula:
% H = (1/N) * sum(H_n) 
%       * N = The sequence length
%       * Hn = The hydrophobicity of the n-th amino acid in the sequence
% The H value ranges from -1.01 to 2.25

% for shorter run time, the function that will actually be used is:
% H = 1/N * sum(#AA_n* Hn)  -->  for each AA, count how many time it
% appears, multiple by Hn_value for the current AA and sum

% Input:
%       * seq = amino acid sequence (protein) -> string
%               The amino acids in the sequence shuld be represented by the
%               Code (A, R, N...) ('aminolookup' matlab function can be used for conversio)
% Output:
%       * H = Mean hydrophobicity

% Hn_values:
%       loaded from 'Hn_values.mat':       
%           Ala(A):0.310	Arg(R):-1.010	Asn(N):-0.600   Asp(D):-0.770	
%           Cys(C):1.540    Gln(Q):-0.220   Glu(E):-0.640   Gly(G):0.000	
%           His(H):0.130    Ile(I):1.800    Leu(L):1.700	Lys(K):-0.990  
%           Met(M):1.230	Phe(F):1.790	Pro(P):0.720    Ser(S):-0.040
%           Thr(T):0.260	Trp(W):2.250    Tyr(Y):0.960    Val(V):1.220

%%

Hn_values = load('Hn_values.mat'); % load Hn_values
Hn_values = Hn_values.Hn_values;

AAcount = aacount(seq); % count how many times each AA appears
AAlist = fieldnames(AAcount); % list of AA in 'AAcount' struct 

H = 0; % init H value
for n = 1:length(AAlist) % for each AA in AAlist
    curAA = AAlist{n};
    H = H + AAcount.(curAA) * Hn_values.(curAA); % Hn for the n-th amino acid
end

H = H / length(seq); % divide the sum by N ('seq' length)

end

