clear all; clc;

% load data file (csv) without header (no column names in the file)
data = readtable('protein seq.csv', 'ReadVariableNames', 0);

% Init columns for output_data table:
Var1 = {};
Var2 = {};
First_D_E_position = {};
uH = {};
R_count = {};
H_count = {};
K_count = {};

for i = 1:size(data, 1)
    curSeq = data.Var3{i}(1:end-1); % remove stop codon ('*') at the end of the sequence
	Var1 = [Var1; data.Var1{i}];
	Var2 = [Var2; data.Var2{i}];
    
    D_E_1st_pos = find((curSeq == 'E') | (curSeq == 'D'),1);
    if ~isempty(D_E_1st_pos) % if there is an E or D in 'curSeq'
        First_D_E_position = [First_D_E_position; D_E_1st_pos];
    else % if there is no D or E in 'curSeq'
        First_D_E_position = [First_D_E_position; {' '}];
    end
    
    if ~isempty(D_E_1st_pos) % if there is an E or D in 'curSeq'
        if D_E_1st_pos <= 30
            uH = [uH; uH_mean_amphipathic_moment(curSeq(1:D_E_1st_pos - 1))];
            AAcount = aacount(curSeq(1:D_E_1st_pos - 1));
        else
            uH = [uH; uH_mean_amphipathic_moment(curSeq(1:30))];
            AAcount = aacount(curSeq(1:30));
        end
    else  % if there is no E or D in 'curSeq'
        if length(curSeq) >= 30
            uH = [uH; uH_mean_amphipathic_moment(curSeq(1:30))];
            AAcount = aacount(curSeq(1:30));
        else 
            uH = [uH; uH_mean_amphipathic_moment(curSeq)];
            AAcount = aacount(curSeq);
        end
    end

    R_count = [R_count; AAcount.R];
    H_count = [H_count; AAcount.H];
    K_count = [K_count; AAcount.K];
end

output_data = table(Var1, Var2, First_D_E_position, uH, R_count, H_count, K_count);

writetable(output_data,'Mean Amphipathic Moment.csv');





