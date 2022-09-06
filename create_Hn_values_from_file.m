

% import 'Eisenberg consensus normalized.xlsx' file as cell array and rename to data

for n = 1:20
    AA = data{n,1};
    val = data{n,2};
    Hn_values.(AA) = val;
end