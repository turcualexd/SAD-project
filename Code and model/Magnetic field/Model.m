clear, clc, close all;

%% Import IGRF 2020 coefficients from file

m_ind = @(i) i+1;
g = zeros(13, m_ind(13));
h = zeros(13, m_ind(13));
table = readtable('IGRF13coeffs.xls');
for i = 1:size(table,1)
    switch char(table{i,1})
        case 'g'
            g(table{i,2}, m_ind(table{i,3})) = table{i,4};
        case 'h'
            h(table{i,2}, m_ind(table{i,3})) = table{i,4};
    end
end

%% Summatory to obtain V

r = 21232;
phi = 0;
theta = 0;
k = 8;      % order of the magnetic field

sum_on_n = 0;
for n = 1:k

    sum_on_m = 0;
    for m = m_ind(0:n)
        sum_on_m = sum_on_m + (g(n,m)*cos(m*phi) + h(n,m)*sin(m*phi)) * 235235;
        % da sistemare
    end


end