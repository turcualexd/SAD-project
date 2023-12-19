function import_IGRF2020(filename)

m_ind = @(i) i+1;
g = zeros(13, m_ind(13));
h = zeros(13, m_ind(13));
table = readtable(filename);
for i = 1:size(table,1)
    switch char(table{i,1})
        case 'g'
            g(table{i,2}, m_ind(table{i,3})) = table{i,4};
        case 'h'
            h(table{i,2}, m_ind(table{i,3})) = table{i,4};
    end
end

save("IGRF13coeffs.mat", "g", "h")