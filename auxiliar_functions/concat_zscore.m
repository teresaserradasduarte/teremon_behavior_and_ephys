function zFR = concat_zscore(FR_all)

zFR = cell(size(FR_all));

for i = 1:length(FR_all) 
    data = FR_all{i};

    data_flat = data(:);
    zs_data_flat = zscore(data_flat);
    zs_data_rs = reshape(zs_data_flat, size(data));

    zFR{i} = zs_data_rs;
end

end

