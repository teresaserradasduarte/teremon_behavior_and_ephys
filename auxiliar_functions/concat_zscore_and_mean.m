function FR_mean = concat_zscore_and_mean(FR_all)

nr_bins = size(cell2mat(FR_all(1)),1);
FR_mean = nan(nr_bins,length(FR_all));

for i = 1:length(FR_all) 
    data = cell2mat(FR_all(i));

    data_flat = data(:);
    zs_data_flat = zscore(data_flat);
    zs_data_rs = reshape(zs_data_flat, size(data));

    mean_zs_data = mean(zs_data_rs,2,"omitnan");

    FR_mean(:,i) = mean_zs_data;
end

end

