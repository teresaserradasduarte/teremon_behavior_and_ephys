function cEgs = readClusterExtraGood(filename,cids)
%function cEgs = readClusterExtraGood(filename,,n_cgs)
% input: ksfilename, cids (cluster ids excluding noise excluding noise)
% output: is cluster extra good
% - 0 = good / mua / unsorted
% - 1 = extragood


fid = fopen(filename);
C = textscan(fid, '%s%s');
fclose(fid);

cids_eg = cellfun(@str2num, C{1}(2:end), 'uni', false);
isExtraGood = cellfun(@(x)strcmp(x,'y'),C{2}(2:end));
extraGoodClusters = cell2mat(cids_eg(isExtraGood))';

cEgs = zeros(size(cids));

for i = 1:length(cids)
    cEgs(i) =  sum(extraGoodClusters==cids(i));
end
    
