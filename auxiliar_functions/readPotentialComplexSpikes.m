function pCS = readPotentialComplexSpikes(filename,cids)
%function pCS = readClusterCS(filename,n_cgs)
% input: ksfilename, cids (cluster ids excluding noise excluding noise)
% output: is cluster extra good
% - 0 = other cell type
% - 1 = potential complex spike


fid = fopen(filename);
C = textscan(fid, '%s%s');
fclose(fid);

cids_cs = cellfun(@str2num, C{1}(2:end), 'uni', false);
isCS = cellfun(@(x)strcmp(x,'y'),C{2}(2:end));
potentialCS = cell2mat(cids_cs(isCS))';

pCS = zeros(size(cids));

for i = 1:length(cids)
    pCS(i) =  sum(potentialCS==cids(i));
end
    
