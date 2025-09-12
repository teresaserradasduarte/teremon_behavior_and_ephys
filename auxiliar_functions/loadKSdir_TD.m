function spikeStruct = loadKSdir_TD(ksDir, varargin)

if ~isempty(varargin)
    params = varargin{1};
else
    params = [];
end

if ~isfield(params, 'excludeNoise')
    params.excludeNoise = true;
end

if exist(fullfile(ksDir,'pc_features.npy'),'file')==0
    params.loadPCs = false;
else
    params.loadPCs = true;
end

% load spike data
% struct
spikeStruct = loadParamsPy(fullfile(ksDir, 'params.py'));

% spike times
ss = readNPY(fullfile(ksDir, 'spike_times.npy'));
st = double(ss)/spikeStruct.sample_rate;

% templates idx
spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); %0idx

% unwhiten templates waveforms
templateWaveforms_whitened = readNPY(fullfile(ksDir, 'templates.npy'));
winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));
templateWaveforms = zeros(size(templateWaveforms_whitened));
for t = 1:size(templateWaveforms,1)
    templateWaveforms(t,:,:) = squeeze(templateWaveforms_whitened(t,:,:))*winv;
end
template_phyIDs = 0:size(templateWaveforms,1)-1;
[~,templatePeakChn] = max(abs(max(templateWaveforms(:,:,:),[],2)-min(templateWaveforms(:,:,:),[],2)),[],3);

% clusters
if exist(fullfile(ksDir, 'spike_clusters.npy'))
    clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
else
    clu = spikeTemplates;
end

% amplitudes
tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));

% PC features
if params.loadPCs
    pcFeat = readNPY(fullfile(ksDir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
    pcFeatInd = readNPY(fullfile(ksDir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
else
    pcFeat = [];
    pcFeatInd = [];
end

% Channels positions
coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
ycoords = coords(:,2); xcoords = coords(:,1);

% spike positions
if  exist(fullfile(ksDir, 'spike_positions.npy'),"file")
    spikePositions = readNPY(fullfile(ksDir, 'spike_positions.npy')); %0idx
else
    spikePositions = [];
end

if ~isempty(spikePositions)
        spikeDepths = spikePositions(:,2);

elseif ~isempty(pcFeat)
    pcFeat = squeeze(pcFeat(:,1,:)); % take first PC only
    pcFeat(pcFeat<0) = 0; % some entries are negative, but we don't really want to push the CoM away from there.
    % which channels for each spike?
    spikeFeatInd = pcFeatInd(spikeTemplates+1,:);
    

    % ycoords of those channels?
    spikeFeatYcoords = ycoords(spikeFeatInd+1); % 2D matrix of size #spikes x 12

    % center of mass is sum(coords.*features)/sum(features)
    spikeDepths = sum(spikeFeatYcoords.*pcFeat.^2,2)./sum(pcFeat.^2,2);
else
    spikeDepths= [];
end

% Quality of clusters: exclude noise
cgsFile = '';
cExtraGsFile = '';
if exist(fullfile(ksDir, 'cluster_groups.csv'),"file") 
    cgsFile = fullfile(ksDir, 'cluster_groups.csv');
end
if exist(fullfile(ksDir, 'cluster_group.tsv'),"file") 
   cgsFile = fullfile(ksDir, 'cluster_group.tsv');
end 
if exist(fullfile(ksDir, 'cluster_Extragood.tsv'),"file") 
   cExtraGsFile = fullfile(ksDir, 'cluster_Extragood.tsv');
elseif exist(fullfile(ksDir, 'cluster_extragood.tsv'),"file") 
   cExtraGsFile = fullfile(ksDir, 'cluster_extragood.tsv');
else
    cegs = zeros(size(unique(spikeTemplates)));
end 

if exist(fullfile(ksDir, 'cluster_CS.tsv'),"file")
   potential_CSFile = fullfile(ksDir, 'cluster_CS.tsv');
elseif exist(fullfile(ksDir, 'cluster_potential_CS.tsv'),"file")
    potential_CSFile = fullfile(ksDir, 'cluster_potential_CS.tsv');
else
    potential_CSFile = [];
end

if ~isempty(cgsFile)
    [cids, cgs] = readClusterGroupsCSV(cgsFile);
    cids_wnoise = cids;

    if params.excludeNoise
        noiseClusters = cids(cgs==0);

        st = st(~ismember(clu, noiseClusters));
        spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
        tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters)); 
        
        if params.loadPCs
            pcFeat = pcFeat(~ismember(clu, noiseClusters), :,:);
            %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
        end

        if ~isempty(spikePositions)
            spikePositions = spikePositions(~ismember(clu, noiseClusters),:);
        end
        if ~isempty(spikeDepths)
            spikeDepths = spikeDepths(~ismember(clu, noiseClusters));
        end 
        
        clu = clu(~ismember(clu, noiseClusters));
        cgs = cgs(~ismember(cids, noiseClusters));
        cids = cids(~ismember(cids, noiseClusters));

        % extra good
        if ~isempty(cExtraGsFile)
            cegs = readClusterExtraGood(cExtraGsFile,cids);
        end 
        
        % potential complex spikes
        if ~isempty(potential_CSFile)
            pCS = readPotentialComplexSpikes(potential_CSFile,cids);
        else
            pCS = [];
        end
    end
else
    clu = spikeTemplates;
    cids = unique(spikeTemplates);
    cgs = 3*ones(size(cids));
    cegs = zeros(size(cids));
end


% save in strucutre
spikeStruct.st = st;
spikeStruct.spikeTemplates = spikeTemplates;
spikeStruct.clu = clu;
spikeStruct.tempScalingAmps = tempScalingAmps;
spikeStruct.cgs = cgs;
spikeStruct.cegs = cegs;
spikeStruct.cids = cids;
spikeStruct.xcoords = xcoords;
spikeStruct.ycoords = ycoords;
spikeStruct.templateWaveforms = templateWaveforms;
spikeStruct.template_phyIDs = template_phyIDs;
spikeStruct.templatePeakChn = templatePeakChn;
spikeStruct.winv = winv;
spikeStruct.pcFeat = pcFeat;
spikeStruct.pcFeatInd = pcFeatInd;
spikeStruct.cids_wnoise = cids_wnoise;
spikeStruct.spikeDepths=spikeDepths;
spikeStruct.spikePositions=spikePositions;
spikeStruct.potential_CS=pCS;




