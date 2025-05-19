%% Geenrate figures - MI landscapes
close all

%% Select
N_tmp = struct2table(N_BG_all);
reg = 'BG';
% exclude the only session with just neurons from the cerebellum
exclude_idx = strcmp(N_tmp.mouse, '1_CoteDor') & strcmp(N_tmp.sess, 'R1');
N = N_tmp(~exclude_idx, :);

N_small = N(:,{'mouse','sess','phyID'});
clear N_tmp N N_BG_all

idx_milka = find(strcmp(N_small.mouse, '4_Milka'));

MI_landscaps_BG = MI_landscaps_all;
%%
% Figure params
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};

nn = [415, 417, 421, 460];
save_CB_neu = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_and_behavior\out_files\20230801_ChocolateGroup\group\BG\MI2';
n_neu = size(MI_landscaps_BG,3);
figure
ff = tiledlayout(1,4)
for n=nn
    nexttile
    imagesc(x_im,y_im,MI_landscaps_BG(:,:,n));
    axis xy;        % Invert the y-axis direction
    hold on
    %plot(stop_n,start_n,'ko','LineWidth',.5)
    %imagesc(x_im,y_im,flipud(MI_paw_neu_smooth));
    title({
        'Mutual information: I(paw endpoint;neuron)';
        sprintf('region: %s | mouse: %s | session: %s | neuron phyID: %d', ...
        reg, ...
        char(N_small(n,:).mouse), ...
        char(N_small(n,:).sess), ...
        N_small(n,:).phyID)
        }, ...
        'Interpreter', 'none');

    alpha_data = ~isnan(MI_landscaps_BG(:,:,n));
    set(gca, 'ALim', [0 1]);
    hImg = findobj(gca, 'Type', 'image');  % get image object only
    set(hImg, 'AlphaData', alpha_data);
    mi_tmp = MI_landscaps_BG(:,:,n);
    caxis([min(mi_tmp(:), [], 'omitnan') max(mi_tmp(:), [], 'omitnan')]);

    %set(gca, 'YTickDir', 'reverse')
    c=colorbar;
    ylabel(c,'Bits');
    xline(0,'--','color',[1 1 1],'LineWidth',1.5)
    yline(0,'--','color',[1 1 1],'LineWidth',1.5)
    ylabel('start time (s)'); xlabel('stop time (s)')
    set(gca,axeOpt{:})
    axis square
    set(gcf,'position',[2025         321        1602         340],'color','w')
%    set(gcf,'position',[ 3165         321         462         340],'color','w')



end

    saveas(gcf,strcat(save_CB_neu,filesep,'MIland_4NEU_examples_BG',...
        num2str(n),'.png'),'png');
print(gcf,strcat(save_CB_neu,filesep,'MIland_4NEU_examples_BG.pdf'), '-dpdf', '-painters');