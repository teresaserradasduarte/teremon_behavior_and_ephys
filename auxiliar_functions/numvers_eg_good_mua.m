%% Numbers
clear; close all; 

%%
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\';
group = '20230801_ChocolateGroup';
animals = {...
    'CoteDor',...
    'Lindt',...
    'Toblerone',...
    'Milka',...
    'FerreroRocher'};
n_mice = length(animals);
n_sess = 7;
%%
datapath = fullfile(rootdir,"ephys_and_behavior","mat_files",group);
load(fullfile(datapath,'eg_neurons.mat'),'n_eg_neu','n_good_neu','n_mua_neu');

%%
sorter = zeros(n_mice,n_sess,2);
sorter(isnan(n_eg_neu))=nan;
sorter(1,2,1)=1;
sorter(1,6,1)=1;
sorter(1,1,2)=1;
sorter(1,6,2)=1;

sorter(2,1,1)=1;
sorter(2,2,1)=1;
sorter(2,4,1)=1;
sorter(2,1,2)=1;
sorter(2,2,2)=1;
sorter(2,4,2)=1;

sorter(3,3,1)=1;
sorter(3,6,1)=1;
sorter(3,3,2)=1;
sorter(3,6,2)=1;

sorter(4,4,1)=1;
sorter(4,4,2)=1;

sorter(5,1,1)=1;
sorter(5,4,1)=1;
sorter(5,6,1)=1;
sorter(5,4,2)=1;
sorter(5,6,2)=1;
sorter(5,1,2)=1;

%%
save_out = fullfile(rootdir,"ephys_and_behavior","out_files",group,'group','numbers');
if ~exist(save_out,"dir"), mkdir(save_out); end

%%
n_eg_bg = n_eg_neu(:,:,1);
n_eg_cb = n_eg_neu(:,:,2);

sorter_bg = sorter(:,:,1);
sorter_cb = sorter(:,:,2);

n_eg_bg_ks = n_eg_bg(sorter_bg==1);
n_eg_bg_ib = n_eg_bg(sorter_bg==0);

n_eg_cb_ks = n_eg_cb(sorter_cb==1);
n_eg_cb_ib = n_eg_cb(sorter_cb==0);

max_ks = max(length(n_eg_bg_ks),length(n_eg_cb_ks));
max_ib = max(length(n_eg_bg_ib),length(n_eg_cb_ib));

vect_eg_ks = nan(max_ks,2);
vect_eg_ks(1:length(n_eg_bg_ks),1) = n_eg_bg_ks;
vect_eg_ks(1:length(n_eg_cb_ks),2) = n_eg_cb_ks;

vect_eg_ib = nan(max_ib,2);
vect_eg_ib(1:length(n_eg_bg_ib),1) = n_eg_bg_ib;
vect_eg_ib(1:length(n_eg_cb_ib),2) = n_eg_cb_ib;


vec_eg = cat(1,vect_eg_ks,vect_eg_ib);

clrs = lines(2);
clrs_a = cat(1,repmat(clrs(1,:),[max_ks 1]),repmat(clrs(2,:),[max_ib 1]));
legth_vec = size(vec_eg,1);

%%
n_g_bg = n_good_neu(:,:,1);
n_g_cb = n_good_neu(:,:,2);

n_g_bg_ks = n_g_bg(sorter_bg==1);
n_g_bg_ib = n_g_bg(sorter_bg==0);

n_g_cb_ks = n_g_cb(sorter_cb==1);
n_g_cb_ib = n_g_cb(sorter_cb==0);

vect_g_ks = nan(max_ks,2);
vect_g_ks(1:length(n_g_bg_ks),1) = n_g_bg_ks;
vect_g_ks(1:length(n_g_cb_ks),2) = n_g_cb_ks;

vect_g_ib = nan(max_ib,2);
vect_g_ib(1:length(n_g_bg_ib),1) = n_g_bg_ib;
vect_g_ib(1:length(n_g_cb_ib),2) = n_g_cb_ib;

vec_g = cat(1,vect_g_ks,vect_g_ib);

%%
n_m_bg = n_mua_neu(:,:,1);
n_m_cb = n_mua_neu(:,:,2);

n_m_bg_ks = n_m_bg(sorter_bg==1);
n_m_bg_ib = n_m_bg(sorter_bg==0);

n_m_cb_ks = n_m_cb(sorter_cb==1);
n_m_cb_ib = n_m_cb(sorter_cb==0);

vect_m_ks = nan(max_ks,2);
vect_m_ks(1:length(n_m_bg_ks),1) = n_m_bg_ks;
vect_m_ks(1:length(n_m_cb_ks),2) = n_m_cb_ks;

vect_m_ib = nan(max_ib,2);
vect_m_ib(1:length(n_m_bg_ib),1) = n_m_bg_ib;
vect_m_ib(1:length(n_m_cb_ib),2) = n_m_cb_ib;

vec_m = cat(1,vect_m_ks,vect_m_ib);



vex_eg_g_m = cat(2,vec_eg,nan(legth_vec,1),vec_g,nan(legth_vec,1),vec_m);

%%
%nxh_plot=8;

figure(10)
boxplot(vex_eg_g_m,'Symbol', 'k.','Color','k','Widths',0.8);
hold on
for pos=1:size(vex_eg_g_m,2)
    %if ismember(pos,1:3:nxh_plot*3), clrs = clrs_BG; else, clrs = clrs_CB; end
    if pos ==1
        f=scatter((1+(rand(01,legth_vec)-0.5)/5),vex_eg_g_m(:,pos),50,clrs_a,'filled','LineWidth',1.5);
    else
        f=scatter(ones(1,legth_vec)+(pos-1).*(1+(rand(01,legth_vec)-0.5)/20),vex_eg_g_m(:,pos),50,clrs_a,'filled','LineWidth',1.5);
    end
end
hold off
xticks([1,2,4,5,7,8])
xticklabels({'BG';'CB';'BG';'CB';'BG';'CB'})
ylabel('# clusters')
set(gca,'linewidth',1.5,'box','off','ticklength',[1,1]*.01)
set(gcf,'position',[680   495   860   483],'color','w')

saveas(gcf,fullfile(save_out,'eg_good_mua.png'),'png')

%%

figure(10)
boxplot(vec_eg,'Symbol', 'k.','Color','k','Widths',0.8);
hold on
for pos=1:size(vec_eg,2)
    %if ismember(pos,1:3:nxh_plot*3), clrs = clrs_BG; else, clrs = clrs_CB; end
    if pos ==1
        f=scatter((1+(rand(01,legth_vec)-0.5)/5),vec_eg(:,pos),50,clrs_a,'filled','LineWidth',1.5);
    else
        f=scatter(ones(1,legth_vec)+(pos-1).*(1+(rand(01,legth_vec)-0.5)/20),vec_eg(:,pos),50,clrs_a,'filled','LineWidth',1.5);
    end
end
hold off
xticks([1,2,])
xticklabels({'BG';'CB'})
ylabel('# clusters')
set(gca,'linewidth',1.5,'box','off','ticklength',[1,1]*.01)
%set(gcf,'position',[680   495   860   483],'color','w')

saveas(gcf,fullfile(save_out,'eg_only.png'),'png')


















