function ax_figCount = figBehCountTrials(behavior, figProp)
% FIGURE: BEHAVIOR COUNTS OF TRIALS

%Figure properties
if nargin<2
    figProp.axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
        0.05,'ticklength',[1,1]*.01,'fontsize',10};
    figProp.lw = 2;
    figProp.transitionBarClr = [1 1 1 .5];
    figProp.lw_ppT = 3;
    figProp.lw_lcrT = 1;
    figProp.pa_pp_alpha = 0.3;
    figProp.pa_lcr_alpha = 0.6;
end

% plot counts
val = behavior.counts.count_valPull_per_unit + behavior.counts.count_valPush_per_unit;
inval = behavior.counts.count_invalPull_per_unit + behavior.counts.count_invalPush_per_unit;
max_yLim = round(max(val))+3;

%patches
pus=behavior.init.timeof.init_time(behavior.init.push_bounds)./60;
pul=behavior.init.timeof.init_time(behavior.init.pull_bounds)./60;
ll=behavior.init.timeof.init_time(behavior.reach.left_bounds)./60;
rr=behavior.init.timeof.init_time(behavior.reach.right_bounds)./60;
cc=behavior.init.timeof.init_time(behavior.reach.center_bounds)./60;

%transition bara
xlin_push = behavior.init.timeof.init_time(behavior.init.transition_2push)./60;
xlin_pull = behavior.init.timeof.init_time(behavior.init.transition_2pull)./60;
xlin_l = behavior.init.timeof.init_time(behavior.reach.moved_left)./60;
xlin_r = behavior.init.timeof.init_time(behavior.reach.moved_right)./60;
xlin_c = behavior.init.timeof.init_time(behavior.reach.moved_center)./60;

% FIGURE
%figure()

% patch push/pull
patch_y = [0 max_yLim max_yLim 0];
for l=1:size(pus,2)
    pa_pus = patch([pus(1,l) pus(1,l) pus(2,l) pus(2,l)],...
        patch_y,behavior.colors.push_clr,'FaceAlpha',figProp.pa_pp_alpha,'EdgeColor','none'); hold on
end
for l=1:size(pul,2)
    pa_pul = patch([pul(1,l) pul(1,l) pul(2,l) pul(2,l)],...
        patch_y,behavior.colors.pull_clr,'FaceAlpha',figProp.pa_pp_alpha,'EdgeColor','none'); hold on
end
ax = gca();
% trasition bars
plot([xlin_push,xlin_push]',ylim(ax), 'Color',figProp.transitionBarClr,'linewidth',figProp.lw_ppT);
plot([xlin_pull,xlin_pull]',ylim(ax), 'Color',figProp.transitionBarClr,'linewidth',figProp.lw_ppT);

% patch left/center/right
patch_yM= [max_yLim max_yLim+3 max_yLim+3 max_yLim];
for l=1:size(ll,2)
    pa_l = patch([ll(1,l) ll(1,l) ll(2,l) ll(2,l)],patch_yM,behavior.colors.left_color,...
        'FaceAlpha',figProp.pa_lcr_alpha ,'EdgeColor','none'); hold on
end
for l=1:size(rr,2)
    pa_r = patch([rr(1,l) rr(1,l) rr(2,l) rr(2,l)],patch_yM,behavior.colors.right_color,...
        'FaceAlpha',figProp.pa_lcr_alpha ,'EdgeColor','none'); hold on
end
for l=1:size(cc,2)
    pa_c = patch([cc(1,l) cc(1,l) cc(2,l) cc(2,l)],patch_yM,behavior.colors.center_color,'FaceAlpha',...
        figProp.pa_lcr_alpha ,'EdgeColor','none'); hold on
end

% transition bars
plot([xlin_l,xlin_l]',ylim(ax),'-.', 'Color',figProp.transitionBarClr,'linewidth',figProp.lw_lcrT);
plot([xlin_r,xlin_r]',ylim(ax),'-.', 'Color',figProp.transitionBarClr,'linewidth',figProp.lw_lcrT);
plot([xlin_c,xlin_c]',ylim(ax),'-.', 'Color',figProp.transitionBarClr,'linewidth',figProp.lw_lcrT);

% plot counts
p1=plot(behavior.counts.units_slide_min,val, 'k','LineWidth',figProp.lw); hold on
p2=plot(behavior.counts.units_slide_min,inval, '--','color',[.7 .7 .7],'LineWidth',figProp.lw);
axis ([0 behavior.session_duration/60 0 max_yLim+3])
ax_figCount = [p1,p2,pa_pus,pa_pul,pa_l,pa_r,pa_c];
set(gca,figProp.axeOpt{:})

% names and labels
% legend([p1,p2,pa_pus,pa_pul,pa_l,pa_r,pa_c],...
%     'valid counts','invalid counts','push','pull','left','right','center',...
%     'edgeColor','w','AutoUpdate','off','Location','eastoutside')
xlabel('time across session (min)'); ylabel('trial counts (3 min^{-1})');

end
