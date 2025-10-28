function landmarks_plt=behavior_landmarks(behavior, figProp,patchY)
% FIGURE: BEHAVIOR COUNTS OF TRIALS

%Figure properties
if nargin<3
patchY.pp = [0 50 50 0];
patchY.rlc = [50 53 53 50];
patchY.extend = 1;
patchY.showPatch = 1;

elseif nargin<2
    figProp.axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
        0.05,'ticklength',[1,1]*.01,'fontsize',10};
    figProp.lw = 2;
    figProp.transitionBarClr = [1 1 1 .5];
    figProp.lw_ppT = 3;
    figProp.lw_lcrT = 1;
    figProp.pa_pp_alpha = 0.3;
    figProp.pa_lcr_alpha = 0.6;
end

if patchY.extend
    %patches
    pus=behavior.init.timeof.init_time(behavior.init.push_bounds)./60;
    pul=behavior.init.timeof.init_time(behavior.init.pull_bounds)./60;
    if (~isempty(pul) && ~isempty(pus))
        if behavior.init.pull_bounds(end)>behavior.init.push_bounds(end)
            pul(end)=behavior.session_duration/60;
        elseif behavior.init.pull_bounds(end)<behavior.init.push_bounds(end)
            pus(end)=behavior.session_duration/60;
        end
    end
end

ll=behavior.init.timeof.init_time(behavior.reach.left_bounds)./60;
rr=behavior.init.timeof.init_time(behavior.reach.right_bounds)./60;
cc=behavior.init.timeof.init_time(behavior.reach.center_bounds)./60;
if patchY.extend
    if (ll(end)>rr(end) && ll(end)>cc(end))
        ll(end)=behavior.session_duration/60;
    elseif (rr(end)>ll(end) && rr(end)>cc(end))
        rr(end)=behavior.session_duration/60;
    elseif (cc(end)>ll(end) && cc(end)>rr(end))
        cc(end)=behavior.session_duration/60;
    end
end

%transition bara
xlin_push = behavior.init.timeof.init_time(behavior.init.transition_2push)./60;
xlin_pull = behavior.init.timeof.init_time(behavior.init.transition_2pull)./60;
xlin_l = behavior.init.timeof.init_time(behavior.reach.moved_left)./60;
xlin_r = behavior.init.timeof.init_time(behavior.reach.moved_right)./60;
xlin_c = behavior.init.timeof.init_time(behavior.reach.moved_center)./60;

% FIGURE
%figure()
ax = gca();
% transition bars
plot([xlin_l,xlin_l]',ylim(ax),'-.', 'Color',figProp.transitionBarClr,'linewidth',figProp.lw_lcrT);
plot([xlin_r,xlin_r]',ylim(ax),'-.', 'Color',figProp.transitionBarClr,'linewidth',figProp.lw_lcrT);
plot([xlin_c,xlin_c]',ylim(ax),'-.', 'Color',figProp.transitionBarClr,'linewidth',figProp.lw_lcrT);

% trasition bars
if ~isempty(xlin_push), plot([xlin_push,xlin_push]',ylim(ax), 'Color',figProp.transitionBarClr,'linewidth',figProp.lw_ppT); end
if ~isempty(xlin_pull), plot([xlin_pull,xlin_pull]',ylim(ax), 'Color',figProp.transitionBarClr,'linewidth',figProp.lw_ppT); end

if patchY.showPatch
    % patch push/pull
    if ~isempty(pus)
        for l=1:size(pus,2)
            pa_pus = patch([pus(1,l) pus(1,l) pus(2,l) pus(2,l)],...
                patchY.pp,behavior.colors.push_clr,'FaceAlpha',figProp.pa_pp_alpha,'EdgeColor','none'); hold on
        end
    else
        pa_pus=[];
    end
    if ~isempty(pul)
        for l=1:size(pul,2)
            pa_pul = patch([pul(1,l) pul(1,l) pul(2,l) pul(2,l)],...
                patchY.pp,behavior.colors.pull_clr,'FaceAlpha',figProp.pa_pp_alpha,'EdgeColor','none'); hold on
        end
    else
        pa_pul=[];
    end
    % patch left/center/right
    for l=1:size(ll,2)
        pa_l = patch([ll(1,l) ll(1,l) ll(2,l) ll(2,l)],patchY.rlc,behavior.colors.left_color,...
            'FaceAlpha',figProp.pa_lcr_alpha ,'EdgeColor','none'); hold on
    end
    for l=1:size(rr,2)
        pa_r = patch([rr(1,l) rr(1,l) rr(2,l) rr(2,l)],patchY.rlc,behavior.colors.right_color,...
            'FaceAlpha',figProp.pa_lcr_alpha ,'EdgeColor','none'); hold on
    end
    for l=1:size(cc,2)
        pa_c = patch([cc(1,l) cc(1,l) cc(2,l) cc(2,l)],patchY.rlc,behavior.colors.center_color,'FaceAlpha',...
            figProp.pa_lcr_alpha ,'EdgeColor','none'); hold on
    end
else
    pa_pus= nan;
    pa_pul= nan;
    pa_l= nan;
    pa_r= nan;
    pa_c = nan;
end
landmarks_plt = [pa_pus,pa_pul,pa_l,pa_r,pa_c];

set(gca,figProp.axeOpt{:})


end
