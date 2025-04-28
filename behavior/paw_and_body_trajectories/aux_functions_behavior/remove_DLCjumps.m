function s_reachX = remove_DLCjumps(xd_reach,th_jmp,flag_interpol)
% removes trajectory jumps typical from deeplabcut and interpolates (optional)
% (interpolation with pchip)
% input: 1D signal to clean
%        threshold for considering jump (optional, default = 100)
%        tinterpolate flag (optional, default = false)
% output: jump free, (optionally) interpolated 1D signal
% teresa, last update 4/03/2023


if nargin<3
    th_jmp = 100;
    flag_interpol = 0;
elseif nargin<2
    th_jmp = 100;
end

my_xd_reach = xd_reach; % input
no_nans_loc=find(~isnan(my_xd_reach)); % input without nans

% Correct for start/end/all points = nan
if numel(find(~isnan(my_xd_reach(:))))==0,my_xd_reach=ones(size(my_xd_reach)).*150; end % if x is all nan
if isnan(my_xd_reach(1)), my_xd_reach(1)=my_xd_reach(find(~isnan(my_xd_reach),1,'first')); end % if first point is nan
if isnan(my_xd_reach(end)), my_xd_reach(end)=my_xd_reach(find(~isnan(my_xd_reach),1,'last')); end % if last poins is nan
nan_xd_reach = my_xd_reach; % corrected (chnanged at end: jmps=nans)

% if the first dot is a jump, ignore it
range_for_med = 1:30;
if length(range_for_med)>length(my_xd_reach)
    med =nanmedian(my_xd_reach);
else

    med = nanmedian(my_xd_reach(range_for_med));
end

next_nearMed = find(abs(med-my_xd_reach)<th_jmp, 1);
if abs(med-my_xd_reach(1))>th_jmp
    my_xd_reach(1:next_nearMed-1)=nan;
end
diff_xd = diff(my_xd_reach); % derivative
[ind_jmp] = find(abs(diff_xd) > th_jmp); % find bigger jumps
siz_xd = numel(my_xd_reach);

% check if there are jumps that didn't come back
noBack_jmps = find(diff(sign(diff_xd(ind_jmp)))==0);

if ~isempty(noBack_jmps) % if that's the case...
    % if jmp is due to end of data, add jmp back at the end
    if ind_jmp(noBack_jmps(end))>=no_nans_loc(end-4)
        ind_jmp=[ind_jmp;no_nans_loc(end)]; %
        noBack_jmps = find(diff(sign(diff_xd(ind_jmp)))==0); % repeat after excluding end of data
    end
    
    % add come back jump, if it didn't return due to nans
    ind_jmp_initial = ind_jmp;
    flag_unknown_jmp = 0;
    added_stops = [];
    for jj=noBack_jmps'
        [find_nans_loc,~] = find(isnan(my_xd_reach(ind_jmp_initial(jj):ind_jmp_initial(jj+1))));
        if ~isempty(find_nans_loc)
            find_fist_jmp_nans = find(diff(find_nans_loc)>1,1,'first');
            if isempty(find_fist_jmp_nans)
                next_nonNan = ind_jmp_initial(jj)+find(diff(find_nans_loc)>=1,1,'first')+find_nans_loc(1)-1;
            else
                next_nonNan = ind_jmp_initial(jj)+find_fist_jmp_nans+find_nans_loc(1)-1;
            end
%             if abs(my_xd_reach(ind_jmp(ii))-my_xd_reach(next_nonNan))>25
%                 added_stops = cat(1,added_stops,next_nonNan);
%             end
        else
            warning('Warning! Jump')
            ind_jmp(jj)=[];
            flag_unknown_jmp = 1;
        end
    end
    % add returns
    ind_jmp_tmp = cat(1,ind_jmp,added_stops);
    ind_jmp = sort(ind_jmp_tmp);
end

% Start and stop
ind_start = ind_jmp(1:2:end);
ind_stop = ind_jmp(2:2:end);

if ~isempty(ind_stop)
    % Ignore jumps too large
    flag_jmp_too_long = 0;
    remove_start=[];
    remove_stop=[];
    
    for jj = 1:numel(ind_start)
        if numel(ind_stop)<jj
            remove_start=[remove_start;jj];
            warning('Warning! jump without return, uknown why... remove it')
        elseif numel(ind_start(jj):ind_stop(jj))>20
            remove_start=[remove_start;jj];
            remove_stop=[remove_stop;jj];
            flag_jmp_too_long = 1;
            warning('Warning! jump too large, uknown why... remove it')
        end
    end
    ind_start(remove_start)=[];
    ind_stop(remove_stop)=[];
    
    % replace jumpd by nan
    for jj = 1:numel(ind_start)
        nan_xd_reach(ind_start(jj):ind_stop(jj))=nan;
    end
end

% Interpolate jumps with pchip if interpolation asked for
ff=1:siz_xd;
if (mod(numel(ind_start)+numel(ind_stop),2)==1) % if after everything, is still odd, just don't do it
    s_reachX = my_xd_reach;
    warning('Warning! odd number of jumps, unknown reason. Ignored jump removal')
else
    if flag_interpol
        s_reachX = pchip(ff,nan_xd_reach,ff);
    else
        s_reachX = nan_xd_reach;
    end
end
% %
% % Plot
% figure
% plot(my_xd_reach,'.'); shg
% hold on
% plot(repmat(med,[1 10]));
% plot(diff_xd)
% %plot(abs(diff_xd))
% plot(ind_jmp,diff_xd(ind_jmp),'*')
% %plot(nan_xd_reach,'.')
% plot(s_reachX,'.')
% hold off

end




