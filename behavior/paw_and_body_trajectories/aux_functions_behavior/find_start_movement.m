function [rising_from_stop,onset_stop,binary_moving] = find_start_movement(reach_speed,close_zero,min_frames_down_toStop,reaching_frames_minimum)
% return the moment of movement start
% input: 1D speed signal (across time)
%        close_zero: minumum speed to consider stopped (int)
%        min_frames_down_toStop: number of frames below close_zero (int)
%        min_frames_down_toStop: number of frames below close_zero (int)
%        min_reach_time: minimum time given to reach, can't be stopped
% teresa, 18/03/2023


if nargin<4
    close_zero = 1;
    min_frames_down_toStop = 10;
    reaching_frames_minimum = 130;
elseif nargin<3
    close_zero = 1;
    min_frames_down_toStop = 10;
elseif nargin<2
    close_zero = 1;
end

% ONSET FROM STOPPED POSITION
at2 =reach_speed;
azt2=zeros(size(at2)); azt2(at2>close_zero)=1; azt2(reaching_frames_minimum:end)=1;
afzt2 = flip(azt2,1); afzt2(end)=1; % Flip zt2 to find the first stable rise timepoint
if isempty(find(diff(afzt2)~=0,1))
    onset_stop = 0;
    rising_from_stop = nan;
else
    s_fjumps_loc=find(diff(afzt2)~=0); s_fjumps_down_loc=find(diff(afzt2)==-1); % find the position of jumps down
    s_diff_jumps=diff(s_fjumps_loc); s_inv_rise_jumps_diff = s_diff_jumps(1:2:end); % find number of frames after jump
    stop_inv_stable_rising_start=s_fjumps_down_loc(s_inv_rise_jumps_diff>min_frames_down_toStop); % Impose a min number of frames to consider stalbe rise
    if isempty(stop_inv_stable_rising_start)
        onset_stop = 0;
        rising_from_stop = nan;
        
    else
        rising_from_stop=numel(afzt2)-stop_inv_stable_rising_start(1); % Invert back
        onset_stop = 1;
    end
end
binary_moving = azt2;

end