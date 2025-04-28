function [dist_trav_3D,speed_3D] = distance_travelled_and_speed(speed_partial,dt,px2mm)
% returns the distance travelled and speed from partial derivatives of
% position
% input: DD position derivative (across time)
%        delta time: tm(2)-tm(1) (int)
% output: distance travelled (in mm) and speed in 3D (1D vec, mm/s)
% teresa, 17/04/2023

% speed in 3D
% Default inputs
if nargin<3
    px2mm = [0.1 0.07 0.1]; %[px2mm_x px2mm_yz]
end
speed_partial_mm = speed_partial;
speed_partial_mm(:,1,:)=speed_partial(:,1,:).*px2mm(1);
speed_partial_mm(:,2,:)=speed_partial(:,2,:).*px2mm(2);
speed_partial_mm(:,2,:)=speed_partial(:,2,:).*px2mm(3);

speed_partial_dt = speed_partial_mm./dt;
    speed_3D = squeeze(sqrt(...
        speed_partial_dt(:,1,:).^2 + ...
        speed_partial_dt(:,2,:).^2 + ...
        speed_partial_dt(:,3,:).^2 ...
        ));
    
    % distance travelled
    dist_points = squeeze(sqrt(...
        speed_partial_mm(:,1,:).^2 + ...
        speed_partial_mm(:,2,:).^2 + ...
        speed_partial_mm(:,3,:).^2 ...
        ));
    dist_trav_3D = cumsum(dist_points,'omitnan');
end