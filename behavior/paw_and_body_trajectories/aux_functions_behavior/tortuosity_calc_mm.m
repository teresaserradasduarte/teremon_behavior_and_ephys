function [dist_trav,displacement,tortuosity] = tortuosity_calc_mm(xyz_vec,range,px2mm)
% returns the power sspectrum of a vector
% input: 3D signal (across time)
%        range = [start stop], optional, default = 1:end
% output: dist_trav,displacement,tortuosity
% teresa, 23/03/2023 (last update: 5/8/2023)

% Default inputs
if nargin<3
    px2mm = [0.1 0.07 0.1]; %[px2mm_x px2mm_yz]
    range = [1 size(xyz_vec,1)];
elseif nargin<2
    range = [1 size(xyz_vec,1)];
end

% Convert to mm - distance is not he same in all directions
xyz_vec_mm = xyz_vec;
xyz_vec_mm(:,1,:)=xyz_vec(:,1,:).*px2mm(1);
xyz_vec_mm(:,2,:)=xyz_vec(:,2,:).*px2mm(2);
xyz_vec_mm(:,3,:)=xyz_vec(:,3,:).*px2mm(3);

% distance travelled
start = range(1);
stop = range(2);
diff_r = diff(xyz_vec_mm(start:stop,:));
dist_points = squeeze(sqrt(...
    diff_r(:,1).^2 + ...
    diff_r(:,2).^2 + ...
    diff_r(:,3).^2 ...
    ));
dist_trav = sum(dist_points,'omitnan');
% displacement
displacement = sqrt(...
    (xyz_vec_mm(stop,1)-xyz_vec_mm(start,1)).^2+...
    (xyz_vec_mm(stop,2)-xyz_vec_mm(start,2)).^2+...
    (xyz_vec_mm(stop,3)-xyz_vec_mm(start,3)).^2);
% tortuosity
tortuosity = dist_trav/displacement;
end