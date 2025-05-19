function cmap = bwr_map(n)
    if nargin < 1
        n = 256; % default number of colors
    end

    % Define control points: blue → white → red
    cdict = [
        0,   0,   1;   % Blue
        1,   1,   1;   % White
        1,   0,   0    % Red
    ];

    % Interpolate over red, green, blue
    x = linspace(0, 1, size(cdict,1));
    xi = linspace(0, 1, n);

    cmap = [ ...
        interp1(x, cdict(:,1), xi)', ...
        interp1(x, cdict(:,2), xi)', ...
        interp1(x, cdict(:,3), xi)' ...
    ];
end