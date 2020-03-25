function [background, masks] = getMasks(vidpath, varargin)

%% getMasks returns background and set of masks for a specified video
%
%  First of all, getMasks build background, by averaging images taken each
%  10 frames from the video (using VideoReader). Then the circle border of
%  the tank is detected to have its center. Using different radii, we
%  finally compute a set of masks that will be useful to detect the fish.
%
%  Inputs:
%  - vidpath [string]: link to video.
%
%  Outputs:
%  - background [2D matrix]: averaged image from video.
%  - masks [3D matrix]: 2D matrices built with 0s and 1s, and stacked along
%    third dimension.


    %% Define video reader and get numframes
    
    % Define videoreader
    vid = VideoReader(vidpath);
    width = vid.Width;
    height = vid.Height;
    
    % Get numframes
    defaultNumFrames = floor(vid.Duration * vid.FrameRate);
    
    % Input parser
    p = inputParser;
    addRequired(p, 'vidpath');
    addOptional(p, 'numframes', defaultNumFrames);
    parse(p, vidpath, varargin{:});
    
    
    %% Compute background

    % Initialize background image
    background = 0; 

    % Add frame every dframe
    dframe = 10;
    for frame = 1:dframe:p.Results.numframes
        im = read(vid, frame);
        im = flipud(double(im)); % image is upside down on matlab
        background = background + im;
    end

    % Compute final background
    background = (background./(p.Results.numframes/dframe));
    
    
    %% Identify circle with highest radius

    % Get circle with highest radius
    [centers, radii] = imfindcircles(background, [round((height)/2.3), round((height)/2)], 'ObjectPolarity', 'dark', 'Sensitivity', 0.99);
    radius = max(radii);
    center = centers(radius == radii, :);
    
    
    %% Compute list of masks
    
    % Define meshgrid for masks
    [x, y] = meshgrid(-(center(2)-1):(width-center(2)), -(center(1)-1):(height-center(1))); % meshgrid for mask
    % ATTENTION: there is a confusion between center(1) and center(2) in the
    % above equation. Lower, when defining another mask, the y-center was in
    % the left part of meshgrid, and the x-center was in the right part of the
    % meshgrid, so I switched it as well here to be sure...

    % Computing averaged pixel value for different radii
    mult_coeffs = 0.80:0.001:1.2;
    pixradius = zeros(size(mult_coeffs));
    for cmul = 1:length(mult_coeffs)
        cmasktemps = ((x.^2 + y.^2) <= (mult_coeffs(cmul)*radius)^2);
        pixradius(cmul) = sum(sum(background .* cmasktemps)) / sum(sum(cmasktemps));
    end

    % Get radius before minimum differential of pixradius
    [~, mindiff] = min(diff(pixradius));
    radius_coeff = mult_coeffs(mindiff);

    % Create masks
    radius_coeff_inc = 0.01;
    radii = ((radius_coeff-10*radius_coeff_inc):radius_coeff_inc:(radius_coeff+5*radius_coeff_inc)) .* radius;
    masks = zeros(height, width, length(radii));
    for r = 1:length(radii)
        maskr = ((x.^2 + y.^2) <= radii(r)^2);
        masks(:, :, r) = maskr;
    end
    

end