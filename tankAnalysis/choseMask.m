function outmask = choseMask(im, background, masks, formerdp, max_diff_disp)

%% choseMask returns highest pixel and selected mask for image im
%
%  If it is not the first frame, choseMask checks for every mask the point
%  with the highest pixel value. It then compares the distance between this
%  new point for each mask to the former distance. The new highest pixel
%  with a distance lower than the threshold gets the new used mask. Instead
%  of taking the exact mask returning this pixel, we take a margin of 3
%  masks, in case fish tail it outside mask.
%
%  Inputs:
%  - im [2D matrix]: image to analyze.
%  - background [2D matrix]: background of video.
%  - masks [3D matrix]: stack of masks (along 3rd dimension).
%  - formerdp [2D vector]: coordinates of former highest pixel.
%  - max_diff_disp [number]: maximum distance allowed between highest
%    pixels at two time periods.
%
%  Outputs:
%  - newdp [2D vector]: coordinates of new highest pixel.
%  - outmask [integer]: index of mask retained.


    %% Convert image to double
    
    im = double(im);
        
        
    %% Else go through masks and check distance between highest pixels
        
    % Try all masks
    nmask = size(masks, 3);
    dpmaskx = zeros(nmask, 1);
    dpmasky = zeros(nmask, 1);
    pixmask = zeros(nmask, 1);
    for mask = 1:nmask      
        % Get highest pixel
        new_im = (background - im) .* masks(:, :, mask);
        [dpframex, dpframey] = find(new_im == max(new_im(:)));
        dpmaskx(mask) = dpframex(1);
        dpmasky(mask) = dpframey(1);
        pixmask(mask) = max(new_im(:));
    end

    % Take mask with highest pixel value, and displacement not above threshold
    [~, pixind] = sort(pixmask, 'descend');
    for pixi = 1:length(pixind)
        pix = pixind(pixi);
        diff_disp = sqrt((dpmaskx(pix)-formerdp(1)).^2 + (dpmasky(pix)-formerdp(2)).^2);
        if abs(diff_disp) < max_diff_disp
            outmask = min([pix+1, length(pixind)]); % take mask just above just in case
            break
        end
        if pixi == length(pixind)
            pix = pixind(1);
            outmask = pix;
        end
    end
    
    
end