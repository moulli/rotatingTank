clear; close all; clc
addpath(genpath('/home/ljp/Science/balancoire_data/balancoire_matlab_files'))
addpath(genpath('/home/ljp/Science/Hippolyte/rotatingTank'))


%% Save all required information from the videos

% Define parameters
fish = 1:18;
light = {' Off', ' On'};
rot = 0:2:12;
folderpath = '/home/ljp/Science/rotatingTank_vidMatrices';

% Loop over all fish, light and rotation
tic
for f = fish
    for l = light
        for r = rot
            
            %% Get path to video
            
            vidpath = strcat('/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish', ...
                              num2str(f), '/Light', l, '/Films/Fish' , l, '_', ...
                              num2str(f), '_', num2str(r), '.avi');
            vidpath = vidpath{1};
    
                          
            %% Prepare output of file

            washingOut = struct;
            washingOut.parameters.fish = f;
            washingOut.parameters.light = l{1}(2:end); % get rid of space
            washingOut.parameters.rotation = r;


            %% Create video reader

            vid = VideoReader(vidpath);
            numframes = floor(vid.Duration * vid.FrameRate);
            width = vid.Width;
            height = vid.Height;

            % Add to output structure
            washingOut.parameters.path = vidpath;
            washingOut.parameters.numframes = numframes;
            washingOut.duration = vid.Duration;
            washingOut.framerate = vid.FrameRate;
            washingOut.parameters.width = width;
            washingOut.parameters.height = height;


            %% Get background and masks

            [background, masks] = getMasks(vidpath, 'numframes', numframes);

            % Add to output structure
            washingOut.preprocessing.background = background;
            washingOut.preprocessing.masks = masks;


            %% Loop over video and save neighbourhood of fish

            % Redefine video reader
            vid = VideoReader(vidpath);

            % Define maximum displacement differential to switch mask in needed
            max_diff_disp = 20;
            % Define neighbourhood to consider around fish
            rneighb = 45;

            % Define outputs
            dp = zeros(numframes, 2);
            maskused = zeros(numframes, 1);

            % Fish neighbourhood
            fish_neighb = zeros(2*rneighb+1, 2*rneighb+1, numframes);

            % Loop over all frames
            for frame = 1:numframes

                % Read frame
                im = readFrame(vid);
                im = flipud(double(im));

                % Chose right mask
                if frame == 1
                    outmask = 8; % arbitrary number for first mask
                else
                    outmask = choseMask(im, background, masks, dp(frame-1, :), max_diff_disp);
                end
                maskused(frame) = outmask;

                % Get highest pixel
                new_im = (background - im) .* masks(:, :, outmask);
                [dpframex, dpframey] = find(new_im == max(new_im(:)));
                newdp = [dpframex(1), dpframey(1)];
                dp(frame, :) = newdp;
                % Get wider mask, in case tail is outside ROI
                outmask = min([outmask+2, size(masks, 3)]);
                new_im = (background - im) .* masks(:, :, outmask);

                % Make a mask to isolate neighbourhood
                [xnei, ynei] = meshgrid(-(newdp(2)-1):(width-newdp(2)), -(newdp(1)-1):(height-newdp(1))); 
                masknei = ((xnei.^2 + ynei.^2) <= rneighb^2);
                new_im = masknei .* new_im;

                % Pad new_im
                new_im = padarray(new_im, [rneighb, rneighb], 0, 'both');
                % Save neighbourhood
                fish_neighb(:, :, frame) = new_im(newdp(1):newdp(1)+2*rneighb, newdp(2):newdp(2)+2*rneighb);

            end

            % POST TREATMENT: x axis is y axis on matlab
            dp = fliplr(dp);

            % Add to output structure
            washingOut.images.parameters.max_diff_disp = max_diff_disp;
            washingOut.images.parameters.rneighb = rneighb;
            washingOut.images.dp = dp;
            washingOut.images.maskused = maskused;
            washingOut.images.fish_neighb = fish_neighb;


            %% Loop over video and save neighbourhood of fish

            % Segment tracking parameters
            num_segments = 5;
            tail_length = 30;
            body_length = 0;
            inertia = 0.5;
            num_pix1 = 10;
            num_pix2 = 150;
            initial_box = 0.6;

            % Define outputs
            fpts = cell(numframes, 1);

            % Loop over all frames
            for frame = 1:numframes

                % Read frame
                new_im = fish_neighb(:, :, frame);

                % Identifying important points
                new_im = (255 - new_im);
                segment_pts = segmentTracking_WM(new_im, 'num_segments', num_segments, ...
                                    'tail_length', tail_length, 'body_length', body_length, ...
                                    'inertia', inertia, 'num_pix1', num_pix1, ...
                                    'num_pix2', num_pix2, 'initial_box', initial_box);
                fpts{frame} = segment_pts;

            end

            % POST TREATMENT: x axis is y axis on matlab
            for frame = 1:numframes
                fpts{frame} = fliplr(fpts{frame});
            end

            % POST TREATMENT: remove anomaly points
            % If orientation for x-coordinate and y-coordinate is inverted for only one
            % time point, we will say this time point is the mean of the former and the
            % later position
            modified = zeros(numframes, 1);
            for frame = 2:(numframes-1)
                % Orientation from before
                ori_m1 = [fpts{frame-1}(2, 1)-fpts{frame-1}(1, 1), fpts{frame-1}(2, 2)-fpts{frame-1}(1, 2)];
                ori_m1 = (ori_m1 > 0);
                % Orientation from now
                ori_n = [fpts{frame}(2, 1)-fpts{frame}(1, 1), fpts{frame}(2, 2)-fpts{frame}(1, 2)];
                ori_n = (ori_n < 0);
                % Orientation from before
                ori_p1 = [fpts{frame+1}(2, 1)-fpts{frame+1}(1, 1), fpts{frame+1}(2, 2)-fpts{frame+1}(1, 2)];
                ori_p1 = (ori_p1 > 0);
                % Correct if necessary
                if isequal(ori_m1, ori_n) && isequal(ori_n, ori_p1)
                    modified(frame) = 1;
                    fpts{frame}(2:end, :) = (fpts{frame-1}(2:end, :) + fpts{frame+1}(2:end, :)) / 2;
                end    
            end

            % Get info from fpts
            ptss = zeros(num_segments+1, 2, numframes);
            for frame = 1:numframes
                for pt = 1:num_segments+1
                    for coord = 1:2
                        ptss(pt, coord, frame) = fpts{frame}(pt, coord) + dp(frame, coord);
                        % NB: equivalent to ptss(:, :, frame) = fpts{frame} + dp(frame, :);
                    end
                end
            end

            % Add to output structure
            washingOut.tracking.parameters.num_segments = num_segments;
            washingOut.tracking.parameters.tail_length = tail_length;
            washingOut.tracking.parameters.body_length = body_length;
            washingOut.tracking.parameters.inertia = inertia;
            washingOut.tracking.parameters.num_pix1 = num_pix1;
            washingOut.tracking.parameters.num_pix2 = num_pix2;
            washingOut.tracking.parameters.initial_box = initial_box;
            washingOut.tracking.fpts = fpts;
            washingOut.tracking.modified = modified;
            washingOut.tracking.ptss = ptss;


            %% Save file in corresponding folder

            savepath = fullfile(folderpath, strcat('washingOut_', num2str(f), '_', l{1}(2:end), '_', num2str(r)));
            save(savepath, 'washingOut');
            
            
            %% Print indication
            
            message = strcat('Fish:', num2str(f), ', lights:', l{1}(2:end), ', rotation:', num2str(r), ', done in:', num2str(toc), 's \n');
            fprintf(message);
            
        end
    end
end
            
            
            
            