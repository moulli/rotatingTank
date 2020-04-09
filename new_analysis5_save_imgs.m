clear; close all; clc
folderpath = '/home/ljp/Science/Hippolyte/rotatingTank_imgs';


%% Save all required information from the videos

% Define parameters
fish = 1:18;
light = {' Off', ' On'};
rot = 0:2:12;

% Loop over all fish, light and rotation
tic
for f = fish
    for l = light
        for r = rot
            
            %% Get path to video
            
            vidpath = strcat('/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish', ...
                              num2str(f), '/Light', l{1}, '/Films/Fish' , l{1}, '_', ...
                              num2str(f), '_', num2str(r), '.avi');
            
            
            %% Define video reader
            
            vid = VideoReader(vidpath);
            numframes = floor(vid.Duration * vid.FrameRate);
            
            
            %% Save each frame

            for frame = 1:numframes
                % get frame and make it square
                im = readFrame(vid);
                im = im(5:(size(im, 1)-4), :);
                % define name of figure
                fstr = num2str(f); if length(fstr) == 1; fstr = strcat('0', fstr); end
                if isequal(l{1}, ' On'); lstr = 'ON'; elseif isequal(l{1}, ' Off'); lstr = 'OF'; end
                rstr = num2str(r); if length(rstr) == 1; rstr = strcat('0', rstr); end
                framestr = num2str(frame);
                framestr = strcat(repmat('0', 1, 4-length(framestr)), framestr);
                imgname = strcat(fstr, lstr, rstr, framestr, '.png');
                imgpath = fullfile(folderpath, imgname);
                % save figure
                imwrite(im, imgpath, 'png')
            end
            
            
            %% Print indication
            
            message = strcat('Fish:', num2str(f), ', lights:', l{1}(2:end), ', rotation:', num2str(r), ', done in:', num2str(toc), 's \n');
            fprintf(message);
            
            
        end
    end
end