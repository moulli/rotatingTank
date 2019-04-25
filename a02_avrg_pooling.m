function out = a02_avrg_pooling(in, stride)
    [xin, yin] = size(in);
    xout = ceil(xin / stride);
    yout = ceil(yin / stride);
    out = zeros(xout, yout);
    for i = 1:xout
        for j = 1:yout
            if i == xout && j == yout
                outemp = in((i-1)*stride+1:end, (j-1)*stride+1:end);
            elseif i == xout
                outemp = in((i-1)*stride+1:end, (j-1)*stride+1:j*stride);
            elseif j == yout
                outemp = in((i-1)*stride+1:i*stride, (j-1)*stride+1:end);
            else
                outemp = in((i-1)*stride+1:i*stride, (j-1)*stride+1:j*stride);
            end
            out(i, j) = mean(outemp(:));
        end
        out(i, j) = mean(outemp(:));
    end
end