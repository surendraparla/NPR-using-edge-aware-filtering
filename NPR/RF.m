%  RF  Domain transform recursive filter edge-preserving filter.
%  F = NC(img, sigma_s, sigma_r, num_iter)

%  Parameters:
%    image           Input image to be filtered.
%    sigma_s         Filter spatial standard deviation.
%    sigma_r         Filter range standard deviation.
%    num_iter        Number of iterations to perform (default: 3).
function F = RF(image, sigma_s, sigma_r, num_iter)

    I = double(image);

    if ~exist('num_iter', 'var')
        num_iter = 3;
    end
    [h, w, num_channels] = size(I);

    % Compute the domain transform (Equation 11 of our paper)
    % Estimate horizontal and vertical partial derivatives using finite
    % differences.
    dIcdx = diff(I, 1, 2);
    dIcdy = diff(I, 1, 1);
    
    dIdx = zeros(h,w);
    dIdy = zeros(h,w);
    
    % Compute the l1-norm distance of neighbor pixels.
    for c = 1:num_channels
        dIdx(:,2:end) = dIdx(:,2:end) + abs( dIcdx(:,:,c) );
        dIdy(2:end,:) = dIdy(2:end,:) + abs( dIcdy(:,:,c) );
    end
    
    % Compute the derivatives of the horizontal and vertical domain transforms.
    dHdx = (1 + sigma_s/sigma_r * dIdx);
    dVdy = (1 + sigma_s/sigma_r * dIdy);
    
    % We do not integrate the domain transforms since our recursive filter
    % uses the derivatives directly.
    %ct_H = cumsum(dHdx, 2);
    %ct_V = cumsum(dVdy, 1);
    
    % The vertical pass is performed using a transposed image.
    dVdy = dVdy';
    
    % Perform the filtering.
    N = num_iter;
    F = I;
    
    sigma_H = sigma_s;
    
    for i = 0:num_iter - 1
    
        % Compute the sigma value for this iteration (Equation 14 of our paper).
        sigma_H_i = sigma_H * sqrt(3) * 2^(N - (i + 1)) / sqrt(4^N - 1);
    
        F = TransformedDomainRecursiveFilter_Horizontal(F, dHdx, sigma_H_i);
        F = image_transpose(F);
    
        F = TransformedDomainRecursiveFilter_Horizontal(F, dVdy, sigma_H_i);
        F = image_transpose(F);
        
    end
    
    F = cast(F, class(image));

end

% Recursive filter
function F = TransformedDomainRecursiveFilter_Horizontal(I, D, sigma)

    % Feedback coefficient (Appendix of our paper).
    a = exp(-sqrt(2) / sigma);
    
    F = I;
    V = a.^D;
    
    [ h , w, num_channels] = size(I);
    
    % Left -> Right filter.
    for i = 2:w
        for c = 1:num_channels
            F(:,i,c) = F(:,i,c) + V(:,i) .* ( F(:,i - 1,c) - F(:,i,c) );
        end
    end
    
    % Right -> Left filter.
    for i = w-1:-1:1
        for c = 1:num_channels
            F(:,i,c) = F(:,i,c) + V(:,i+1) .* ( F(:,i + 1,c) - F(:,i,c) );
        end
    end

end

function T = image_transpose(I)

    [h, w, num_channels] = size(I);
    
    T = zeros([w h num_channels], class(I));
    
    for c = 1:num_channels
        T(:,:,c) = I(:,:,c)';
    end
    
end
