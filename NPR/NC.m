%  NC  Domain transform normalized convolution edge-preserving filter.
%  F = NC(img, sigma_s, sigma_r, num_iterations)

%  Parameters:
%    image             Input image to be filtered.
%    sigma_s         Filter spatial standard deviation.
%    sigma_r         Filter range standard deviation.
%    num_iter  Number of iterations to perform (default: 3).

function F = NC(image, sigma_s, sigma_r, num_iter)

    I = double(image);

    if ~exist('num_iterations', 'var')
        num_iter = 3;
    end   
    [h, w, num_channels] = size(I);

    % Compute the domain transform 
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
    
    % Integrate the domain transforms.
    ct_H = cumsum(dHdx, 2);
    ct_V = cumsum(dVdy, 1);
    
    % The vertical pass is performed using a transposed image.
    ct_V = ct_V';
    
    % Perform the filtering.
    N = num_iter;
    F = I;

    sigma_H = sigma_s;
    
    for i = 0:num_iter - 1
        
        % Compute the sigma value for this iteration (Equation 14 of our paper).
        sigma_H_i = sigma_H * sqrt(3) * 2^(N - (i + 1)) / sqrt(4^N - 1);
        
        % Compute the radius of the box filter with the desired variance.
        box_radius = sqrt(3) * sigma_H_i;
        
        F = TransformedDomainBoxFilter_Horizontal(F, ct_H, box_radius);
        F = image_transpose(F);

        F = TransformedDomainBoxFilter_Horizontal(F, ct_V, box_radius);
        F = image_transpose(F);
        
    end
    
    F = cast(F, class(image));

end
%Box filter normalized convolution in the transformed domain.
function F = TransformedDomainBoxFilter_Horizontal(I, xform_domain_position, box_radius)

    [h, w, num_channels] = size(I);

    % Compute the lower and upper limits of the box kernel in the transformed domain.
    l_pos = xform_domain_position - box_radius;
    u_pos = xform_domain_position + box_radius;

    % Find the indices of the pixels associated with the lower and upper limits
    % of the box kernel.
    l_idx = zeros(size(xform_domain_position));
    u_idx = zeros(size(xform_domain_position));

    for row = 1:h
        xform_domain_pos_row = [xform_domain_position(row,:) inf];

        l_pos_row = l_pos(row,:);
        u_pos_row = u_pos(row,:);

        local_l_idx = zeros(1, w);
        local_u_idx = zeros(1, w);

        local_l_idx(1) = find(xform_domain_pos_row > l_pos_row(1), 1, 'first');
        local_u_idx(1) = find(xform_domain_pos_row > u_pos_row(1), 1, 'first');

        for col = 2:w
            local_l_idx(col) = local_l_idx(col-1) + ...
                find(xform_domain_pos_row(local_l_idx(col-1):end) > l_pos_row(col), 1, 'first') - 1;

            local_u_idx(col) = local_u_idx(col-1) + ...
                find(xform_domain_pos_row(local_u_idx(col-1):end) > u_pos_row(col), 1, 'first') - 1;
        end

        l_idx(row,:) = local_l_idx;
        u_idx(row,:) = local_u_idx;
    end

    % Compute the box filter using a summed area table.
    SAT            = zeros([h w+1 num_channels]);
    SAT(:,2:end,:) = cumsum(I, 2);
    row_indices    = repmat((1:h)', 1, w);
    F              = zeros(size(I));

    for c = 1:num_channels
        a = sub2ind(size(SAT), row_indices, l_idx, repmat(c, h, w));
        b = sub2ind(size(SAT), row_indices, u_idx, repmat(c, h, w));
        F(:,:,c) = (SAT(b) - SAT(a)) ./ (u_idx - l_idx);
    end
    
end

function T = image_transpose(I)

    [h, w, num_channels] = size(I);
    
    T = zeros([w h num_channels], class(I));
    
    for c = 1:num_channels
        T(:,:,c) = I(:,:,c)';
    end
    
end
