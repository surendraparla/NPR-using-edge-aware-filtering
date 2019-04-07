%  Domain transform interpolated convolution (IC) edge-preserving filter.
% 
%  F = IC(image, sigma_s, sigma_r, num_iter)
%
%  Parameters:
%    image           Input image to be filtered.
%    sigma_s         Filter spatial standard deviation.
%    sigma_r         Filter range standard deviation.
%    num_iter  Number of iterations to perform (default: 3).


function F = IC(image, sigma_s, sigma_r, num_iter)

    I = double(image);

    if ~exist('num_iter', 'var') %if variable doesnot exist 
        num_iter = 3;
    end
    [h, w, num_channels] = size(I);

    % Compute the domain transform 

    dIcdx = diff(I, 1, 2);     % partial differentiation along y
    dIcdy = diff(I, 1, 1);       % partial differentiation along x
    
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
    
    F = cast(F, class(img));

end

% Box filter interpolated convolution in the transformed domain.
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
    areas = bsxfun(@times, ...
                0.5 .* (I(:,2:end,:) + I(:,1:end-1,:)), ...
                xform_domain_position(:,2:end,:) - xform_domain_position(:,1:end-1,:) ...
            );
    
    SAT            = zeros([h w num_channels]);
    SAT(:,2:end,:) = cumsum(areas, 2);
    row_indices    = repmat((1:h)', 1, w);
    F              = zeros(size(I));

    I                     = padarray(I,                     [0 1 0], 'replicate');
    SAT                   = padarray(SAT,                   [0 1 0]);
    xform_domain_position = padarray(xform_domain_position, [0 1 0], 'replicate');
    
    % Pixel values outside the bounds of the image are assumed to equal the
    % nearest pixel border value.
    xform_domain_position(:,1)   = xform_domain_position(:,1)   - 1.2 * box_radius;
    xform_domain_position(:,end) = xform_domain_position(:,end) + 1.2 * box_radius;
    
    l_idx = l_idx + 1;
    
    for c = 1:num_channels
        
        l1_c = sub2ind(size(SAT), row_indices, l_idx, repmat(c, h, w));
        u0_c = sub2ind(size(SAT), row_indices, u_idx, repmat(c, h, w));
        
        l0_c = sub2ind(size(SAT), row_indices, l_idx - 1, repmat(c, h, w));
        u1_c = sub2ind(size(SAT), row_indices, u_idx + 1, repmat(c, h, w));
        
        l1 = sub2ind(size(SAT), row_indices, l_idx);
        u0 = sub2ind(size(SAT), row_indices, u_idx);
        
        l0 = sub2ind(size(SAT), row_indices, l_idx - 1);
        u1 = sub2ind(size(SAT), row_indices, u_idx + 1);

        % Full (center) areas.
        C = SAT(u0_c) - SAT(l1_c);
        
        % Left fractional areas.
        alpha = (l_pos - xform_domain_position(l0)) ./ (xform_domain_position(l1) - xform_domain_position(l0));
        yi    = I(l0_c) + alpha .* ( I(l1_c) - I(l0_c) );
        L     = 0.5 .* (yi + I(l1_c)) .* (1-alpha) .* (xform_domain_position(l1) - xform_domain_position(l0));
        
        % Right fractional areas.
        alpha = (u_pos - xform_domain_position(u0)) ./ (xform_domain_position(u1) - xform_domain_position(u0));
        yi    = I(u0_c) + alpha .* ( I(u1_c) - I(u0_c) );
        R     = 0.5 .* (yi + I(u0_c)) .* (alpha) .* (xform_domain_position(u1) - xform_domain_position(u0));
        
        F(:,:,c) = (L + C + R) / (2 * box_radius);

    end
    
end

function T = image_transpose(I)

    [h, w, num_channels] = size(I);
    
    T = zeros([w h num_channels], class(I));
    
    for c = 1:num_channels
        T(:,:,c) = I(:,:,c)';
    end
    
end
