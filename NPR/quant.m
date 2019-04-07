%image -- matrix values of image in the range of 1 to 100
%q_level -- quantization level
%I     -- 
function I = quant(image, q_level)
        [x, y] = size(image);
        I = zeros(x,y);
        for j = 1:y
            for i = 1:x
                I(i,j) = converter((image(i,j)))*q_level;
            end
        end
end

function val = converter(x)
        val = int16(x/10);
end