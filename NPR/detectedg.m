function I = detectedg(image, sigE, T, pieE)
        I1 = imgaussfilt(image,  sigE ,'FilterSize',5);
        I2 = imgaussfilt(image,sqrt(1.6)*sigE ,'FilterSize',5 );
        I = I1 - T*I2;
        [x,y] = size(image);
        for j = 1:y
            for i= 1:x
                if I(i,j) > 0
                    I(i,j) = 1;
                else 
                    I(i,j) = 1 + tanh(pieE*I(i,j));
                end
            end
        end
end
