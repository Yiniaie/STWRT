function [I_out]=Judge_infrared_image(Img)
y3=0;
[y1,y2,y3]=size(Img);
if y3>1 %判断图像是否为红外图像
    for j1=1:y1        
        for j2=1:y2
            if (Img(j1,j2,1)==Img(j1,j2,2)) & (Img(j1,j2,1)==Img(j1,j2,3))
                
            else
                warning('不是红外图像');
                    I_out=rgb2gray(Img);
                    return;
            end
        end

    end
    I_out=rgb2gray(Img);
    return;
end
I_out =Img;
end

% 
% Img = imread('D:\Dataset\1\data1\17.BMP');
% y3=0;
% [y1,y2,y3]=size(Img);
% if y3>1 %判断图像是否为红外图像
%     for j1=1:y1        
%         for j2=1:y2
%             if (Img(j1,j2,1)==Img(j1,j2,2)) & (Img(j1,j2,1)==Img(j1,j2,3))
%                 
%             else
%                 warning('不是红外图像');
%             end
%         end
% 
%     end
% end