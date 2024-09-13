function ST=Prior_ST

ST.D2 = @D2;
ST.D3 = @D3;
ST.ECA_STT = @ECA_STT;


end

%%         the Original_prior structure 
function output2 = D2(input)
    output2 =[];
    for i=1:size(input,3)
        %      %step 1: calculate two eigenvalues from structure tensor
        [lambda1, lambda2] = structure_tensor_lambda(input(:,:,i), 3);
        %      %step 2: calculate corner strength function
        cornerStrength = (((lambda1.*lambda2)./(lambda1 + lambda2)));
        %      %step 3: obtain final weight map
        maxValue = (max(lambda1,lambda2));
        priorWeight0 = cornerStrength .* maxValue;
        priorWeight0(isnan(priorWeight0)) = 0;
        priorWeight0(isinf(priorWeight0)) = 0;
        priorWeight1 = mat2gray(priorWeight0);
        output2 = cat(3,output2,priorWeight1);
    end 
end


function output = D3(input)
    output=[];
    [n1, n2, total_images]  = size(input);
    num_frames =6;
    
    % 计算最后一批图像的数量
    last_batch_size = mod(total_images, num_frames);

    % 计算需要向前取几帧来凑足 6 帧
    num_frames_to_forward = num_frames - last_batch_size;
    % 循环处理图像
    for i = 1:num_frames:total_images
        % 获取当前批次的图像文件
        batch_files = input(:,:,i:min(i+num_frames-1, total_images));

        % 如果最后一批图像不足 6 帧，则向前取几帧来凑足 6 帧
        if size(batch_files,3) < num_frames
            forward_start_index = max(1, i - num_frames_to_forward);
            batch_files_input = cat(3,input(:,:,forward_start_index:i-1), batch_files);
            [batch_files_output] = D3_sub(batch_files_input);
            output = cat(3,output,batch_files_output(:,:,end-last_batch_size+1:end));
        else
            [batch_files_output] = D3_sub(batch_files);
            output = cat(3,output,batch_files_output);
           
        end
    end
    patch_max=max(output,[],[1 2]);
    output = output./patch_max; %sorted_X = sort(D(:,:,1), 'descend');
end

function [output] = D3_sub(input0)
    output0=[];
    [~, ~, n3] = size(input0);
        for i=1:n3
            usigma(:,:,i) = preprocess_gauss(input0(:,:,i));
        end
        num_end=3;
        for j =1: n3-num_end
    %%     the 3D piror structure 
%             [ux, uy, uz] = gradient(usigma);
            ux=derivatives(usigma(:,:,j:j+num_end),'x');
            uy=derivatives(usigma(:,:,j:j+num_end),'y');
            uz=derivatives(usigma(:,:,j:j+num_end),'z');
            [Jxx, Jxy, Jxz, Jyy, Jyz, Jzz] = StructureTensor3D(ux,uy,uz);          
            [mu1,mu2,mu3,v3x,v3y,v3z,v2x,v2y,v2z,v1x,v1y,v1z]=EigenVectors3D(Jxx, Jxy, Jxz, Jyy, Jyz, Jzz); 
            D=(mu1.*mu2.*mu3)./(mu1+mu2+mu3);%.*maxValue; 
            D(isnan(D)) = 0;        %sorted_X = sort(mu1(:,:,10), 'descend');lambda1
            D(isinf(D)) = 0;        %sorted_X = sort(lambda2(:,:,3), 'descend');lambda1
            patch_max=max(D,[],[1 2]);
            D = D./patch_max; %sorted_X = sort(D(:,:,1), 'descend');
            output0(:,:,j)=D(:,:,1);
            output0(:,:,j+num_end)=D(:,:,num_end);
        end
        output=output0;
end

function output = ECA_STT(input)
%%     the 3D piror structure 
        output =[];
        beta = 0.1; 
        for i=1:size(input,3)
            %1   Original
            [lambda11, lambda22] = structure_tensor_lambda(input(:,:,i),3);
            reCI = (lambda11+lambda22)./(lambda11.*lambda22);          % reciprocal of coner awareness indicator
            reEI = 1./(lambda11-lambda22);                             % reciprocal of edge awareness indicator
            WeTE = ((1+beta.^2).*reCI.*reEI)./((beta.^2.*reCI)+reEI);  % target enhancement weight
            output = cat(3,output,WeTE);
        end
end
