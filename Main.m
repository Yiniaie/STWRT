clear ;
close all;
clc; 
 
%%  addpath(genpath('../InfraredLib-too;box'));       %rmpath(genpath(cd)); 
rmpath(genpath(cd));
addpath('functions')
addpath('functions/proximal_operators/')
addpath('functions\tensor_tools\')
addpath('functions\norm\')
opts.model = eval(['@','STWRT']); method_path='STWRT';frame_num=30; opts.method_name='Ours'; opts.method_time='2024';
save_path =  './result';

SamplePath1 = 'pictures';      %�洢ͼ���·��
SamplePath1=char(SamplePath1);
files = dir(fullfile(SamplePath1));
files = files(3:end);
len1 = size(files,1);%�ļ���Ŀ
files =natsortfiles(files);%��������
img_multi=[];
opts.frame_num = len1;
f = 0;
for i=1:len1  
    f = f+1;
    fileName = strcat(SamplePath1,'/',files(i).name); 
    img = imread(fileName);
    %  �����Ƿ�Ϊ����ͼ�� 
    img = Judge_infrared_image(img);
    %Ԥ���� 
%         img = imresize(img,[256,256]);
    result(f).image =img;   
    img = double(img);
    img_multi= cat(3,img_multi,img);
end

%   execuate
tic
[B_out, T_out, W_out] = STWRT(img_multi,opts);
times1=toc;
disp(['bath������ʱ�䣺',num2str(times1)]);
times2 = times1/len1;
disp(['Per����ʱ�䣺',num2str(times2)]);

for m = 1:f 
    T_out0 = uint8( normalize2(T_out(:,:,m)));
    B_out0 = uint8( normalize2(B_out(:,:,m)));
    W_out0 = uint8( normalize2(W_out(:,:,m)));
     % �����ļ���
    names_T = sprintf('result/T_%d.png', m);
    names_B = sprintf('result/B_%d.png', m);
    names_W = sprintf('result/W_%d.png', m);
    
    % ����ͼ��
    imwrite(T_out0, names_T);
    imwrite(B_out0, names_B);
    imwrite(W_out0, names_W);
end
%%    save result and opts


% end  %% Hyperparamet
% rmpath(genpath(cd)); 

