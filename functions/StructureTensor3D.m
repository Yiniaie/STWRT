function [Jxx, Jxy, Jxz, Jyy, Jyz, Jzz]=StructureTensor3D(ux,uy,uz)
% This function calculates the 3D Structure Tensor
% Jp ( grad(u) ) = conv ( Kp  , grad(u) * grad(u)^T )
% 
% Function is written by D.Kroon University of Twente (September 2009)

% J(grad u_sigma)
Jxx = ux.^2;
Jxy = ux.*uy;
Jxz = ux.*uz;
Jyy = uy.^2;
Jyz = uy.*uz;
Jzz = uz.^2;

% Do the gaussian smoothing

% Jxx = imgaussian(Jxx,rho,2*rho);
% Jxy = imgaussian(Jxy,rho,2*rho);
% Jxz = imgaussian(Jxz,rho,2*rho);
% Jyy = imgaussian(Jyy,rho,2*rho);
% Jyz = imgaussian(Jyz,rho,2*rho);
% Jzz = imgaussian(Jzz,rho,2*rho);
K = fspecial('gaussian', [3 3], 9); % Gaussian kernel

Jxx = imfilter(Jxx, K, 'symmetric');
Jxy = imfilter(Jxy, K, 'symmetric');
Jxz = imfilter(Jxz, K, 'symmetric');
Jyy = imfilter(Jyy, K, 'symmetric');
Jyz = imfilter(Jyz, K, 'symmetric');
Jzz = imfilter(Jzz, K, 'symmetric');