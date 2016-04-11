clear all
close all
clc

load phiNum.dat
load phiAna.dat

N.z = length(phiNum);
plot(log(phiNum),0:N.z-1,'b',log(phiAna),0:N.z-1,'r--');
% plot(log(phiAna),0:N.z-1,'r--');

% 
% 
% load phi.dat
% figure;
% imagesc(phi')