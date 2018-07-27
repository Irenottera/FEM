clc 
clear all 


eL2 = [ 0.0171 ; 0.0079 ; 0.005 ; 0.0037 ; 0.0029 ; 0.0023 ; 0.0019; 0.0011; 8.4593e-04];

totelem = [4*4 ; 8*8 ; 12*12 ; 16*16 ; 20*20 ; 25*25 ; 30*30; 50*50 ; 64*64];
h = [1/4 ; 1/8 ; 1/12 ; 1/16 ; 1/20 ; 1/25 ; 1/30; 1/50; 1/64];
h2 = h.^2;
% loglog (h, eL2)
% xlim([0,0.3]);
% ylim([0,0.3]);
loglog(totelem, eL2)
hold on 
loglog (totelem, h)
legend('error' , 'h')
xlabel ('number of elements')


