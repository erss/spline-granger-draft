
% clear all;
% close all;
% load('/Users/erss/Documents/MATLAB/ar_model/Paper_figures/alldataetc___/single_node_spectral_peak.mat')
fig1_single_node_low_freq 

model_tester = model_standard;
 %load('/Users/erss/Documents/MATLAB/ar_model/Paper_figures/alldataetc___/single_node_low_freq.mat')
fig2_single_node_spectral_peak
 close all
 %%
 %[x, y] =  grstat( model_true,model_spline,model_tester );
for i = 1:100
  [x, y] =  grstat1( model_true,model_spline,model_tester );
 xx(i)= x.stat;
 yy(i)= y.stat;
 fprintf(num2str(i))
  
end
barplot(Labels,yy,xx);