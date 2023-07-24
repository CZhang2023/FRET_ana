% Step response analysis to get Rd, Ru and Radapt.
%
i = 1;
plot(t,R,'ko','markersize',4,'markerfacecolor','k')
[~,y] = ginput(2);RLdown(i) = min(y);RLup(i) = max(y);clear y;
[x,~] = ginput(2);Radapt(i) = mean(R(t>min(x)&t<max(x)));clear x;
% L = L/1000;
clc;
clear i;
%}

% Get a0 via saturated stimulus.
%
clear all;
files = dir('R_*.mat');
a_all = [];
L_all = [];
for index_file = 1:length(files)
    load(files(index_file).name);
    for index_L = 1:length(L)
        if L(index_L)>=0.05
            temp = (1-RLdown(index_L))/(RLup(index_L)-RLdown(index_L));
            a_all = [a_all;temp];
            L_all = [L_all;L(index_L)];
        end
    end
clearvars -except *all files index*
end
mean(a_all)
min(a_all)
%}

% Transform Radapt to Aadapt.
%
clear all;
load('DR_HCB1288_all.mat');
files = dir('R_*.mat');
Aadapt_all = [];
Ladapt_all = [];
for index_file = 1:length(files)
    load(files(index_file).name);
    for index_L = 1:length(L)
        Aadapt = (Radapt(index_L) - RLdown(index_L))/(1-RLdown(index_L))...
        *(DRfitForA0(0)-DRfitForA0(L(index_L)))+DRfitForA0(L(index_L));
        Aadapt_all = [Aadapt_all;Aadapt];
        Ladapt_all = [Ladapt_all;L(index_L)];
    end
clearvars -except *all files index* DRfitForA0;
end
clearvars -except Aadapt_all Ladapt_all DRfitForA0;
hold on;
subplot(221);semilogx(Ladapt_all,Aadapt_all,'bo');

%}

% the average Aadapt for different L.
%
Ladapt_ave = unique(Ladapt_all);
for i = 1:length(Ladapt_ave)
    Aadapt_ave(i,1) = mean(Aadapt_all(Ladapt_all==Ladapt_ave(i)));
    Aadapt_std(i,1) = std(Aadapt_all(Ladapt_all==Ladapt_ave(i)),1);
end

EL_ave = log((1+Ladapt_ave./0.0182)./(1+Ladapt_ave./3));
EM_ave = (log(1./Aadapt_ave-1))/DRfitForA0.N-EL_ave;
M_ave = -EM_ave/1.875+1;
Kr_b_ave = Aadapt_ave./(1-Aadapt_ave);

Kr_b_all = Aadapt_all./(1-Aadapt_all);
Kr_b_all3 = Aadapt_all.^3./(1-Aadapt_all);
for i = 1:length(Ladapt_ave)
    Kr_b_all_ave(i) = mean(Kr_b_all(Ladapt_all==Ladapt_ave(i)));
    Kr_b_all3_ave(i) = mean(Kr_b_all3(Ladapt_all==Ladapt_ave(i)));
    Kr_b_all_std(i) = std(Kr_b_all(Ladapt_all==Ladapt_ave(i)),1);
    Kr_b_all3_std(i) = std(Kr_b_all3(Ladapt_all==Ladapt_ave(i)),1);
end

subplot(223);plot(EL_ave,EM_ave,'bo');
subplot(224);plot(M_ave,Kr_b_ave,'bo');
hold off;
errorbar(Ladapt_ave,Aadapt_ave,Aadapt_std);
%}

%
EL_all = log((1+Ladapt_all./0.0182)./(1+Ladapt_all./3));
EM_all = (log(1./Aadapt_all-1))/DRfitForA0.N-EL_all;
M_all = -EM_all/1.875+1;
Kr_b_all = Aadapt_all./(1-Aadapt_all);

a*exp(-3.*(1.095*x+2.4)^2/2/11.4/l)
%}
