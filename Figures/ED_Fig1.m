% Code for "" Supplementary Fig.1 
% Please run "Fig0_Load_data.m" first.
% enter your directories here:
data_folder = 'Z:/TzuChiao/Manuscript/Codes/data/';
save_folder = "U:\Manuscripts\TzuChiao Paper\NatureGenetics\Revision1\Revised_images/"; 

%% a: hi-C and contact freq. from ORCA
% Load cHi-C data
hL_chic = load([data_folder, 'table_mm10HL-E115-Wt-Mm-cHiC-merged.hicup.MAPQ30.KR_5kb.WashU_5kb.mat']);
hL_chic = hL_chic.mat;
hL_chic = imresize(hL_chic(1:150,1:150), [75,75], "bilinear");
% Produce plot
figure(); imagesc(hL_chic);
rm = GetColorMap('redToWhiteK');
colormap(flipud(rm));
clim([25 250]);
axis square;

saveas(gcf, strcat(save_folder, "fig1_hL_hiC_K.epsc"));


fL_chic = load([data_folder, 'table_mm10FL-E115-Wt-Mm-cHiC-merged.hicup.MAPQ30.KR_5kb.WashU_5kb.mat']);
fL_chic = fL_chic.mat;
fL_chic = imresize(fL_chic(1:150,1:150), [75,75], "bilinear");
% Produce plot
figure(2); imagesc(fL_chic);
rm = GetColorMap('redToWhiteK');
colormap(flipud(rm));
clim([25 250]);
axis square;

saveas(gcf, strcat(save_folder, "fig1_fL_hiC_K.epsc"));

%% b,c: median distance comparison  
ca = [300,500];
rm = GetColorMap('redToWhite');
im1 = nanmedian(hLmaps,3);
im1 = InterpMapNans(im1,'badHybes',badHybes,'badPixels',[12,24]);
im2 = nanmedian(fLmaps,3);
im2 = InterpMapNans(im2,'badHybes',badHybes,'badPixels',[12,24]);
f4 = figure(4); clf;
s1 = subplot(1,2,1); imagesc(im1 );caxis(ca);  colorbar; axis square; colormap(rm);
title(['HL n=',num2str(size(hLmaps,3))]);
s2 = subplot(1,2,2); imagesc( im2 ); caxis(ca);  colorbar; axis square; colormap(s2, rm);
title(['FL n=',num2str(size(fLmaps,3))]);
saveas(gcf, strcat(save_folder, "Suppfig1_HF_medianDists.epsc"));
%% Difference map and statistical test (Takes a while)
allL = cat(3, hLmaps, fLmaps);
[subMap, pMap, rawSubMap] = pTresholdSubMap(allL, 1:size(hLmaps,3), (size(hLmaps,3)+1):size(allL,3), ...
           'normalization', 'none', 'type', 'dist','plotOption','', 'correction', 'none');
%% d: raw distance difference
load cutoffBluewhitered2.mat B;
interpSubMap = InterpMapNans(rawSubMap,'badHybes',badHybes,'badPixels',[12,24]);
figure();   imagesc(interpSubMap);    colormap(B);    colorbar(); axis square;
caxis([-40,40]);
title(strcat("HL vs FL median distance, raw, n= ", ...
                num2str([size(hLmaps,3), size(fLmaps,3)])));
saveas(gcf, strcat(save_folder, "SFig1_raw.epsc"));
%% e: pvalues of distance differences
cmaps = makeCMAPS();
filteredPMap = pMap;
filteredPMap(badHybes, :) = 1;
filteredPMap(:, badHybes) = 1;

% Create truncated cmap
cmap = cmaps{5};
cmap = [cmap;[200,200,200]./255];
figure(); 
imagesc(filteredPMap);
colormap(cmap); colorbar; axis square;
caxis([-10 0]);

saveas(gcf, strcat(save_folder, "SFig1_p.epsc"));
