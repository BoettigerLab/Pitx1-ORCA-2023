%% Fig S3: Group individual experiment
% Enter where you want to save the figures.
save_folder = 'U:\Manuscripts\TzuChiao Paper\NatureGenetics\Revision1\Revised_images/'; 
%%
individual_contact = cell(4,2);
individual_distance = cell(4,2);
for i = 1:3
    for j = 1:2
        [contact,~] = ContactFrac(limbDists{i,j},'threshold',th);
        individual_contact{i,j} = InterpMapNans(contact,'badHybes',badHybes,'badPixels',[12,24]);
        individual_distance{i,j} = InterpMapNans(median(limbDists{i,j}, 3, 'omitnan'),'badHybes',badHybes,'badPixels',[12,24]);
    end
end
for j = 1:2
    comb4 = cat(3, limbDists{4,j}, limbDists{5,j});
    [contact,~] = ContactFrac(comb4,'threshold',th);
    individual_contact{4,j} = InterpMapNans(contact,'badHybes',badHybes,'badPixels',[12,24]);
    individual_distance{4,j} = InterpMapNans(median(comb4, 3, 'omitnan'),'badHybes',badHybes,'badPixels',[12,24]);
end
%%
tiledlayout(4,2);
for i = 1:4
    for j = 1:2
        nexttile;
        map = individual_contact{i,j};
        imagesc(map);
        rm = GetColorMap('redToWhiteK');
        colormap(flipud(rm));
        clim([0.1 0.4]);
        % colorbar;
        axis square;
    end
end
saveas(gcf, strcat(save_folder, "SuppFig4_contacts.epsc"));

%% Correlation among experiments

corrD = zeros(0, 8);
corrC = zeros(0, 8);
for d = 1:75
    dVecH = [];
    cVecH = [];
    dVecF = [];
    cVecF = [];
    for i = 1:4
        cVecH = horzcat(cVecH, diag(individual_contact{i,1}, d));
        cVecF = horzcat(cVecF, diag(individual_contact{i,2}, d));
        dVecH = horzcat(dVecH, diag(individual_distance{i,1}, d));
        dVecF = horzcat(dVecF, diag(individual_distance{i,2}, d));
    end
    corrC = cat(1, corrC, horzcat(cVecH, cVecF));
    corrD = cat(1, corrD, horzcat(dVecH, dVecF));
end
%
rm = GetColorMap('redToWhite');
[R1,P1] = corrcoef(corrC, 'Rows', 'pairwise');
[R2,P2] = corrcoef(corrD, 'Rows', 'pairwise');
% plotting
figure();
imagesc(R1); axis square; colormap(flipud(rm));
clim([0.9 1]);
saveas(gcf, strcat(save_folder, "SuppFig4_contactCorr.epsc"));
figure();
imagesc(R2); axis square; colormap(flipud(rm)); clim([0.7 1]);
saveas(gcf, strcat(save_folder, "SuppFig4_distCorr.epsc"));
%
figure();
HL_ind = 1:4;
FL_ind = 5:8;
imagesc(R1(HL_ind,HL_ind)); axis square; colormap(flipud(rm));
clim([0.9 1]);  colorbar;
saveas(gcf, strcat(save_folder, "SuppFig4_cHCorr.epsc"));
figure();
imagesc(R1(FL_ind,FL_ind)); axis square; colormap(flipud(rm));
clim([0.9 1]);  colorbar;
saveas(gcf, strcat(save_folder, "SuppFig4_cFCorr.epsc"));


