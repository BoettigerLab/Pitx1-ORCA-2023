%%
sim_folder = 'T:\2023-05-05_Pitx1Sims\';
save_folder = 'U:\Manuscripts\TzuChiao Paper\NatureGenetics\Revision1\Revised_images\'; 

%% List out folder names, load the basic model:
    th = 5; % threshold for contact;
    
    mapStep = 5;
    hiFolders = { 
    [sim_folder, 'Model_v4\'];
    [sim_folder, '\Model_v4_sep400_hi\']; %
    [sim_folder, '\Model_v4_sep800_hi\']; %
    [sim_folder, '\Model_v4_life100_hi\']; % 1/2 life time
    [sim_folder, '\Model_v4_2111_hi\'];  % uni-directional
    [sim_folder, '\Model_v4\']};   %
    
    loFolders = { [sim_folder, '\Model_v4_lo\'];    
    [sim_folder, '\Model_v4_sep400_lo\'];
    [sim_folder, '\Model_v4_sep800_lo\'];
    [sim_folder, '\Model_v4_life100_lo\'];
    [sim_folder, '\Model_v4_2111_lo\']; 
    [sim_folder, '\Model_v4_mid\']};
    
k = 1;

polyStep =5;
[hi_pols,hi_maps,hi_LEFs,~,~] = LoadOpenPolySim(hiFolders{k},'mapStep',polyStep,'polyStep',polyStep,'timeStep',1,'mapStep',1);
[lo_pols,lo_maps,lo_LEFs,~,~] = LoadOpenPolySim(loFolders{k},'mapStep',polyStep,'polyStep',polyStep,'timeStep',1,'mapStep',1);
%% a, example polymers
ctcf = 41:40:161;
% find example with 3 way contact (hindlimb)
r = 100;
t = 24;
%
currPoly = hi_pols{r}(:,:,t);
lefs = hi_LEFs{r,t};
ctcf3D = currPoly(ctcf,:);
%
plotSimPolymer(currPoly, ctcf3D, lefs, 'thin', 'dumbbell');
%
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);
axis off;
%% Find which coordinates are linked by LEF.
lefCoords = zeros(size(lefs,1), 2);
for i = 1:size(lefs, 1)
    for t = 1:2
        lefCoords(i,t) = find(currPoly(:,1) == lefs(i,1,t) & ...
            currPoly(:,2) == lefs(i,2,t) & ...
            currPoly(:,3) == lefs(i,3,t));
    end
end
lefCoords
%% Exporting individual plots
exportgraphics(gcf,strcat(save_folder,"4_01_hi_example.eps"),...
   'BackgroundColor','white','ContentType','vector');

%% b, example lo-polymer
r = 4;
t = 26;

currPoly = lo_pols{r}(:,:,t);
lefs = lo_LEFs{r,t};
ctcf3D = currPoly(ctcf,:);

plotSimPolymer(currPoly, ctcf3D, lefs, 'thin', 'dumbbell');
%
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);
axis off;
%% Find which coordinates are linked by LEF.
lefCoords = zeros(size(lefs,1), 2);
for i = 1:size(lefs, 1)
    for t = 1:2
        lefCoords(i,t) = find(currPoly(:,1) == lefs(i,1,t) & ...
            currPoly(:,2) == lefs(i,2,t) & ...
            currPoly(:,3) == lefs(i,3,t));
    end
end
lefCoords
%% Exporting individual plots
exportgraphics(gcf,strcat(save_folder,"4_01_lo_example.eps"),...
   'BackgroundColor','white','ContentType','vector');

%% Extended Data Fig. 8: Analysis of centrality
hiPols = cat(3, hi_pols{:});
loPols = cat(3, lo_pols{:});
[hDist, hrG] = getCDistRg(hiPols);
[lDist, lrG] = getCDistRg(loPols);
hiMaps = cat(3, hi_maps{:});
loMaps = cat(3, lo_maps{:});

%% ED Fig. 8a: median # contact

hContact = squeeze(sum(hiMaps < th))';
lContact = squeeze(sum(loMaps < th))';

% Boot strapping for median estimate
boot_time = 1000;
contactMeds = zeros(boot_time, 200, 2);

% 95% interval
% 1000 times. %resampling N = N
rng(26);
inds = randi(size(hDist, 1), boot_time, size(hDist, 1));

for i = 1:boot_time
    contactMeds(i,:,1) = median(hContact(inds(i,:), :), 1, 'omitnan');
    contactMeds(i,:,2) = median(lContact(inds(i,:), :), 1, 'omitnan');
end
% median plot
contactMedSum = zeros(200, 3, 2); %[L, M, U] * [hi, lo]
Prts = [5, 50, 95];
for i = 1:3
    contactMedSum(:,i,:) = prctile(contactMeds, Prts(i));
end
% writematrix(contactMedSum, [save_folder, 'SuppFig8_MultiContact.txt']) %
% output with Rscripts

%% ED Fig. 8b: Boot strapping for median distance estimate
boot_time = 1000;
distMeds = zeros(boot_time, 200, 2);

rng(26);
inds = randi(size(hDist, 1), boot_time, size(hDist, 1));

for i = 1:boot_time
    distMeds(i,:,1) = median(hDist(inds(i,:), :), 1, 'omitnan');
    distMeds(i,:,2) = median(lDist(inds(i,:), :), 1, 'omitnan');
end
% median plot
distMedSum = zeros(200, 3, 2); %[L, M, U] * [hi, lo]
Prts = [5, 50, 95];
for i = 1:3
    distMedSum(:,i,:) = prctile(distMeds, Prts(i));
end
% writematrix(distMedSum, [save_folder, 'SuppFig7_CenterDist.txt']);

%% ED Fig. 8c: Produce graph for the scatter plot
th = 5;
goodCen = reshape(hDist,1,[]);
ContactRate = reshape(squeeze(sum(hiMaps < th))',1,[]);
% Good old scatter
figure();
dscatter(goodCen', ContactRate', 'MARKER', 's', ...
    'smoothing', 'None', 'MSIZE', 10, 'FILLED', true, ...
    'BINS', [200, length(unique(ContactRate))]);
%ylim([0 0.6]);
set(gcf,'position',[200,200,500,300]);
colorbar;

% Get the log-log Pearson correlation
[R,P, RL, RU] = corrcoef(log(goodCen),log(ContactRate'));

%%
exportgraphics(gcf,strcat(save_folder,"SuppFig8_scatter.eps"),...
   'BackgroundColor','white','ContentType','vector');


%% c-f and Extended Data Fig. 9: different parameters
K = length(hiFolders);
load('cutoffBluewhitered2');
t = 5;
for k = 1:K
%
    polyStep =5;
    [hi_pols,hi_maps,hi_LEFs,~,~] = LoadOpenPolySim(hiFolders{k},'mapStep',polyStep,'polyStep',polyStep,'timeStep',1,'mapStep',1);
    [lo_pols,lo_maps,lo_LEFs,~,~] = LoadOpenPolySim(loFolders{k},'mapStep',polyStep,'polyStep',polyStep,'timeStep',1,'mapStep',1);
    hiMaps = cat(3,hi_maps{:});
    loMaps = cat(3,lo_maps{:});
    contMapHi = ContactFrac(hiMaps,'threshold',t);
    contMapLo = ContactFrac(loMaps,'threshold',t);
    figure(1); subplot(K,2,2*(k-1)+1); 
    imagesc(log2(contMapHi)); 
    axis square; clim([-8 0]);
    set(gca,'XTick',[], 'YTick', []);
    nameparts = strsplit(hiFolders{k},filesep);
    ylabel(nameparts{end-1},'interpreter','none');
    figure(1); subplot(K,2,2*k); 
    imagesc(log2(contMapLo)); 
    axis square; clim([-8 0]);
    rm = GetColorMap('redToWhiteK');
    colormap(flipud(rm));
    set(gca,'XTick',[], 'YTick', []);


    figure(3); subplot(K,1,k);
    imagesc(log2(contMapHi./contMapLo)); clim([-1.5 1.5]);
     colormap(B);
     ylabel(nameparts{end-1},'interpreter','none');
     axis square;
     set(gca,'XTick',[], 'YTick', []);
% % 
    ctcfSites = [41, 81, 121, 161];
    norm_hiMaps = hiMaps./NormMap(nanmedian(hiMaps,3));

    sep1 = getSepScores2(norm_hiMaps, 81, 40, true);
    sep2 = getSepScores2(norm_hiMaps, 121, 40, true);
    
    isMerge_h = sep1 < 1 & sep2 < 1;
    isStack_h = squeeze(hiMaps(ctcfSites(1),ctcfSites(2),:) < t | ...
        hiMaps(ctcfSites(2),ctcfSites(3),:) < t | ...
        hiMaps(ctcfSites(3),ctcfSites(4),:) < t | ...
        hiMaps(ctcfSites(2),ctcfSites(4),:) < t  );  
    isStack_h = isStack_h & ~isMerge_h;
    isLoop_h = ~(isStack_h | isMerge_h);
    isEP_h = squeeze(hiMaps(ctcfSites(1),ctcfSites(4),:) < t);

    rr_merge_h = RelativeRisk3(isEP_h,isMerge_h,ones(size(isEP_h)));
    rr_stack_h = RelativeRisk3(isEP_h,isStack_h,ones(size(isEP_h)));
    rr_loop_h = RelativeRisk3(isEP_h,isLoop_h,ones(size(isEP_h)));


    norm_loMaps = loMaps./NormMap(nanmedian(loMaps,3));
    sep1 = getSepScores2(norm_loMaps, 81, 40, true);
    sep2 = getSepScores2(norm_loMaps, 121, 40, true);
    isMerge_f = sep1 < 1 & sep2 < 1;
    isStack_f = squeeze(hiMaps(ctcfSites(1),ctcfSites(2),:) < t | ...
        hiMaps(ctcfSites(2),ctcfSites(3),:) < t | ...
        hiMaps(ctcfSites(3),ctcfSites(4),:) < t | ...
        hiMaps(ctcfSites(2),ctcfSites(4),:) < t  ); 
    isStack_f = isStack_f & ~isMerge_f;
    isLoop_f = ~(isStack_f | isMerge_f);
    isEP_f = squeeze(hiMaps(ctcfSites(1),ctcfSites(4),:) < t);

    rr_merge_f = RelativeRisk3(isEP_f, isMerge_f, ones(size(isEP_f)));
    rr_stack_f = RelativeRisk3(isEP_f, isStack_f, ones(size(isEP_f)));
    rr_loop_f = RelativeRisk3(isEP_f, isLoop_f, ones(size(isEP_f)));


    % c, d: cell proportions
    outCell = zeros(4,2,3); %
    % The outCell contains such:
    % 1st dimension: [merge, stack, loop, any]
    % 2nd dimension: [hindlimb, forelimb]
    % 3rd dimension: [proportion of all chromatin fitting a model, 
    %  proportion of chromatin where Pen-Pitx1 is contacting and fitting a model, 
    %  number of all chromatin fitting a model];

    % 3 model RR pooled
    outCellRR = zeros(3,2,4);
    % The outCell contains such:
    % 1st dimension: [merge, stack, loop]
    % 2nd dimension: [hindlimb, forelimb]
    % 3rd dimension: [Risk ratio, confidence interval left end, 
    %  confidence interval right end, number of all chromatin fitting a model];

    for l = 1:2
        if l == 1
            isEP = isEP_h;
            ism = isMerge_h;
            iso = isLoop_h;
            iss = isStack_h;
            has = ones(size(isMerge_h));

            %Dealing with Risk ratio
            outCellRR(1,l,:) = rr_merge_h;
            outCellRR(2,l,:) = rr_stack_h;
            outCellRR(3,l,:) = rr_loop_h;
        else
            isEP = isEP_f;
            ism = isMerge_f;
            iso = isLoop_f;
            iss = isStack_f;
            has = ones(size(isMerge_f));

            %Dealing with Risk ratio
            outCellRR(1,l,:) = rr_merge_f;
            outCellRR(2,l,:) = rr_stack_f;
            outCellRR(3,l,:) = rr_loop_f;
        end

        %merge domain
        outCell(1,l,1) = sum(ism)/sum(has);
        outCell(1,l,2) = sum(isEP & ism)/sum(has);
        outCell(1,l,3) = sum(has);

        outCell(4,l,1) = sum(isEP & has) / sum(has);
        outCell(4,l,2) = sum(isEP & has) / sum(has);
        outCell(4,l,3) = sum(has);    

        % stack domain
        outCell(2,l,1) = sum(iss)/sum(has);
        outCell(2,l,2) = sum(isEP & iss)/sum(has);
        outCell(2,l,3) = sum(has);

        %other
        outCell(3,l,1) = sum(iso)/sum(has);
        outCell(3,l,2) = sum(isEP & iso)/sum(has);
        outCell(3,l,3) = sum(has);

    end
 
% The outCell can be exported and ploted using ggplot2 in R.
% The error bars are calculated as (p * (1-p)/N)^0.5,
% where p is the proportion and N is the number.

    fname = [num2str(k), nameparts{end-1}];

     fileID = fopen([save_folder,'o3modelPercent_pooled_', fname, '.txt'],'w');
     for a = 1:size(outCell, 1)
         for l = 1:size(outCell,2)
            fprintf(fileID, "%i\t%i\t%f\t%f\t%i\n", cat(1, [a;l], squeeze(outCell(a,l,:))));
         end
     end
     fclose(fileID);

    %
     fileID = fopen([save_folder,'o3modelRR_pooled_', fname, '.txt'],'w');
     for a = 1:size(outCellRR, 1)
         for l = 1:size(outCellRR,2)
            fprintf(fileID, "%i\t%i\t%f\t%f\t%f\t%i\n", cat(1, [a;l], squeeze(outCellRR(a,l,:))));
         end
     end
     fclose(fileID);
end
%% Exporting individual plots
exportgraphics(gcf,strcat(save_folder,"4_02_contactMaps.eps"),...
   'BackgroundColor','white','ContentType','vector');
%% Exporting individual plots
exportgraphics(gcf,strcat(save_folder,"4_02_diffs.eps"),...
   'BackgroundColor','white','ContentType','vector');
