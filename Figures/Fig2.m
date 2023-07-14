%%
save_folder = 'U:\Manuscripts\TzuChiao Paper\NatureGenetics\Revision1\Revised_images\'; 
%% Categorize the polymers into the 3 different models
TAD1 = (pitx+1):ra3;
TAD2 = (ra3+1):ra4;
TAD3 = (ra4+1):(pen);

%% Categorizing traces. Needed for Fig. 3b, d as well.
normMap = NormMap(nanmedian(hLmaps,3));
normHLmaps = hLmaps./normMap;
% 11 and 9 are the smaller of the two TAD sizes across B1 and B2, respectively.
sep1 = getSepScores2(normHLmaps, ra3, 11, true);
sep2 = getSepScores2(normHLmaps, ra4, 9, true);
normMapF = NormMap(nanmedian(fLmaps,3));
normFLmaps = fLmaps./normMapF;
sep1f = getSepScores2(normFLmaps, ra3, 11, true);
sep2f = getSepScores2(normFLmaps, ra4, 9, true);
sepTh = 1.0;

hasData = findHasData(hLmaps);
hasDataF = findHasData(fLmaps);

isMerged = sep1<sepTh & sep2<sepTh & hasData;
isMergedF = sep1f<sepTh & sep2f<sepTh & hasDataF;

% Define stack
isStack = findStack(hLmaps, th) & ~isMerged & hasData;
isStackF = findStack(fLmaps, th) & ~isMergedF & hasDataF;

% Define other
isOther = hasData & ~isStack & ~isMerged;
isOtherF = hasDataF & ~isStackF & ~isMergedF;

%% e, f, and Extended Data Fig. 5
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
        map = hLmaps;
        ism = isMerged;
        iso = isOther;
        iss = isStack;
        has = hasData;
    else
        map = fLmaps;
        ism = isMergedF;
        iso = isOtherF;
        iss = isStackF;
        has = hasDataF;
    end
    %merge domain
    isEP = squeeze(map(pitx,pen,:)<th);
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


    %Dealing with Risk ratio
    outCellRR(1,l,:) = RelativeRisk3(isEP & has, ism, has);
    
    % stack domain
    outCellRR(2,l,:) = RelativeRisk3(isEP & has, iss, has);

    %loopOut
    outCellRR(3,l,:) = RelativeRisk3(isEP & has, iso, has);
end

% The outCell can be exported and ploted using ggplot2 in R.
% The error bars are calculated as (p * (1-p)/N)^0.5,
% where p is the proportion and N is the number.
%% Outputting outCells
fname = '20230711';

 fileID = fopen(['o3modelPercent_pooled_', fname, '.txt'],'w');
 for a = 1:size(outCell, 1)
     for l = 1:size(outCell,2)
        fprintf(fileID, "%i\t%i\t%f\t%f\t%i\n", cat(1, [a;l], squeeze(outCell(a,l,:))));
     end
 end
 fclose(fileID);

%
 fileID = fopen(['o3modelRR_pooled_', fname, '.txt'],'w');
 for a = 1:size(outCellRR, 1)
     for l = 1:size(outCellRR,2)
        fprintf(fileID, "%i\t%i\t%f\t%f\t%f\t%i\n", cat(1, [a;l], squeeze(outCellRR(a,l,:))));
     end
 end
 fclose(fileID);

%% d: individual traces
shift = -11607;
% Needs changing if total traces shifted
mergeExp = 19200 + shift;
stackedExp = 12702 + shift;
otherExp = 286;

cmap = GetColorMap('redToWhite');
    
%%
axis_dim = [0,400];
cmap(1,:) = [200,200,200]./255;

for showID = [mergeExp, stackedExp, otherExp]
    plotPolymerWrapper(hLpol(:,:,showID));
    % Separate nan from 0
    sc_map = hLmaps(:,:,showID);
    sc_map(sc_map < axis_dim(1)) = axis_dim(1) + 5;
    nanInd = find(sum(isnan(sc_map)) >= 74);
    for i = 1:size(sc_map,1)
        if ~ismember(i,nanInd)
            sc_map(i,i) = axis_dim(1) + 5;
        end
    end
    figure(showID); imagesc(sc_map); colormap(cmap);
    colorbar;axis square; clim(axis_dim);
end
%% Exporting individual plots
exportgraphics(gcf,strcat(save_folder,"2_03_OtherMap.eps"),...
   'BackgroundColor','white','ContentType','vector');

