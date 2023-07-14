% scanning the threshold
data_folder = 'Z:/TzuChiao/Manuscript/Codes/data/';
save_folder = 'U:\Manuscripts\TzuChiao Paper\NatureGenetics\Revision1\Revised_images\';
%% Determine the 3 models using Fig2's code.

%% Thresholds
thresholds = [25, 50, 100,150, 200:100:400, 600, 800, 1000];

% What to do for each threshold?
% 1. make contact freq map
% a. HL, FL, dif.
% b. that contact freq diff at specific locus? no?
% c. correlation with hi-C?
% 2. Analyze the 3 prevalences
% What to record for the 3 prevalences?
% 1. HL-FL for 3 models
% 2. % EP belonging to 3 models
% 3. RR for 3 models.
for th = thresholds

    [hLcontact,~] = ContactFrac(hLmaps,'threshold',th);
    hLcontact = InterpMapNans(hLcontact,'badHybes',badHybes,'badPixels',[12,24]);
    %
    % f1 = figure(1); clf;
    % low_c = prctile(hLcontact, 18, "all");
    % high_c = prctile(hLcontact, 94, "all");
    % imagesc(hLcontact); colorbar; caxis([low_c, high_c]); 
    % rm = flipud(GetColorMap('redToWhiteK'));
    % axis square;
    % colormap(rm);
    % set(gcf,'color','w');
    % saveas(f1, strcat(save_folder, num2str(th, '%04.f'), "_hLmap.epsc"));
    % 
    %
    [fLcontact,n2] = ContactFrac( fLmaps,'threshold',th);
    fLcontact = InterpMapNans(fLcontact,'badHybes',badHybes,'badPixels',[12,24]);
    hLcontact(hLcontact<0) = 0; fLcontact(fLcontact<0)=0;

    % 
    % % record key loci diffs
    % ens = [pel, ra3, ra4, pen, neurog];
    % enNames = {'pel','ra3','ra4','pen','neurog'};
    % 
    % en_dist_out = zeros(4,5);%[hL_EP, hL_N, fL_EP, fL_N; pelB, RA3, RA4, Pen, Neurog1];
    % en_dist_out(1,:) = sum(squeeze(hLmaps(ens, pitx, :) < th), 2);
    % en_dist_out(2,:) = sum(squeeze(hLmaps(ens, pitx, :) < inf), 2);
    % en_dist_out(3,:) = sum(squeeze(fLmaps(ens, pitx, :) < th), 2);
    % en_dist_out(4,:) = sum(squeeze(fLmaps(ens, pitx, :) < inf), 2);
    % 
    % 
    % fileID = fopen(strcat(save_folder, num2str(th, '%04.f'), "_enhancerDiffs.txt"),'w');
    % for a = 1:size(en_dist_out, 1)
    %     fprintf(fileID, "%f\t%f\t%f\t%f\t%f\r\n", en_dist_out(a,:));
    % end
    % fclose(fileID);
    % 
    % 
    % % Dif map
    % f2 = figure(3); clf; 
    % imagesc(log2(hLcontact./fLcontact)); colorbar;
    % %title('log_2( hindlimb / forelimb ) contact freq.');
    % set(gcf,'color','w'); 
    % axis image;  
    % load cutoffBluewhitered2.mat B;
    % colormap(B);
    % caxis([-.3 .3]);
    % saveas(f2, strcat(save_folder, num2str(th, '%04.f'), "_hfDiffmap.png"));
    % saveas(f2, strcat(save_folder, num2str(th, '%04.f'), "_hfDiffmap.epsc"));
    % 
    % Correlation with hi-C
    hL_chic = load([data_folder, 'table_mm10HL-E115-Wt-Mm-cHiC-merged.hicup.MAPQ30.KR_5kb.WashU_5kb.mat']);
    hL_chic = hL_chic.mat;
    fL_chic = load([data_folder, 'table_mm10FL-E115-Wt-Mm-cHiC-merged.hicup.MAPQ30.KR_5kb.WashU_5kb.mat']);
    fL_chic = fL_chic.mat;
    chic_FL = imresize(fL_chic(1:150,1:150), [75,75], "nearest");
    chic_HL = imresize(hL_chic(1:150,1:150), [75,75], "nearest");
    corrS = zeros(0,4); %chicF, chicH, contF, contH
    for d = 1:75
        dVec = horzcat(diag(chic_FL, d), diag(chic_HL, d), diag(fLcontact, d), diag(hLcontact, d));
        corrS = vertcat(corrS, dVec);
    end
    [hL_corr,hL_p] = corr(log10(corrS(:,2)+eps),log10(corrS(:,4)+eps));
    [fL_corr,fL_p] = corr(log10(corrS(:,1)+eps),log10(corrS(:,3)+eps));
    %
    fileID = fopen(strcat(save_folder, num2str(th, '%04.f'), "_hicCorr.txt"),'w');
    fprintf(fileID, "%f\t%f\t%f\t%f", [hL_corr,hL_p, fL_corr,fL_p]);
    fclose(fileID);

    continue
    
    % Categorize the polymers into the 3 different models
    isStack = findStack(hLmaps, th) & ~isMerged & hasData;
    isStackF = findStack(fLmaps, th) & ~isMergedF & hasDataF;
    
    % Define other
    isOther = hasData & ~isStack & ~isMerged;
    isOtherF = hasDataF & ~isStackF & ~isMergedF;
    
    %
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


    
    fileID = fopen(strcat(save_folder, num2str(th, '%04.f'), "_3modelPercent.txt"),'w');
    for a = 1:size(outCell, 1)
        for l = 1:size(outCell,2)
            fprintf(fileID, "%i\t%i\t%f\t%f\t%i\r\n", cat(1, [a;l], squeeze(outCell(a,l,:))));
        end
    end
    fclose(fileID);
    
   
    % The outCell can be exported and ploted using ggplot2 in R.
    %
    fileID = fopen(strcat(save_folder, num2str(th, '%04.f'), "_RR.txt"),'w');
    for a = 1:size(outCellRR, 1)
        for l = 1:size(outCellRR,2)
            fprintf(fileID, "%i\t%i\t%f\t%f\t%f\t%i\r\n", cat(1, [a;l], squeeze(outCellRR(a,l,:))));
        end
    end
    fclose(fileID);
end
