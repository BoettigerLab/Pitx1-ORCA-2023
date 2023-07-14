% enter your directories here:
NAS02_v4 = 'Z:/';
data_folder = [NAS02_v4,'TzuChiao/Manuscript/Codes/data/'];
save_folder = 'U:\Manuscripts\TzuChiao Paper\NatureGenetics\Revision1\Revised_images/'; 
%% load the Erez Hi-C borders from GM
borderTable = readtable([data_folder, 'GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt'],...
    'delimiter','\t','Format','%s%u%u%s%u%u%s%f%f%f%f%f');

fname = [data_folder, 'GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_with_motifs.txt'];
opts = detectImportOptions(fname);
opts.VariableTypes{1} = 'char';  % This will change column 15 from 'double' or whatever it is to 'char', which is what you want.
opts.VariableTypes{4} = 'char';  % This will change column 15 from 'double' or whatever it is to 'char', which is what you want.
loopTable = readtable(fname,opts);

% load some Hi-C data
c = 12;
res = 25e3; 
hicStart = 44e6;
hicDomainSize = 20e6;
hicLocus = WriteLocusName(['chr',num2str(c)],hicStart,hicStart+hicDomainSize);
hoxc = ReadJuiceboxMatrix([data_folder,'hg19_GM_hoxC_balanced_25kb.txt'],'locus',hicLocus);

% load enhancer promoter pairs
% different table for each cell type
enTable = readtable([data_folder, 'fantom5_hg19_human.associations.hdr.txt'],'delimiter','\t');
[enChr,enStart,enStop] = cellfun(@ParseLocusName, enTable.enhancer,'UniformOutput',false);
prStart = enTable.distance + cat(1,enStart{:});
prName = enTable.promoter;
epTable = table(enChr, cat(1,enStart{:}), prStart, prName);
epTable.Properties.VariableNames = {'chr','enhancer_pos','promoter_pos','gene_name'};

%% a
c=12;
hicStart = 44e6;
hicDomainSize = 20e6;
cstr = num2str(c);
mk_sz = 5;
buf = 20e3; % buffer - E and P are at least this distance apart from the hopped border. 

%% Plot TADs and TAD calls
x = hicStart:res*8:hicStart+hicDomainSize;
isChrE = strcmp(epTable.chr,['chr', cstr]);
isChrB = strcmp(borderTable.chr1, cstr)   & borderTable.f1 > 1; %  
isChrL = strcmp(loopTable.chr1, cstr) & loopTable.o ./loopTable.e_donut > 0;

e = epTable.enhancer_pos(isChrE);
p = epTable.promoter_pos(isChrE);
bs = unique([borderTable.x1(isChrB); borderTable.x2(isChrB); ...
    loopTable.centroid1(isChrL); loopTable.centroid2(isChrL)]); %

figure(2); clf; imagesc(x,x, hoxc); colorbar;
GetColorMap('gray', 'flip', true); caxis([0,1400]);
figure(2); hold on;
sz = 100e3;
bsX = [bs-sz, bs];
bsY = [bs, bs];
plot((bsX+sz)', bsY','Color', 'w'); hold on;
plot(bsY', bsX','Color', 'w');
%
nEP = length(e);
epCross = nan(nEP,1);
closest_border = zeros(nEP, 2); % E, P.
% Closest borders are not the border hopped
for n=1:nEP
    % Calculate closest border for E and P
    [M, e_I] = min(abs(e(n) - bs));
    [M, p_I] = min(abs(p(n) - bs));
    closest_border(n, 1) = bs(e_I);
    closest_border(n, 2) = bs(p_I);
    % See if there are more borders in between.
    border = ((e(n) < (bs -buf)) & (p(n) > (bs + buf))) | ((e(n) > (bs+buf)) & (p(n) < (bs -buf)));
    if ~any(border)
        epCross(n) = false;
    else
        border(e_I) = false;
        border(p_I) = false;
        if any(border)
            epCross(n) = true;
        end
    end
end

e = e(~isnan(epCross));
p = p(~isnan(epCross));
epCross = epCross(~isnan(epCross));

epCross = logical(epCross);
% add E-P pairs 
figure(2); hold on; % the genome coord version
promoterFirst = e > p;
% Plot the EP-pairs--------------------------
hold on; plot(e(epCross == 1 & promoterFirst),p(epCross == 1& promoterFirst),'+','markerSize',mk_sz,'Color',[237 76 155]./255);
hold on; plot(p(epCross == 1 & ~promoterFirst),e(epCross == 1 & ~promoterFirst),'+','markerSize',mk_sz,'Color',[237 76 155]./255);
hold on; plot(e(epCross == 0 & promoterFirst),p(epCross == 0 & promoterFirst),'x','markerSize',mk_sz*1.2,'Color',[226 151 193]./255);
hold on; plot(p(epCross == 0 & ~promoterFirst),e(epCross == 0 & ~promoterFirst),'x','markerSize',mk_sz*1.2,'Color',[226 151 193]./255);

axis square;

% add promoters and enhancers on the diagonal 
figure(2); hold on; % the genome coord version
p_unique = unique(p);
e_unique = unique(e);
hold on; plot(p_unique,p_unique,'o','markerSize',mk_sz,'Color',[240,98,34]./255);
hold on; plot(e_unique,e_unique,'s','markerSize',mk_sz*1.2,'Color',[64 177 197]./255);

set(gcf,'color','w');
xlim([4.7e7,4.85e7]); ylim([4.7e7,4.85e7]); % a nice with-in border region
%%
exportgraphics(gcf,strcat(save_folder,"fig5_example.eps"),...
    'BackgroundColor','white','ContentType','vector');
%% Get the EP pair only plot

% add promoters and enhancers on the diagonal 
figure(2); hold on; % the genome coord version
axis square;
p_unique = unique(p);
e_unique = unique(e);
hold on; plot(p_unique,p_unique,'o','markerSize',mk_sz,'Color',[240,98,34]./255);
hold on; plot(e_unique,e_unique,'s','markerSize',mk_sz*1.2,'Color',[64 177 197]./255);
Xs_pfirst = [p(promoterFirst) e(promoterFirst) e(promoterFirst)]';
Ys_pfirst = [p(promoterFirst) p(promoterFirst) e(promoterFirst)]';
Xs_efirst = [e(~promoterFirst) p(~promoterFirst) p(~promoterFirst)]';
Ys_efirst = [e(~promoterFirst) e(~promoterFirst) p(~promoterFirst)]';
hold on; plot(Xs_pfirst, Ys_pfirst, 'Color', "#C8C8C8", 'LineStyle', '--');
hold on; plot(Xs_efirst, Ys_efirst, 'Color', "#C8C8C8", 'LineStyle','--');

set(gcf,'color','w');
set(gca,'YDir', 'reverse');

xlim([4.7e7,4.85e7]); ylim([4.7e7,4.85e7]); % a nice with-in border region
%%
exportgraphics(gcf,strcat(save_folder,"fig5_example_EPonly.eps"),...
    'BackgroundColor','white','ContentType','vector');

%% Register enhancer-promoter pairs
distEPs = cell(23,1); %23 chromosomes
epCrosses = cell(23,1);
epToBorders = cell(23,1);

cumNEP = 0;
for c =1:23
    if c == 23
        cstr = 'X';
    else
        cstr = num2str(c);
    end
    isChrE = strcmp(epTable.chr,['chr',cstr]);
    % enStart(isChr)
    e = epTable.enhancer_pos(isChrE);
    p = epTable.promoter_pos(isChrE);
    geneName = epTable.gene_name(isChrE);
    
    nEP = length(e);
    distEP = abs(e-p);
    epCross = nan(nEP,1);
    epToBorder = nan(nEP,2);
    
    isChrB = strcmp(borderTable.chr1, cstr)   & borderTable.f1 > 1; %  
    isChrL = strcmp(loopTable.chr1, cstr) & loopTable.o ./loopTable.e_donut > 0;
    bs = unique([borderTable.x1(isChrB); borderTable.x2(isChrB); loopTable.centroid1(isChrL); loopTable.centroid2(isChrL)]); 
    bs = cast(bs, "double");
    for n=1:nEP % n = 10
        % not consider things outside the borders
        if e(n) < bs(1) || e(n) > bs(end) || p(n) < bs(1) || p(n) > bs(end)
            continue
        end
        
        % The boundary they crossed must not be their nearest borders.
        % If no border --> epCross(n) = false
        % if yes border:
        
        % Calculate closest border for E and P
        [e_M, e_I] = min(abs(e(n) - bs));
        [p_M, p_I] = min(abs(p(n) - bs));
        closest_border(n, 1) = bs(e_I);
        closest_border(n, 2) = bs(p_I);
        % See if there are more borders in between.
        border = ((e(n) < (bs -buf)) & (p(n) > (bs + buf))) | ((e(n) > (bs+buf)) & (p(n) < (bs -buf)));
        if ~any(border)
            epCross(n) = false;
        else
            border(e_I) = false;
            border(p_I) = false;
            if any(border)
                epCross(n) = true;
            end
        end
        epToBorder(n,1) = e_M;
        epToBorder(n,2) = p_M;
    end
    distEPs{c} = distEP;
    epCrosses{c} = epCross;
    epToBorders{c} = epToBorder;
end

epToBorder = cat(1,epToBorders{:});
epCross = cat(1,epCrosses{:});
distEP = cat(1,distEPs{:});

%% b
noCrossesL = epToBorder(epCross==0 & distEP>100e3,:)/1e3;
yeCrossesL = epToBorder(epCross==1 & distEP>100e3,:)/1e3;

figure(6); clf;  % max of long range
histogram(max(yeCrossesL,[],2),0:1e1:1e3,'normalization','probability','EdgeColor','none','FaceColor',[237 76 155]./255); hold on;
histogram(max(noCrossesL,[],2),0:1e1:1e3,'normalization','probability','EdgeColor','none','FaceColor',[226 151 193]./255); hold on;
histogram(max(yeCrossesL,[],2),0:1e1:1e3,'normalization','probability','DisplayStyle','stairs', 'EdgeColor',[237 76 155]./255); 

legend({['yes cross. mean dist.=',num2str(nanmean(max(yeCrossesL,[],2)),3 ),' kb'],...
    ['no cross. mean dist.=',num2str(nanmean(max(noCrossesL,[],2)),3 ),' kb']}); 
title('max EP distance to a border');
set(gcf,'color','w'); 
xlabel('distance (kb)');
%%
noCrosses = epToBorder(epCross==0,:)/1e3;
yeCrosses = epToBorder(epCross==1,:)/1e3;

figure(7); clf;  % max of long range
histogram(max(yeCrosses,[],2),0:1e1:1e3,'normalization','probability','EdgeColor','none','FaceColor',[237 76 155]./255); hold on;
histogram(max(noCrosses,[],2),0:1e1:1e3,'normalization','probability','EdgeColor','none','FaceColor',[226 151 193]./255); hold on;
histogram(max(yeCrosses,[],2),0:1e1:1e3,'normalization','probability','DisplayStyle','stairs', 'EdgeColor',[237 76 155]./255); 

legend({['yes cross. mean dist.=',num2str(nanmean(max(yeCrosses,[],2)),3 ),' kb'],...
    ['no cross. mean dist.=',num2str(nanmean(max(noCrosses,[],2)),3 ),' kb']}); 
title('max EP distance to a border');
set(gcf,'color','w'); 
xlabel('distance (kb)');
%%
exportgraphics(gcf,strcat(save_folder,"fig5_histo.eps"),...
    'BackgroundColor','white','ContentType','vector');
%%
[h,p,ks2stat] = kstest2(max(yeCrossesL,[],2), max(noCrossesL, [], 2));
