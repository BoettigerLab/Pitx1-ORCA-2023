function outCell = RelativeRisk2(condition, exposure, varargin)
%  [rr,rr_CI] = RelativeRisk(condition,exposure)
%  [rr,rr_CI] = RelativeRisk(condition,exposure,'cI',.95)
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
pars = inputParser;
addParameter(pars, 'cI', .95);
addParameter(pars, 'oneLimbPolys', {});
addParameter(pars, 'lengthArray', []);
parse(pars, varargin{:});
% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
cI = pars.Results.cI;
oneLimbPolys = pars.Results.oneLimbPolys;
lengthArray = pars.Results.lengthArray;

Z = norminv([(1 - cI)/2, 0.5 + 0.5*cI]);
Z = Z(2);


if ~isempty(oneLimbPolys)
    backCondition = mapBackSelection(condition, oneLimbPolys);
    backExposure = mapBackSelection(exposure, oneLimbPolys);
elseif ~isempty(lengthArray)
    backCondition = mapBackSelection2(condition, lengthArray);
    backExposure = mapBackSelection2(exposure, lengthArray);
else
    disp("Must input oneLimbPolys or lengthArray!");
end
outCell = zeros(size(backExposure, 1), 4); %(i,1) = rr, (i,2) = rrCI lower bound, (i,4) = n
for i = 1:size(backCondition, 1)
    rr = (sum(backCondition{i} & backExposure{i})/sum(backExposure{i})) /...
        (sum(backCondition{i} & ~backExposure{i})/sum(~backExposure{i}));
    logCI = Z*(sum(~backCondition{i} & backExposure{i})/...
        (sum(backCondition{i} & backExposure{i}) * sum(backExposure{i})) + ...
        sum(~backCondition{i} & ~backExposure{i})/...
        (sum(backCondition{i} & ~backExposure{i}) * sum(~backExposure{i})))^0.5;
    rr_CI = [exp(log(rr) - logCI), exp(log(rr) + logCI)];
    outCell(i,1) = rr;
    outCell(i,2:3) = rr_CI;
    outCell(i,4) = length(backCondition{i});
end

end

function backCell = mapBackSelection(selection, oneLimbPolys)
    backCell = cell(size(oneLimbPolys));
    cumu = 1;
    for i = 1:size(oneLimbPolys, 1)    
        backCell{i} = selection(cumu:(cumu + size(oneLimbPolys{i}, 3) - 1));
        cumu = cumu + size(oneLimbPolys{i}, 3);
    end
end

function backCell = mapBackSelection2(selection, lengthArray)
    backCell = cell(size(lengthArray));
    cumu = 1;
    for i = 1:size(lengthArray, 1)    
        backCell{i} = selection(cumu:(cumu + lengthArray(i) - 1));
        cumu = cumu + lengthArray(i);
    end
end


%% Debug examples.
% n = 100;
% cond = 34;
% condExp = 23;
% exp = 50;
% condition = zeros(n, 1);
% exposure = zeros(n,1);
% condition(1:cond) = 1;
% exposure(1:condExp) = 1;
% exposure(cond+1:cond + exp-condExp) = 1;
% 
% outCell = RelativeRisk23(condition, exposure, {zeros(3,3,n)})
