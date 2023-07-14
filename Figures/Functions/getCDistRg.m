function [cDist, rG] = getCDistRg(hPol)
    [nB,~,nCells] = size(hPol);
    cDist = nan(nCells,nB);
    rG = nan(nCells,1);
    for c=1:nCells
        currPol = hPol(:,1:3,c);
        % Dealing with bad hybs. 
        % Note that this part is lazily hard-coded for interogating either
        % the whole region (size == 75) or the orange TAD (size == 14)
        % If other segments are desired, please change the code
        % accordingly.
        if size(hPol, 1) == 75
            currPol = removeBadHybs_pol(currPol, [26, 49]);
        elseif size(hPol, 1) == 14
            currPol = removeBadHybs_pol(currPol, [26-22+1]);
        end
        newPol = CenterPolymer(currPol);
        cDist(c,:) = vecnorm(newPol, 2, 2);
        rG(c) = RadiusOfGyration(currPol);
    end
end