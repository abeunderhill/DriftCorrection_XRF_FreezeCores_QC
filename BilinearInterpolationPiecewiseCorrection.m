%% ==========================================================
% BATCH: Space-time QC correction using 3-time QC (Pre / Full / Post)
% and piecewise linear spatial weighting between QC locations.
%
% For each element i:
%  1) For each QC location j, compute drift anchors:
%       D_j(t_pre)=1
%       D_j(t_mid)=M_j(mid)/M_j(pre)
%       D_j(t_post)=M_j(post)/M_j(pre)
%     then interpolate in time => D_j(t)
%  2) For each row at position z between QC z_j and z_{j+1}:
%       D(t,z) = (1-w)D_j(t) + w D_{j+1}(t)
%  3) Correct: x_corr = x / D(t,z)
%
% Output includes Pre/Full/Post rows (optional).
%% ==========================================================

clear; clc; close all;


inFiles = {
    'FILENAME.csv'
    
};



outSuffix = '_Corr_BilinearPiecewise_260216';     
outFolder = '';                     

% Column names
scanCol = 'Scan';                   % Pre / Full / Post
timeCol = 't_hr';                   % hours
posCol  = 'position_mm_';          % mm 


elements = {'ELEMENTS TO BE CORRECTED};

% QC detection + selection tolerances (mm)
clusterTol_mm = 1.5;   
posTol_mm     = 1.0;   

% QC value summary at each QC 
qcSummary = 'mean';    

% Time interpolation for each QC location's drift function D_j(t)
timeInterp = 'linear'; 

% Apply correction to which rows?
correctWhich = 'AllPhases';  % 'FullOnly' or 'AllPhases'

% Safety
enforcePositive = true;   
minD = 0.01;               
maxD = 20;

%% ---------------- MAKE OUTPUT FOLDER  ----------------
if ~isempty(outFolder) && ~exist(outFolder,'dir')
    mkdir(outFolder);
end

%% ---------------- RUN BATCH ----------------
for f = 1:numel(inFiles)
    inFile = inFiles{f};

    fprintf('\n============================================================\n');
    fprintf('Processing %d/%d: %s\n', f, numel(inFiles), inFile);

    try
        Tout = processOneFile_SpaceTime(inFile, scanCol, timeCol, posCol, elements, ...
            clusterTol_mm, posTol_mm, qcSummary, timeInterp, correctWhich, ...
            enforcePositive, minD, maxD);

        % Output path
        [p, name, ext] = fileparts(inFile);
        if isempty(outFolder)
            outPath = fullfile(p, [name outSuffix ext]);
        else
            outPath = fullfile(outFolder, [name outSuffix ext]);
        end

        writetable(Tout, outPath);
        fprintf('Wrote: %s\n', outPath);

    catch ME
        fprintf(2, 'FAILED: %s\nReason: %s\n', inFile, ME.message);
    end
end

fprintf('\nDone.\n');

% FUNCTIONS

function Tout = processOneFile_SpaceTime(inFile, scanCol, timeCol, posCol, elements, ...
    clusterTol_mm, posTol_mm, qcSummary, timeInterp, correctWhich, ...
    enforcePositive, minD, maxD)

    T = readtable(inFile);

    % Required columns
    req = {scanCol, timeCol, posCol};
    for i = 1:numel(req)
        if ~ismember(req{i}, T.Properties.VariableNames)
            error("Missing required column '%s'.", req{i});
        end
    end

    scanStr = string(T.(scanCol));
    isPre  = scanStr=="Pre";
    isFull = scanStr=="Full";
    isPost = scanStr=="Post";

    if ~any(isPre) || ~any(isFull) || ~any(isPost)
        error("Scan column must include Pre, Full, Post.");
    end

    % Detect QC positions from Pre positions (repeated)
    prePos = T.(posCol)(isPre);
    prePos = prePos(isfinite(prePos));
    qcPos = detectQCPositions(prePos, clusterTol_mm);

    if numel(qcPos) < 3
        error("Detected only %d QC positions; expected 3 or 4.", numel(qcPos));
    end

    qcPos = sort(qcPos(:));
    nQC = numel(qcPos);

    fprintf('QC positions (mm): ');
    fprintf('%.2f ', qcPos);
    fprintf('\n');

    % Decide which rows to correct (All Phases includes Pre/Mid/Post Corr)
    switch lower(correctWhich)
        case 'fullonly'
            doCorr = isFull;
        case 'allphases'
            doCorr = true(height(T),1);
        otherwise
            error("correctWhich must be 'FullOnly' or 'AllPhases'.");
    end

    Tout = T;  % keep all rows to validate 

    tAll = Tout.(timeCol);
    zAll = Tout.(posCol);

    %  Loop elements 
    for e = 1:numel(elements)
        el = elements{e};

        if ~ismember(el, T.Properties.VariableNames)
            fprintf("  %s: missing column, skipping.\n", el);
            continue;
        end

        % Build per-QC drift functions D_j(t)
        % Each QC j provides anchors at (t_pre,1), (t_mid,Dmid), (t_post,Dpost)
        Dj = cell(nQC,1);  % function handles (via interp1) not stored; store anchors instead
        tAnch = nan(nQC,3);
        DAnch = nan(nQC,3);

        goodQC = true(nQC,1);

        for j = 1:nQC
            z0 = qcPos(j);

            mPre  = isPre  & abs(T.(posCol)-z0) <= posTol_mm;
            mMid  = isFull & abs(T.(posCol)-z0) <= posTol_mm; % "during" within Full at this location
            mPost = isPost & abs(T.(posCol)-z0) <= posTol_mm;

            if ~any(mPre) || ~any(mMid) || ~any(mPost)
                goodQC(j) = false;
                continue;
            end

            % Times
            tpre  = median(T.(timeCol)(mPre),  'omitnan');
            tmid  = median(T.(timeCol)(mMid),  'omitnan');
            tpost = median(T.(timeCol)(mPost), 'omitnan');

            % Values
            vpre  = T.(el)(mPre);
            vmid  = T.(el)(mMid);
            vpost = T.(el)(mPost);

            switch lower(qcSummary)
                case 'mean'
                    Mpre  = mean(vpre,  'omitnan');
                    Mmid  = mean(vmid,  'omitnan');
                    Mpost = mean(vpost, 'omitnan');
                case 'median'
                    Mpre  = median(vpre,  'omitnan');
                    Mmid  = median(vmid,  'omitnan');
                    Mpost = median(vpost, 'omitnan');
                otherwise
                    error("qcSummary must be 'mean' or 'median'.");
            end

            if ~isfinite(Mpre) || Mpre==0 || ~isfinite(Mmid) || ~isfinite(Mpost)
                goodQC(j) = false;
                continue;
            end

            % Drift anchors relative to PRE at THIS QC location
            Dpre  = 1;
            Dmid  = Mmid / Mpre;
            Dpost = Mpost / Mpre;

            tAnch(j,:) = [tpre tmid tpost];
            DAnch(j,:) = [Dpre Dmid Dpost];
        end

        % Require at least 2 QC locations to interpolate in space
        if sum(goodQC) < 2
            fprintf("  %s: not enough valid QC locations, skipping.\n", el);
            continue;
        end

        % For any invalid QC, fill by copying nearest valid QC
        [tAnch, DAnch, goodQC] = fillMissingQC(tAnch, DAnch, goodQC, qcPos);

        % Precompute D_j(t) for all rows we will correct 
        % Djt(j,:) = D at each row time
        Djt = nan(nQC, numel(tAll));

        for j = 1:nQC
            tj = tAnch(j,:).';
            Djv = DAnch(j,:).';

            % Sort by time 
            [tj, ord] = sort(tj);
            Djv = Djv(ord);

            % Interpolate drift vs time for this QC location
            Djt(j,:) = interp1(tj, Djv, tAll, timeInterp, 'extrap');
        end

        % Clamp / safety
        if enforcePositive
            bad = ~isfinite(Djt) | (Djt <= 0);
            Djt(bad) = 1;
        end
        Djt = min(max(Djt, minD), maxD);

        % Compute space-weighting D(t,z) for each row
        Dtz = computeSpaceWeightedD(zAll, qcPos, Djt);

        % Apply only to desired rows
        raw = Tout.(el);
        corr = raw;
        corr(doCorr) = raw(doCorr) ./ Dtz(doCorr);

        % Store
        Tout.([el '_Dtz'])   = Dtz;
        Tout.([el '_corrST'])= corr;

        fprintf("  %s: corrected (%s)\n", el, correctWhich);
    end
end

function Dtz = computeSpaceWeightedD(z, qcPos, Djt)
% Piecewise linear in space between QC positions.
% Inputs:
%   z     : Nx1 positions for rows
%   qcPos : nQCx1 QC positions 
%   Djt   : nQCxN drift at each QC location for each row time
% Output:
%   Dtz   : Nx1 space-weighted drift

    z = z(:);
    nQC = numel(qcPos);
    N = numel(z);

    Dtz = nan(N,1);

    % For each row, find segment [qcPos(k), qcPos(k+1)]
    for ii = 1:N
        zi = z(ii);

        if zi <= qcPos(1)
            % above first QC: use first
            Dtz(ii) = Djt(1,ii);
            continue;
        end
        if zi >= qcPos(end)
            % below last QC: use last
            Dtz(ii) = Djt(end,ii);
            continue;
        end

        % find k such that qcPos(k) <= zi < qcPos(k+1)
        k = find(qcPos <= zi, 1, 'last');
        if k >= nQC
            k = nQC-1;
        end

        z1 = qcPos(k);
        z2 = qcPos(k+1);
        w = (zi - z1) / (z2 - z1);

        D1 = Djt(k,ii);
        D2 = Djt(k+1,ii);

        Dtz(ii) = (1-w)*D1 + w*D2;
    end
end

function qcPos = detectQCPositions(prePos, tol_mm)
% Detect distinct QC locations from Pre positions using clustering.
    p = sort(prePos(:));
    p = p(isfinite(p));
    if isempty(p), qcPos = []; return; end

    clusters = {};
    startIdx = 1;
    for i = 2:numel(p)
        if (p(i) - p(i-1)) > tol_mm
            clusters{end+1} = p(startIdx:i-1); 
            startIdx = i;
        end
    end
    clusters{end+1} = p(startIdx:end);

    qcPos = zeros(numel(clusters),1);
    for k = 1:numel(clusters)
        qcPos(k) = median(clusters{k});
    end
    qcPos = sort(qcPos);
end

function [tAnch, DAnch, goodQC] = fillMissingQC(tAnch, DAnch, goodQC, qcPos)
% If a QC anchor location is missing copy anchors from nearest valid QC.
    nQC = numel(qcPos);
    validIdx = find(goodQC);

    for j = 1:nQC
        if goodQC(j), continue; end

        % nearest valid
        [~, k] = min(abs(qcPos(validIdx) - qcPos(j)));
        jNear = validIdx(k);

        tAnch(j,:) = tAnch(jNear,:);
        DAnch(j,:) = DAnch(jNear,:);
        goodQC(j)  = true;
    end
end
