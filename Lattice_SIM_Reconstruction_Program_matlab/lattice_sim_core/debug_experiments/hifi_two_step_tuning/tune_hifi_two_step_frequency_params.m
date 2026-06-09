function tune_hifi_two_step_frequency_params()
%TUNE_HIFI_TWO_STEP_FREQUENCY_PARAMS Search two-step Wiener frequency settings.

scriptDir = fileparts(mfilename('fullpath'));
coreDir = fileparts(fileparts(scriptDir));
addpath(fullfile(coreDir, 'functions'));

resultPath = fullfile(coreDir, 'output', 'diagnostics', 'result.mat');
loaded = load(resultPath, 'result');
combine = loaded.result.diagnostics.combine;

outputDir = fullfile(coreDir, 'output', 'diagnostics', 'hifi_two_step_tuning');
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

grid.wienerW1 = [0.15, 0.20, 0.30, 0.40, 0.60, 0.90, 1.20];
grid.wienerW2 = [0.015, 0.020, 0.030, 0.040, 0.060, 0.100];
grid.apodizationRadius = [0.75, 1.00, 1.25, 1.50, 1.75];
grid.hifiDenominatorScaleW1 = [0.60, 0.80, 1.00, 1.20];
grid.hifiDenominatorScaleW2 = [0.80, 1.00, 1.20];

metrics = evaluateGrid(combine, grid);
metrics = sortrows(metrics, 'score', 'descend');

writetable(metrics, fullfile(outputDir, 'hifi_two_step_frequency_tuning.csv'));
save(fullfile(outputDir, 'hifi_two_step_frequency_tuning.mat'), ...
    'metrics', 'grid', '-v7.3');
writeSummary(outputDir, metrics);

disp(metrics(1:min(10, height(metrics)), :));
fprintf('Best tuning summary: %s\n', fullfile(outputDir, 'hifi_two_step_frequency_tuning_summary.txt'));
end

function metrics = evaluateGrid(combine, grid)
fftDirectlyCombined = combine.fftDirectlyCombined;
wienerDenominatorW1Base = getDiagnosticField(combine, ...
    'wienerDenominatorW1Base', combine.blendDenominator);
wienerDenominatorW2Base = getDiagnosticField(combine, ...
    'wienerDenominatorW2Base', combine.blendDenominator);
hifiMask = combine.hifiMask > 0;
hifiOtf = makeExtendedHifiOtf(hifiMask);
physicalMask = abs(combine.physicalOtfValues) > 0 & hifiMask;
expandedMask = hifiMask & ~physicalMask;
sidebandMasks = makeSidebandMasks(combine, physicalMask, hifiMask);

inputSupportRms = rmsMagnitude(fftDirectlyCombined, hifiMask);
inputTotalEnergy = energyInMask(fftDirectlyCombined, hifiMask);

rows = numel(grid.wienerW1) * numel(grid.wienerW2) * ...
    numel(grid.apodizationRadius) * numel(grid.hifiDenominatorScaleW1) * ...
    numel(grid.hifiDenominatorScaleW2);
result = struct( ...
    'wienerW1', cell(rows, 1), ...
    'wienerW2', cell(rows, 1), ...
    'apodizationRadius', cell(rows, 1), ...
    'hifiDenominatorScaleW1', cell(rows, 1), ...
    'hifiDenominatorScaleW2', cell(rows, 1), ...
    'score', cell(rows, 1), ...
    'expandedToCenterRms', cell(rows, 1), ...
    'sidebandMeanToCenterRms', cell(rows, 1), ...
    'sidebandMinToCenterRms', cell(rows, 1), ...
    'sidebandBalance', cell(rows, 1), ...
    'totalEnergyRatioToInput', cell(rows, 1), ...
    'p99GainRatioToInput', cell(rows, 1));

row = 0;
for w1 = grid.wienerW1
    for w2 = grid.wienerW2
        for apoRadius = grid.apodizationRadius
            apodizationMask = makeGaussianApodization(hifiMask, apoRadius);
            for scaleW1 = grid.hifiDenominatorScaleW1
                wienerW1 = hifiOtf ./ (wienerDenominatorW1Base * scaleW1 + w1 ^ 2);
                spectrumAfterW1 = fftDirectlyCombined .* wienerW1 .* hifiMask;
                for scaleW2 = grid.hifiDenominatorScaleW2
                    wienerW2 = apodizationMask ./ (wienerDenominatorW2Base * scaleW2 + w2 ^ 2);
                    spectrum = spectrumAfterW1 .* wienerW2 .* hifiMask;
                    row = row + 1;
                    result(row) = scoreSpectrum(spectrum, fftDirectlyCombined, ...
                        physicalMask, expandedMask, sidebandMasks, inputSupportRms, ...
                        inputTotalEnergy, w1, w2, apoRadius, scaleW1, scaleW2);
                end
            end
        end
    end
end

metrics = struct2table(result);
end

function value = getDiagnosticField(diagnostics, fieldName, fallbackValue)
if isfield(diagnostics, fieldName)
    value = diagnostics.(fieldName);
else
    value = fallbackValue;
end
end

function row = scoreSpectrum(spectrum, inputSpectrum, physicalMask, expandedMask, ...
    sidebandMasks, inputSupportRms, inputTotalEnergy, w1, w2, apoRadius, scaleW1, scaleW2)
centerRms = rmsMagnitude(spectrum, physicalMask);
expandedRms = rmsMagnitude(spectrum, expandedMask);
sidebandRms = zeros(1, numel(sidebandMasks));
for idx = 1:numel(sidebandMasks)
    sidebandRms(idx) = rmsMagnitude(spectrum, sidebandMasks{idx});
end

expandedToCenter = expandedRms / max(eps, centerRms);
sidebandMeanToCenter = mean(sidebandRms) / max(eps, centerRms);
sidebandMinToCenter = min(sidebandRms) / max(eps, centerRms);
sidebandBalance = min(sidebandRms) / max(eps, mean(sidebandRms));
totalEnergyRatio = energyInMask(spectrum, true(size(spectrum))) / max(eps, inputTotalEnergy);
p99GainRatio = percentileMagnitude(spectrum, 99) / max(eps, inputSupportRms);

targetExpanded = 0.18;
visibility = exp(-0.5 * ((log(max(expandedToCenter, eps)) - log(targetExpanded)) / log(1.8)) ^ 2);
sidebandVisibility = min(sidebandMeanToCenter / targetExpanded, 1.5);
energyPenalty = 1 / (1 + max(0, totalEnergyRatio - 0.45) ^ 2);
gainPenalty = 1 / (1 + max(0, p99GainRatio - 2.5) ^ 2);
score = visibility * (0.55 + 0.45 * sidebandVisibility) * ...
    (0.35 + 0.65 * sidebandBalance) * energyPenalty * gainPenalty;

row.wienerW1 = w1;
row.wienerW2 = w2;
row.apodizationRadius = apoRadius;
row.hifiDenominatorScaleW1 = scaleW1;
row.hifiDenominatorScaleW2 = scaleW2;
row.score = score;
row.expandedToCenterRms = expandedToCenter;
row.sidebandMeanToCenterRms = sidebandMeanToCenter;
row.sidebandMinToCenterRms = sidebandMinToCenter;
row.sidebandBalance = sidebandBalance;
row.totalEnergyRatioToInput = totalEnergyRatio;
row.p99GainRatioToInput = p99GainRatio;
end

function masks = makeSidebandMasks(combine, physicalMask, hifiMask)
masks = cell(1, 4);
for idx = 1:4
    mask = combine.bandTaperMasks{idx + 1} > 0 & ~physicalMask & hifiMask;
    if nnz(mask) == 0
        mask = combine.bandTaperMasks{idx + 1} > 0 & hifiMask;
    end
    masks{idx} = mask;
end
end

function hifiOtf = makeExtendedHifiOtf(supportMask)
[h, w] = size(supportMask);
[x, y] = meshgrid((1:w) - floor(w/2) - 1, (1:h) - floor(h/2) - 1);
radius = hypot(x, y);
supportRadius = max(radius(supportMask), [], 'all');
cutoffRadius = max(1, supportRadius + 1);
rho = radius ./ cutoffRadius;
hifiOtf = zeros(h, w);
inside = rho <= 1;
hifiOtf(inside) = (2 / pi) * (acos(rho(inside)) - ...
    rho(inside) .* sqrt(1 - rho(inside) .^ 2));
hifiOtf = hifiOtf .* double(supportMask);
end

function apo = makeGaussianApodization(supportMask, apoRadius)
[h, w] = size(supportMask);
[x, y] = meshgrid(1:w, 1:h);
center = [floor(h/2) + 1, floor(w/2) + 1];
radius = hypot((x - center(2)) * 2 / w, (y - center(1)) * 2 / h);
apo = exp(-0.5 * (radius ./ apoRadius .* sqrt(2*log(2))).^2);
apo = apo .* double(supportMask);
end

function value = rmsMagnitude(spectrum, mask)
values = abs(spectrum(mask));
if isempty(values)
    value = 0;
else
    value = sqrt(mean(values .^ 2));
end
end

function value = energyInMask(spectrum, mask)
values = abs(spectrum(mask));
value = sum(values .^ 2);
end

function value = percentileMagnitude(spectrum, percentile)
values = abs(spectrum(:));
values = sort(values);
idx = max(1, min(numel(values), round(numel(values) * percentile / 100)));
value = values(idx);
end

function writeSummary(outputDir, metrics)
best = metrics(1, :);
fid = fopen(fullfile(outputDir, 'hifi_two_step_frequency_tuning_summary.txt'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Best HiFi-style two-step frequency tuning\n');
fprintf(fid, 'wienerW1\t%.6g\n', best.wienerW1);
fprintf(fid, 'wienerW2\t%.6g\n', best.wienerW2);
fprintf(fid, 'apodizationRadius\t%.6g\n', best.apodizationRadius);
fprintf(fid, 'hifiDenominatorScaleW1\t%.6g\n', best.hifiDenominatorScaleW1);
fprintf(fid, 'hifiDenominatorScaleW2\t%.6g\n', best.hifiDenominatorScaleW2);
fprintf(fid, 'score\t%.6g\n', best.score);
fprintf(fid, 'expandedToCenterRms\t%.6g\n', best.expandedToCenterRms);
fprintf(fid, 'sidebandMeanToCenterRms\t%.6g\n', best.sidebandMeanToCenterRms);
fprintf(fid, 'sidebandMinToCenterRms\t%.6g\n', best.sidebandMinToCenterRms);
fprintf(fid, 'sidebandBalance\t%.6g\n', best.sidebandBalance);
fprintf(fid, 'totalEnergyRatioToInput\t%.6g\n', best.totalEnergyRatioToInput);
fprintf(fid, 'p99GainRatioToInput\t%.6g\n', best.p99GainRatioToInput);
end
