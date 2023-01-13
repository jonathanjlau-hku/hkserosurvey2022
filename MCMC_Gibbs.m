function chain = MCMC_Gibbs(logLikelihood, logPrior, chain)

if (~isempty(chain.stepSTD))
    chain.stepSTD = chain.stepSTD;
else
    chain.stepSTD = repmat(0.1, 1, chain.numParameters);
end
if (size(chain.stepSTD,1)>1)
    chain.stepSTD = chain.stepSTD';
end

if (~isfolder(chain.dir))
    mkdir(chain.dir);
end
[rsize, csize] = size(chain.startingPoint);
if (rsize~=1 && csize==1)
    chain.startingPoint = chain.startingPoint';
end

backupTimeInterval = 60;   % save work every 60 min

rvStorageSize = min(chain.totalNumSteps,1e6);
chain.probAccept = mean(chain.probAcceptRange)*ones(1, chain.numParameters);

chain.record = zeros(chain.totalNumSteps, chain.numParameters);
chain.logLikelihood = zeros(chain.totalNumSteps, 1);
chain.logPosterior = zeros(chain.totalNumSteps, 1);
chainGood = false;

chain.bestParameters = chain.startingPoint;
chain.bestLogLikelihood = -1e300;
chain.bestLogPosterior = -1e300;

chainNum = 1;
pTarget = mean(chain.probAcceptRange);
pTarget_logit = log(pTarget/(1-pTarget));


while (chainGood==false)
    x = chain.bestParameters;
    chain.startingPoint = x;
    xLogLikelihood = logLikelihood(x, chain.data);
    xLogPosterior = xLogLikelihood+logPrior(x, chain.data);
    
    chain.record(1,:) = x;
    chain.logLikelihood(1) = xLogLikelihood;
    chain.logPosterior(1) = xLogPosterior;
    numAccept = ones(1, chain.numParameters);
    numSteps = ones(1, chain.numParameters);
    
    if (chainNum>2 && chainNum<20)
        stepSTD_target = exp((pTarget_logit-regressionCoef(:,1))./regressionCoef(:,2));
        scalingFactor = min(10, max(0.1, stepSTD_target'./chain.stepSTD));
        fI = find(chain.probAccept>0.9);
        if (~isempty(fI))
            scalingFactor(fI) = 2;
        end
        fI = find(chain.probAccept<0.1);
        if (~isempty(fI))
            scalingFactor(fI) = 0.5;
        end
        chain.stepSTD = chain.stepSTD.*scalingFactor;
    else
        for para=1:chain.numParameters
            if (chain.probAccept(para)<chain.probAcceptRange(1) )
                chain.stepSTD(para) = 0.5*chain.stepSTD(para);
            elseif (chain.probAccept(para)>chain.probAcceptRange(2))
                chain.stepSTD(para) = 2*chain.stepSTD(para);
            else
                chain.stepSTD(para) = (1+unifrnd(-0.01,0.01))*chain.stepSTD(para);
            end
        end
    end
    fI = find(chain.UB==1);
    chain.stepSTD(fI) = min(chain.stepSTD(fI), 0.5);
%     chain.stepSTD(fI) = chain.stepSTD(fI);
    fI = chain.probAccept>chain.probAcceptRange(1) & chain.probAccept<chain.probAcceptRange(2);
    for para=1:chain.numParameters
        chain.num_bounces{para} = [];
    end
    
    numStepsPerUpdate = ones(1,chain.numParameters);
    if (sum(fI)>0)
        numStepsPerUpdate(fI==false) = 2+ceil(sum(fI==false)/sum(fI));
    end
    paraUpdate = 1;
    tic;
    totalSecElapsed = 0;
    
    probStorage = unifrnd(0, 1, rvStorageSize, 1);
    probStorageIndex = 1;
    
    stepStorage = normrnd(0, 1, rvStorageSize*chain.numParameters, 1);
    stepStorageIndex = 1;
    iter = 2;
    tic_backup = tic;
    tic_start = tic;
    displayCounter = chain.numStepsBetweenDisplay;
    

    while (iter<=chain.totalNumSteps)
        for subStep=1:numStepsPerUpdate(paraUpdate)
          
            y = x;
            
            if (chain.logSpace(paraUpdate)==false)
                y(paraUpdate) = x(paraUpdate)+stepStorage(stepStorageIndex)*chain.stepSTD(paraUpdate);
            else
                y(paraUpdate) = x(paraUpdate)*exp(stepStorage(stepStorageIndex)*chain.stepSTD(paraUpdate));
            end
            stepStorageIndex = stepStorageIndex+1;
            if (stepStorageIndex+chain.numParameters>rvStorageSize)
                stepStorage = normrnd(0, 1, rvStorageSize*chain.numParameters, 1);
                stepStorageIndex = 1;
            end
            
            jj = 0;
            while (y(paraUpdate)<chain.LB(paraUpdate) || y(paraUpdate)>chain.UB(paraUpdate))
                if (y(paraUpdate)<chain.LB(paraUpdate))
                    y(paraUpdate) = chain.LB(paraUpdate)+(chain.LB(paraUpdate)-y(paraUpdate));
                else
                    y(paraUpdate) = chain.UB(paraUpdate)-(y(paraUpdate)-chain.UB(paraUpdate));
                end
                jj = jj+1;
            end
            if (jj>0)
                chain.num_bounces{paraUpdate}(end+1) = jj;
            end
            if (jj>1)
                 [paraUpdate jj]
            end
            yLogLikelihood = logLikelihood(y, chain.data);
            yLogPosterior = yLogLikelihood+logPrior(y, chain.data);
            
            if (yLogLikelihood>-1e300 && probStorage(probStorageIndex)<min(1,exp(yLogPosterior-xLogPosterior)))
                x = y;
                xLogLikelihood = yLogLikelihood;
                xLogPosterior = yLogPosterior;
                
                
                numAccept(paraUpdate) = numAccept(paraUpdate)+1;
            end
            
            probStorageIndex = probStorageIndex+1;
            if (probStorageIndex==rvStorageSize)
                probStorage = unifrnd(0, 1, rvStorageSize, 1);
                probStorageIndex = 1;
            end
            
            if (xLogLikelihood>chain.bestLogLikelihood)
                chain.bestLogLikelihood = xLogLikelihood;
            end
            if (xLogPosterior>chain.bestLogPosterior)
                chain.bestLogPosterior =  xLogPosterior;
                chain.bestParameters = x;
            end
            
            numSteps(paraUpdate) = numSteps(paraUpdate)+1;
            
            chain.iter = iter;
            chain.record(iter,:) = x;
            chain.logLikelihood(iter) = xLogLikelihood;
            chain.logPosterior(iter) = xLogPosterior; 
            chain.probAccept(paraUpdate) = numAccept(paraUpdate)/numSteps(paraUpdate);
           
            iter = iter+1;
        end
        if (iter>chain.minNumStepsBeforeRestart && (chain.probAccept(paraUpdate)<chain.probAcceptRange(1) || ...
                chain.probAccept(paraUpdate)>chain.probAcceptRange(2) || var(numStepsPerUpdate)>0))
            chainGood = false;
            
            fI = find(chain.probAccept < chain.probAcceptRange(1) | chain.probAccept > chain.probAcceptRange(2));
            AA = [fI; chain.probAccept(fI); chain.stepSTD(fI)];
            
            if (chain.probAccept(paraUpdate)<chain.probAcceptRange(1) || chain.probAccept(paraUpdate)>chain.probAcceptRange(2))
                fprintf(['\nRestarting MCMC after ' num2str(iter)  ' steps because the following P(Accept) is out of range. ']);
                fprintf('\nParameters\t\tP(Accept)\t\tchain.stepSTD\n');
                fprintf('%5d\t\t\t%5.2f\t\t\t%8.3e\n', AA);
            else
                fprintf(['\nRestarting MCMC after ' num2str(iter)  ' steps after tuning the proposal distributions.\n']);
            end
            
            stepSTD_record(chainNum,:) = chain.stepSTD;
            probAccept_record(chainNum, :) = chain.probAccept;
            if (chainNum>1)
                figure(500);
                clf;
                totalNumPlots = chain.numParameters;
                numFigCol = min(5, totalNumPlots);
                numFigRow = ceil(totalNumPlots/numFigCol);
                for para=1:chain.numParameters
                    subplot(numFigRow, numFigCol, para);
                    hold on;
                    X_regress = [ones(chainNum, 1) log(stepSTD_record(:,para))];
                    Y_regress = log(probAccept_record(:,para)./(1-max(1e-10,probAccept_record(:,para))));
                    [B, ~] = regress(Y_regress, X_regress);
                    plot(log(stepSTD_record(:,para)), Y_regress, 'o');
                    plot(log(stepSTD_record(:,para)), X_regress*B);
                    title(['Slope = ' num2str(B(2))]);
                    regressionCoef(para,:) = B;
                    xlabel('Log of std dev of proposal');
                    ylabel('Logit of P(Accept)');
                    title(['Para ' num2str(para)]);
                end
                plotMCMC(chain, false);
            end
        else
            chainGood = true;
            
        end
        if (chainGood==false)
            chainNum = chainNum+1;
            break;
        end
        
        if (paraUpdate==chain.numParameters)
            paraUpdate = 1;
        else
            paraUpdate = paraUpdate+1;
        end
        
%         totalSecElapsed = totalSecElapsed+toc;
%         tic;
        if (iter>=displayCounter)
            displayCounter = displayCounter+chain.numStepsBetweenDisplay;
            if (~isempty(chain.dir))
                fprintf(['\nExperiment: ' chain.dir '\nAlgorithm: Gibbs']);
            end
            
            fprintf(['\nMCMC progress: ' num2str(iter/chain.totalNumSteps*100, '%.1f') '%% complete.']);
            fprintf(['\nBest log-likelihood: ' num2str(chain.bestLogLikelihood)]);
            fprintf(['\nBest log-posterior: ' num2str(chain.bestLogPosterior)]);
            totalSecElapsed = toc(tic_start)/60;
            fprintf(['\nTime elapsed = ' num2str(totalSecElapsed, '%.1f') ' min. Remaining time = ' ...
                num2str(totalSecElapsed/iter*(chain.totalNumSteps-iter), '%.1f') ' min.\n']);
            plotMCMC(chain, true);
            
            if ((toc-tic_backup)/60 > backupTimeInterval)
                tic_backup = tic;
                fname = [chain.dir '/MCMC_Records'];
                save(fname, 'chain');
            end
        end
    end
end
plotMCMC(chain, true);
end

function funval = plotMCMC(chain, runAutocorrelation)

nn = chain.iter-1;
chain.numParameters = size(chain.record,2);
if (runAutocorrelation)
    lagmax = max(2,round(nn/10));
    for para=1:chain.numParameters
        xx = chain.record(1:nn,para);
        [acf(:,para),lags,bounds] = autocorr(xx, 'NumLags', lagmax);
    end
end

h = figure(100);
set(h, 'Name', 'Traceplots');
clf;
numFigCol = ceil(sqrt(chain.numParameters));
numFigRow = ceil(chain.numParameters/numFigCol);
subplotIndex = 0;
for ii=1:chain.numParameters
    subplotIndex = subplotIndex+1;
    subplot(numFigRow, numFigCol, subplotIndex);
    plot(chain.record(1:nn,ii));
    xlim([0 nn]);
    ylabel(['Para ' num2str(ii)]);
end

h = figure(101);
set(h, 'Name', 'Posterior distributions');
clf;
numFigCol = ceil(sqrt(chain.numParameters));
numFigRow = ceil(chain.numParameters/numFigCol);
nn = chain.iter-1;
subplotIndex = 0;
for ii=1:chain.numParameters
    subplotIndex = subplotIndex+1;
    subplot(numFigRow, numFigCol, subplotIndex);
    hist(chain.record(1:nn,ii), round(sqrt(nn))/2);
    xlabel(['Para ' num2str(ii)]);
end


h = figure(103);
set(h, 'Name', 'Other stats');
numFigRow = 1;
if (runAutocorrelation)
    numFigCol = 3;
else
    numFigCol = 2;
end
subplotIndex = 0;
subplotIndex = subplotIndex+1;
subplot(numFigRow, numFigCol, subplotIndex)
plot(chain.logLikelihood(1:nn));
xlim([0 nn]);
ylabel('logL');
xlabel('Iteration');

subplotIndex = subplotIndex+1;
subplot(numFigRow, numFigCol, subplotIndex)
bar(chain.probAccept);
xlim([0 chain.numParameters+1])
ylim([0 1]);
ylabel('P(Accept)');
xlabel('Parameter');

if (runAutocorrelation)
    subplotIndex = subplotIndex+1;
    subplot(numFigRow, numFigCol, subplotIndex)
    plot(lags,acf);
    xlim([1 max(lags)]);
    lt = {};
    for para=1:chain.numParameters
        lt{end+1} = ['Para ' num2str(para)];
    end
    legend(lt);
    xlabel('Lag');
    ylabel('Autocorrelation');
end

end