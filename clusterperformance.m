function [timeperformance,ARI,orderofdistances,outputinfo,LBKeogh] = clusterperformance(N,repetitions,LBKeogh)

    clusters = 2;

orderofdistances = {'Euclid','LB_Keogh','Cepstral', 'Extended','H2','HInf'};
outputinfo = {'Length of series: ', N; 'LBKeogh',LBKeogh};
timeperformance = zeros(repetitions,6);
ARI = zeros(repetitions,6);

for r = 1:repetitions
    
  disp(['Size 2^' num2str(log2(N)) ' Repetition: ' num2str(r) '/' num2str(repetitions)])
  
  [inputs,outputs,clusterid] = generateelectriccircuits(N);
  
  
    normoutputs = zscore(outputs')';
    norminputs = zscore(inputs')';
    [n,k] = size(outputs);
  
    disp('Start Euclidean Distance')
    
    tic
    DistOutput = squareform(pdist(normoutputs));
    DistInput = squareform(pdist(norminputs));
% Change alpha to a number between 0 and 1 to explicitly add a cost term involving the input to the Euclidean distance.
    alpha = 0;
    VectorEuclid = squareform((1-alpha)*DistOutput + alpha*DistInput);
    
    Zeuclid = linkage(VectorEuclid);
    % ceuclid = cluster(Zeuclid,'cutoff',10^3,'criterion','distance');
    ceuclid = cluster(Zeuclid,'maxclust',clusters);
    timeperformance(r,1) = toc;
    ARI(r,1) = adjrandindex(ceuclid,clusterid);
 
    if LBKeogh
        disp('Start LB_Keogh')
        tic
        DistLBKeogh = zeros(n,n);
        w=floor(k*(5/100));
            for i = 1:n
                for j = 1:i-1
                     DistLBKeogh(i,j) = dLBKeogh(normoutputs(i,:)', normoutputs(j,:)', 1);
                     DistLBKeogh(j,i) = DistLBKeogh(i,j);
                end
                disp(['Size 2^' num2str(log2(N)) ' Repetition: ' num2str(r) '/' num2str(repetitions)  ' LBKeogh: ' num2str(i) '/' num2str(n)])
            end
        VectorLBKeogh = squareform(DistLBKeogh);
        
% Alternatively, change the previous with the following to add the inputs in explicitly.
% This was not done by default as it increases computation time.
%         alpha = 0.5;
%         DistLBKeoghOutputs = zeros(n,n);
%         DistLBKeoghInputs = zeros(n,n);
%         w=floor(k*(5/100));
%             for i = 1:n
%                 for j = 1:i-1
%                      DistLBKeoghOutputs(i,j) = dLBKeogh(normoutputs(i,:)', normoutputs(j,:)', 1);
%                      DistLBKeoghOutputs(j,i) = DistLBKeoghOutputs(i,j);
%                      DistLBKeoghInputs(i,j) = dLBKeogh(norminputs(i,:)', norminputs(j,:)', 1);
%                      DistLBKeoghInputs(j,i) = DistLBKeoghInputs(i,j);
%                 end
%                 disp(['Size 2^' num2str(log2(N)) ' Repetition: ' num2str(r) '/' num2str(repetitions)  ' LBKeogh: ' num2str(i) '/' num2str(n)])
%             end
%         VectorLBKeogh = squareform((1-alpha)*DistLBKeoghOutputs + alpha*DistLBKeoghInputs);

        ZLBKeogh = linkage(VectorLBKeogh);
        % cLBKeogh = cluster(ZLBKeogh,'cutoff',10^3,'criterion','distance');
        cLBKeogh = cluster(ZLBKeogh,'maxclust',clusters);
        timeperformance(r,2) = toc;
        ARI(r,2) = adjrandindex(cLBKeogh,clusterid);
    end
    
    disp('Start Original Cepstral Distance')
    tic
    DistCepstral = zeros(n,n);

     if (k > 2^7)
        OutputCeps = zeros(n,size(ifft(log(pwelch(outputs(1,:),[],[],'twosided')),'symmetric'),1));
        cepstrumcutoff = floor(size(ifft(log(pwelch(outputs(1,:),[],[],'twosided')),'symmetric'),1));
        weights = 1:1:cepstrumcutoff-1;
        for i = 1:n
            OutputCeps(i,:) = ifft(log(pwelch(outputs(i,:),[],[],'twosided')),'symmetric');
        end
        disp('Welch')
     else
        OutputCeps = zeros(n,size(ifft(log(pmtm(outputs(1,:),'twosided')),'symmetric'),1));
        cepstrumcutoff = floor(size(ifft(log(pmtm(outputs(1,:),'twosided')),'symmetric'),1));
        weights = 1:1:cepstrumcutoff-1;
            for i = 1:n
            OutputCeps(i,:) = ifft(log(pmtm(outputs(i,:),'twosided')),'symmetric');
            end
        disp('mtm')
     end

    clear i

    for i = 1:n

        for j = 1:i-1

            DistCepstral(i,j) = weights*((OutputCeps(i,2:cepstrumcutoff)' - OutputCeps(j,2:cepstrumcutoff)').^2);
            DistCepstral(j,i) = DistCepstral(i,j);
        end
        disp(['Size 2^' num2str(log2(N)) ' Repetition: ' num2str(r) '/' num2str(repetitions)  ' Cepstral: ' num2str(i) '/' num2str(n)])
    end

    VectorCepstral = squareform(DistCepstral);
    Zcepstral = linkage(VectorCepstral);
    % ccepstral = cluster(Zcepstral,'cutoff',0.35,'criterion','distance');
    ccepstral = cluster(Zcepstral,'maxclust',clusters);
    timeperformance(r,3) = toc;
    ARI(r,3) = adjrandindex(ccepstral,clusterid);
    
    disp('Start Extended Cepstral Distance')
    tic
    DistExtended = zeros(n,n);
    
     if (k > 2^7)
        InputCeps = zeros(n,size(ifft(log(pwelch(inputs(1,:),[],[],'twosided')),'symmetric'),1));
        OutputCeps = zeros(n,size(ifft(log(pwelch(outputs(1,:),[],[],'twosided')),'symmetric'),1));
        cepstrumcutoff = floor(size(ifft(log(pwelch(outputs(1,:),[],[],'twosided')),'symmetric'),1));
        weights = 1:1:cepstrumcutoff-1;
        for i = 1:n
            InputCeps(i,:) = ifft(log(pwelch(inputs(i,:),[],[],'twosided')),'symmetric');
            OutputCeps(i,:) = ifft(log(pwelch(outputs(i,:),[],[],'twosided')),'symmetric');
        end
        disp('Welch')
     else
        InputCeps = zeros(n,size(ifft(log(pmtm(inputs(1,:),'twosided')),'symmetric'),1));
        OutputCeps = zeros(n,size(ifft(log(pmtm(outputs(1,:),'twosided')),'symmetric'),1));
        cepstrumcutoff = floor(size(ifft(log(pmtm(outputs(1,:),'twosided')),'symmetric'),1));
        weights = 1:1:cepstrumcutoff-1;
            for i = 1:n
            InputCeps(i,:) = ifft(log(pmtm(inputs(i,:),'twosided')),'symmetric');
            OutputCeps(i,:) = ifft(log(pmtm(outputs(i,:),'twosided')),'symmetric');
            end
        disp('mtm')
     end

    for i = 1:n

        for j = 1:i-1

            DistExtended(i,j) = weights*(((OutputCeps(i,2:cepstrumcutoff) - InputCeps(i,2:cepstrumcutoff))' - (OutputCeps(j,2:cepstrumcutoff) - InputCeps(j,2:cepstrumcutoff))').^2);
            DistExtended(j,i) = DistExtended(i,j);
        end
        disp(['Size 2^' num2str(log2(N)) ' Repetition: ' num2str(r) '/' num2str(repetitions)  ' Extended Cepstrum: ' num2str(i) '/' num2str(n)])
    end

    VectorExtended = squareform(DistExtended);
    Zextended = linkage(VectorExtended);
    % cextended = cluster(Zextended,'cutoff',0.45,'criterion','distance');
    cextended = cluster(Zextended,'maxclust',clusters);
    timeperformance(r,4) = toc;
    ARI(r,4) = adjrandindex(cextended,clusterid);

    disp('Start estimating transfer function models')

    tic
    systems = cell(n,1);
    opt = tfestOptions;
    opt.Advanced.StabilityThreshold.z = 1-10^(-2);
    for i = 1:n
        Y = iddata(outputs(i,:)',inputs(i,:)',1);
        systems{i} = tfest(Y,5,[],NaN,'Ts',Y.Ts,opt);
        disp(['Size 2^' num2str(log2(N)) ' Repetition: ' num2str(r) '/' num2str(repetitions)  ' Estimate model: ' num2str(i) '/' num2str(n)])
    end
    systemestimatetime = toc;


    disp('Start Calculating H2-norm Distance')
    tic
    DistH2 = zeros(n,n);

    for i = 1:n
        j = i+1;
        while j<=n

            DistH2(i,j) = norm(systems{i}-systems{j});
            DistH2(j,i) = DistH2(i,j);
            j = j+1;
        end
        disp(['Size 2^' num2str(log2(N)) ' Repetition: ' num2str(r) '/' num2str(repetitions)  ' H2: ' num2str(i) '/' num2str(n)])
    end

    VectorH2 = squareform(DistH2);
    ZH2 = linkage(VectorH2);
    % cH2 = cluster(ZH2,'cutoff',0.25,'criterion','distance');
    cH2 = cluster(ZH2,'maxclust',clusters);
    timeperformance(r,5) = toc + systemestimatetime;
    ARI(r,5) = adjrandindex(cH2,clusterid);
    
    
    disp('Start Calculating HInf-norm Distance')
    tic
    DistHInf = zeros(n,n);

    for i = 1:n
        j = i+1;
        while j<=n

            DistHInf(i,j) = norm(systems{i}-systems{j},Inf);
            DistHInf(j,i) = DistHInf(i,j);
            j = j+1;
        end
        disp(['Size 2^' num2str(log2(N)) ' Repetition: ' num2str(r) '/' num2str(repetitions)  ' HInf: ' num2str(i) '/' num2str(n)])
    end

    VectorHInf = squareform(DistHInf);
    ZHInf = linkage(VectorHInf);
    cHInf = cluster(ZHInf,'maxclust',clusters);
    timeperformance(r,6) = toc + systemestimatetime;
    ARI(r,6) = adjrandindex(cHInf,clusterid);
    
end
end
