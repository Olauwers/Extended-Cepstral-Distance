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
    cepstrumcutoff = floor(size(pwelch(inputs(1,:),[],[],'centered'),1)/2);
  
    disp('Start Euclidean Distance')
    
    tic
    DistEuclid = squareform(pdist(normoutputs));
    VectorEuclid = squareform(DistEuclid);
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
        ZLBKeogh = linkage(VectorLBKeogh);
        % cLBKeogh = cluster(ZLBKeogh,'cutoff',10^3,'criterion','distance');
        cLBKeogh = cluster(ZLBKeogh,'maxclust',clusters);
        timeperformance(r,2) = toc;
        ARI(r,2) = adjrandindex(cLBKeogh,clusterid);
    end
    
    disp('Start Original Cepstral Distance')
    tic
    OutputCeps = zeros(n,size(pwelch(outputs(1,:),[],[],'centered'),1));
    DistCepstral = zeros(n,n);
    weights = 0:1:cepstrumcutoff-1;

    for i = 1:n

        OutputCeps(i,:) = ifftshift(log(pwelch(outputs(i,:),[],[],'centered')));
    end

    clear i

    for i = 1:n

        for j = 1:i-1

            DistCepstral(i,j) = weights*((OutputCeps(i,1:cepstrumcutoff)' - OutputCeps(j,1:cepstrumcutoff)').^2);
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
    InputCeps = zeros(n,size(pwelch(inputs(1,:),[],[],'centered'),1));
    OutputCeps = zeros(n,size(pwelch(outputs(1,:),[],[],'centered'),1));
    DistExtended = zeros(n,n);
    weights = 1:1:cepstrumcutoff-1;

    for i = 1:n

        InputCeps(i,:) = ifftshift(log(pwelch(inputs(i,:),[],[],'centered')));
        OutputCeps(i,:) = ifftshift(log(pwelch(outputs(i,:),[],[],'centered')));
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