function [timeperformancetotal,ARItotal,orderofdistances,outputinfo,sizes] = runExperiments(repetitions)

sizes = [2^6 2^8 2^10 2^12 2^14 2^16];
timeperformancetotal = zeros(6,repetitions,6);
ARItotal = zeros(6,repetitions,6);

[timeperformancetotal(1,:,:),ARItotal(1,:,:),orderofdistances,outputinfo{1},LBKeogh] = clusterperformance(sizes(1),repetitions,true);
[timeperformancetotal(2,:,:),ARItotal(2,:,:),orderofdistances,outputinfo{1},LBKeogh] = clusterperformance(sizes(2),repetitions,true);
[timeperformancetotal(3,:,:),ARItotal(3,:,:),orderofdistances,outputinfo{1},LBKeogh] = clusterperformance(sizes(3),repetitions,true);
[timeperformancetotal(4,:,:),ARItotal(4,:,:),orderofdistances,outputinfo{1},LBKeogh] = clusterperformance(sizes(4),repetitions,true);
[timeperformancetotal(5,:,:),ARItotal(5,:,:),orderofdistances,outputinfo{1},LBKeogh] = clusterperformance(sizes(5),repetitions,true);
[timeperformancetotal(6,:,:),ARItotal(6,:,:),orderofdistances,outputinfo{1},LBKeogh] = clusterperformance(sizes(6),repetitions,true);

end
