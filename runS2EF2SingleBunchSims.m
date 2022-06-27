datapath = append(pwd,'/matfiles/');
setname = 'scan01/';

N = 1; % Number of MC data points

rng('shuffle');

simvals = struct('P1',[-0.2 0.0], 'P2',[0 0.0], 'V1',[0 0.001], 'V2',[0 0.001], ...
              'qi',[0 0.01], 'dx',[0 0], 'dy',[0 0]); % Center vals, excursions

tic
for n = 1:N
    vals = structfun(@samplevals,simvals,'UniformOutput',false); % Randomly samples parameter space
    data = scanFunctionS2EF2SingleBunchSims(vals);
    data.scanvals = vals;
    filename = append(num2str(now,'%4.8f'),'.mat');
    filepath = append(datapath,setname,filename);
   % save(filepath,'data','-v7'); % Saves the simulation data
end
toc

function val = samplevals(in)
low = in(1) - in(2);
high = in(1) + in(2);
val = (high-low)*rand()+low;
end
