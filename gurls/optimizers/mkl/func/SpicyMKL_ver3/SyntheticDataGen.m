function [X,Y] = SyntheticDataGen(datatype,numsample,param)

switch(datatype)
 case 1
  [X,Y]=GenSine(numsample, param);
 case 2
  [X,Y]=GenHalfpipe(numsample, param);
 case 3
  [X,Y]=GenStripes(numsample, param);
end;

function [X,Y]=GenSine(numsample, param)
X = randn(numsample,2);
X(:,2) = X(:,2)*param.noise;
X(:,2) = X(:,2) + sin(X(:,1)*param.freq); 
Y = rand(numsample,1);
Y = (Y>0.5) - (Y<=0.5);
X(:,2) = X(:,2) + Y*param.dist;

function [X,Y]=GenHalfpipe(numsample, param)
pnum = ceil(numsample/2);
Y = [ones(pnum,1); -ones(numsample-pnum,1)];
pind = 1:pnum;
theta = pi*(rand(pnum,1)-0.5);
X(pind,1) = 0.5*sin(theta);
X(pind,2) = 0.5*cos(theta);
nind = (pnum+1):numsample;
nnum = length(nind); 
theta = pi*(rand(nnum,1)-0.5);
X(nind,1) = (0.5 - 2.5*param.noise)*sin(theta);
X(nind,2) = (0.5 - 2.5*param.noise)*cos(theta); 
noise = param.noise*randn(pnum,2);
X(pind,:) = X(pind,:) + noise;
noise = param.noise*randn(nnum,2)*0.5; 
X(nind,:) = X(nind,:) + noise;

function [X,Y]=GenStripes(numsample, param)
XX = rand(numsample,2);
kizami = 1/param.freq;
nori = kizami*0.2;
st = 0; ed = kizami;
X = XX;
sig = 1;
for ii=1:ceil(param.freq)
  ind = intersect(find(XX(:,1) <= ed),find(XX(:,1) > st));
  X(ind,1) = XX(ind,1) + ones(length(ind),1)*nori*(ii-1);
  st = ed; ed = st + kizami;
  Y(ind,1) = sig;
  sig = -sig;
end;
