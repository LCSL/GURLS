function [ finalcode ] = generatesimplexcode(T)

%T is number of classes is the number of vectors in dimension T-1
% code est l'ensemble des codes du simplex jusqu'à la dimension T-1
%copyright Prof de Zacharie

code{1}=[-1;1];
code{2}=[1 0; -0.5 sqrt(3)/2;-0.5 -sqrt(3)/2]';

for dim =3:T-1 
sinalpha= sqrt(1-1/(dim)^2);
previous= code{dim-1};
current= zeros(dim,dim+1);
current(1,1)=1;
current(1,2:end)=-1/dim *ones(1,dim);
current(2:end,2:end)=previous*sinalpha;
code{dim}=current;
clear current;
clear previous;
end 
finalcode=code{T-1};
end