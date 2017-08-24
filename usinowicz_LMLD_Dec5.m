%Notes Dec5
%Using this code mostly to mess around with approximations at the very end.
%They all give close fits. Two versions of the linearized, low-density
%version: simple quadratics (Ni=4*x+4 or 2*x+6, Nm=6*x+2), full
%approximation (Ni=a*x+p, Nm=b*x+p). There may still be some questions with
%Nr: does it work best with x, c*x+p-p. Is c based on 1 or 1/2? By eye,
%best fit right now is with full approximations, Nr=x. 


% This code is a cleaned up version of the most current code tht I am
% using to run the spatial lottery model -- so no stages. Look at the code
% in usinowicz_3stg_FDM.m for a non-spatial, stage-structure model. I've
% tried to make notes throughout to explain what things are. 

%The basic structure of this program is three large iterative loops. The
%two outermost loops have nothing to do with the mechanics of the model,
%they are simply for doing reps to build some stats, or to increment
%through certain parameters (in this program, the size of the starting
%cluster). You can basically ignore them for learning how the model works. 

%There are actually two model imbedded in this code -- the spatial (variable pop)
%and non-spatial (variable pop_lot) versions of the lottery model. The 
%non-spatial (classic) model is really simple and could probably be written in 
%about 12 lines of code. I've tried to comment around it to show you where it is. 

%Most of the code is for the spatial model. The most difficult aspect of
%this code to understand is probably the dispersal step. The easiest way to
%do this in Matlab is to make use of its built-in Fourier transform
%algorithms and to do a convultion of the matrix of individuals, with a
%matrix defining the dispersal kernel (in Fourier space). This is a very
%mathematically complicated thing, but a very easy tool to use in
%Matlab. If you want to learn more about it, I can dig up some web sites
%that I found initially that did a good job explaining this stuff. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%This is used for initial cluster size (1st loop). I have everything set
% to 1 when I am doing constant low-level background invasion (set by 
% the parameter "bet" short for "beta" in the literature. 

delinc= 1;
delri =1;
delre =1;


%This sets the reps of the 2nd loop -- used for stat building when
%wanted. 

runs = 1; 

%Also: Check below: just because this is set, does not mean any individuals are
%actually being place on the lattice. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation variables and species variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lattice size 
np=32;

%Time steps. Sets the size of 3rd loop. 
ngens=4000;  

%Species 1 (invader) repro rate: mean and variance. 
m1=1.0;
vp1=1;

%Species 2 (resident) repro rate: mean and variance. 
m2=1.01;
vp2=1;

%Correlation between species 1 and 2
cf12=0;

%Background invasion rate (1e-4 is my standard) 
bet=1e-3;

%Mortality rate
mu=0.1;

%Sometimes I use the survival rate instead, when I want to make it
%stochastic. "vs1" sets the variance in survival rate. 
survival=1-mu;
vs1=0;

%Once the size of 1st and 3rd loops are known, initiaize matrices to store
%each run of the simulation

pop_all=zeros(ngens+1,runs);
pop_lot_all=zeros(ngens+1,runs);

%Index variable for entering data into matrices
enter=1;

% 1st loop -- steps through different starting cluster sizes, if wanted. 
for cs = delri:delinc:delre,
    
% 2nd loop -- for stat building. Just multiple runs of simulation given a
% specific parameter set and initial cluster size
for its=1:runs,
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets simulation variables that need to be reset every run
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Initialize two square lattices
XY=zeros(np);
XYnxt=XY;    

% For keeping track of and displaying populations

%Lottery model with local dispersal (the full spatial model)
pop=zeros(ngens+1,1);

%Classic lottery model
pop_lot=pop; pop_lotB=pop;

%Nucleation difference equation
nde=pop;

%If doing stochastic mortality
s1=survival+sqrt(vs1)*randn(ngens+1,1);
    
%Keep track of average probability of invader
mpxy=pop;
mneigh=pop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section for generating random reproduction rates with the given mean,
%variance, and covariance structure. Current distribution type: 
%    LOGNORMAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
ms= [m1; m2];
vs = [vp1; vp2];
CorrMat =[1 cf12; cf12 1];
            
%%%%%% Transform these to LOGNORMAL means, variances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

MLog=log((ms.^2)./sqrt(vs+ms.^2));
VLog=log(vs./(ms.^2)+1);
cmv=sqrt(vs)*sqrt(vs)'.*CorrMat;
cmvLog=log(1+(cmv)./abs(ms*ms'));

%This matrix contains a random variable for each time step representing
%repdroduction for each species. 

psn=exp(mvnrnd(MLog,cmvLog,ngens+1));      



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section is for initial placement of occupied sites on the lattice.
%There are   ways to do this. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%First: number of indiviudals of species 1 to place

%ind=(np^2)-1; %Only 1 of invader (species 2)            
%ind= np^2/2;   %50% invade  
%ind = 100;

%To mimic process for metastable stuff (that is, with contstant background 
%invasion rate) make this 0% 

ind=0*np^2; % percentage based

%1)RANDOM placement 
for i=1:ind,
   ax=ceil(np*rand);
   ay=ceil(np*rand);
   XY(ax,ay)=1;
end;

%2) single, CENTERED cluster. Use this when incrementing through different
%initial cluster sizes:

% %cluster neighborhood order
% sqs=(1:2:15).^2;
% a=sqs(sqs>=cs);
% c_ord=a(1);
% cluster=zeros(c_ord,1);
% for i=1:cs,
%    cluster(i)=1;
% end;
% base=sqrt(c_ord);
% cluster=reshape(cluster, base, base);
% 
% %Center the cluster
%    
% strt=(np/2)-base;
% fin =(np/2)-1;
% for xx=strt:fin,
%   for yy=strt:fin
% 
%       XY(xx,yy)=cluster(xx-((np/2)-base)+1,yy-((np/2)-base)+1 );
%   end;
% end;   
   
   
%Initial population values for population matrices:

pop_lot(1)= 1e-4;
pop_lotB(1)= 1e-4;
nde(1) = 4;
%pop_lot(1)=sum(sum(XY))/np^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section gets into the dispersal mechanics. 
%In this code I only have  the simplest CA type dispersal defined. I have 
%code for more realistic kernels elsewhere (e.g. 2D Gaussian). That will be 
%something we can incorporate later on. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This number should be set manually. It should be odd. It is the 
%diameter of a square dispersal neighborhood. For example, rn =3 is a 3 by
%3 square with the dispersing individual at the center. Thus rn=3 defines
%dispersal to nereast neighbors, including corners.
rn=3;
rp=(rn-1)/2;
eq1=0.5;
%This matrix is the same size as the entire lattice and puts the dispersal
%neighborhood at the center. It is needed for the convolution. 
G2D=zeros (np);

% Set up the matrix 
strt=((np/2)+1)-(rn-1)/2;
fin =((np/2)+1)+(rn-1)/2;
for xx=strt:fin,
  for yy=strt:fin
      G2D(xx,yy)=1;
  end;
end;
G2D ((np/2)+1,(np/2)+1) =0;  %Don't count middle cell

%For normal neighborhoods:
sig1=rn^2 -1;

%for von Neumann neighborhood (no dispersal to the corners). 
%G2D(strt,strt)=0; G2D(strt,fin)=0; G2D(fin,strt)=0; G2D(fin,fin)=0;
%sig1=4;

%Take the Fast Fourier Transform of G2D. The convolution is performed in
%Fourier space.

fg2d=fft2(G2D);


%3rd loop. Currently setup with a while statement, but could be implemented
%with for loop. I am using the while just so that I can end a loop early
%under certain circumstances and save time when doing large numbers of
%reps. 

t=1;

   while ( t<=ngens)

        %Determine how many neigbors each cell has by performing a 2-D 
        %convolution with G2D in Fourier space. 

        fr1= fft2(XY);
        neigh = real(fftshift(ifft2(fr1.*fg2d)));

        % Lottery-style transition probabilities defined for each site on
        % the lattice 
        Pxy =(bet+ psn(t,1).*neigh)./(psn(t,1).*neigh + psn(t,2).*(sig1-neigh) + bet); 

        %Births: 
        % 1stMake a matrix of probabilities 
        r1=zeros(np);
        r2=r1;
        repro2=rand(np);

        %Assign births
        r1(repro2<=Pxy)=1;
        r2(repro2>=Pxy)=1;

        %Deaths:
        d1 =(XY.*rand(np)>=survival&XY.*rand(np)<1);
        d2 =((ones(np)-XY).*rand(np)>=survival& (ones(np)-XY).*rand(np)<1);
        dead=d1+d2;    
       
        % Subtract the dead and add the living: 
        % This creates the next generation population by 
        % 1: adding births to deaths, creating a matrix of 2s 
        % where birthsoccurred on recently vacated spaces.
        % 2: making sure this matrix includes only 2s and 0s
        % (clearing 1s). 
        % 3: subtract the dead (d1) as well as this matrix from current population (XY)
        % 4: Values of -2 will represent sites of new capture
        % for species 1. Convert all negatives to a value of 1. 
                
        nxt=r1+dead;
        nxt(nxt==1)=0;
        XYnxt=XY-d1-nxt;
        XYnxt(XYnxt==-2)=1;
                              
        %Summary population for lattice simulation
        pop(t)=sum(sum(XY));
        mneigh(t)=sum(sum(neigh.*XY));
        
        %Replace XY with XYnxt to move forward in time
        XY=XYnxt;
        
        
        mpxy(t)=mean(mean(Pxy(Pxy>0)));
        
%         % This part is for counting cluster transitions
%         if t>1,
%             if pop(t-1)+1<= max_c && pop(t)+1 <=max_c,
% 
%                 clust_data(pop(t-1)+1, pop(t)+1) = clust_data(pop(t-1)+1, pop(t)+1)+1;
% 
%             end;
%         end;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the classic (non-spatial) lottery model. The 2-species model
%effectively reduces to a 1-d model, since the population of species 2 is
%just 1-pop_species_1. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pop_lot(t+1)=s1(t)*pop_lot(t)+(1-s1(t))*(psn(t,1)*pop_lot(t)/(psn(t,1)*pop_lot(t)+psn(t,2)*(1-pop_lot(t))));
        pop_lotB(t+1)=s1(t)*pop_lotB(t)+(1-s1(t))*(psn(t,2)*pop_lotB(t)/(psn(t,2)*pop_lotB(t)+psn(t,1)*(1-pop_lotB(t))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the nucleation difference equation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n2=(pi*(sqrt(nde(t)/pi))+nde(t)/2);
n1=pi*(sqrt(nde(t)/pi)+rp*(rp+1)/2)+((rn-1)*nde(t))-eq1*(rn/np)^2*nde(t).^2;
mn=pi*(sqrt(nde(t)/pi)+rp*(rp+1)/2)+(rn*(rn-1)*nde(t))-(rn*(rn-1))/np^2*nde(t).^2;

kNi=(mn./n1);
%kNi2=(mn./n2);
kNi2=(mn./n2);
if(kNi2>(rn^2-1)) kNi2=(rn^2-1); end;

%So this only works if n1 is > the number used to calculate p1
p1=mu*kNi.*psn(t,1)./(psn(t,1).*kNi+psn(t,2).*(rn^2-kNi));
p2=mu*kNi2.*psn(t,2)./(psn(t,2)*kNi2+psn(t,1)*(rn^2-kNi2));
%n2=nde(t);
Bn=n1*p1;
Dn=n2*p2;
nde(t+1)=nde(t)+(Bn-Dn);

% If you want to see what is happening with XY uncomment this. However, it
% is very time-intensive for large numbers of time steps (ngens)
%         figure(5)
%         pcolor(XY)




 t=t+1; 
  

   end;
   
   
   
%Scale pop
pop=pop./np^2;

%Saving the populations for each run

pop_all(:,its)=pop;       
pop_lot_all(:,its)=pop_lot;

% If you want to see a figure at the end of each run. Turn this off if you
% are doing a large number of runs.

tt=1:ngens;
figure (1)
plot(tt,pop(1:ngens),tt,pop_lot(1:ngens));
legend( 'Simulation', 'Lottery Model');
    
end;

enter=enter+1;
end;

% Below I include code designed as a simple way to compare the relative competition experienced by a
% species trying to invade the other. The slope of the growth curve immediately
% after invasion initiates should reflect relative interspecific competition. 
% The exact region to use if you are going to measure this is somewhat subjective.
% It should be approximately linear between a near-zero lower bound, and an
% upper bound below any significant curviture. But as you'll see, with a stochastic
% model it is difficult to pick out exactly where growth stops being
% linear. As a general rule with the lottery model, though, an upper limit 
% seems to be about 0.1*equilibrium -- which is usually in the range 0.3 to 0.5. 

% Also, with the spatial model things are complicated slighly more
% by the fact that the lower limit to this linear region should be marked
% by the critical cluster size -- the population size at which growth rate
% is no longer controlled by the background invasion rate, beta (variable
% bet) but by the parameters of the models (i.e. repro rates, mortality
% rates). For the initial parameters in this model it is 24. So pop =
% 24/np^2; 

aveslp_pop=zeros(runs,1);
aveslp_rid=zeros(runs,1);

% This number (24) changes if other parameters of the model change

r1=psn(:,1);
r0=psn(:,2);


rcc=abs(mean((-sqrt((256*psn(:,1)+36*psn(:,2)).*(psn(:,1)+psn(:,2)))+20*psn(:,1)+26*psn(:,2))./(psn(:,1)+5*psn(:,2))));
ccrit=rcc./np^2;

for k=1:runs,
    % Find Ccrit in the spatial model. With the & comparison is more time
    % consuming, but is necessary to account for multiple peaks. 
        apB=find (pop_all(:,k) <=ccrit & pop_all(:,k) >=ccrit-0.5*ccrit, 1, 'last');
    
    % Find the time at which population reach 0.1*equilibrium in the
    % spatial model
        apE=find(pop_all(:,k) >=.1*mean(pop_all(:,k)) & pop_all(:,k) <=.1*mean(pop_all(:,k))+0.05, 1, 'last');
        
    % Find the time at which population reach 0.1*equilibrium in the
    % classic model
        alB=find (pop_lot_all(:,k) <=ccrit & pop_lot_all(:,k) >=ccrit-0.5*ccrit, 1, 'last');
        alE=find (pop_lot_all(:,k) >=.1*mean(pop_lot_all(:,k)) & pop_lot_all(:,k) <=.1*mean(pop_lot_all(:,k))+0.05, 1, 'last');
        
    % Fit a straight line through this region for the spatial model
        ps=polyfit(1:((apE-apB)+1),log(pop_all(apB:apE, k)'),1);
        
    % Fit a straight line through this region for the classic model
        pr=polyfit(1:((alE-alB)+1),log(pop_lot_all(alB:alE, k)'),1);
        
        aveslp_pop(k)=ps(1);
        aveslp_rid(k)=pr(1);
end

%This is the analytical approximation for the classic lottery model from 
% e.g. Chesson and Warner 1981
chesson81=mean(log(1+mu*(psn(:,1)./psn(:,2)-1)));
chesson81b=mean(log(1+mu*(psn(:,1)-psn(:,2))./psn(:,2)));

cc=600;
%Approximation from nucleation/binomial distribution 
%Note: These all seem to work, to varying degrees. The first actually works
%the best with kNi=kNi2. This is because sometimes kNi underestimates the 
%ratio in favor of the resident. 
%The second one works well at the critical radius, if the critical radius 
%is actually large enough such that kNi=kNi2=3 (in the case of sig1=9)
%The third one tends to overestimate slightly for low mu, but becomes more
%accurate as mu gets larger. My guess is that clusters are less coherent as
%mu gets larger, translating into higher-than-average values of kNi and
%kNi2, which are both assumed to be 1/2*sig1 in this approximation. 
ic1=1;
me13=zeros(cc,1); me13b=me13; me13D=me13; me13c1=me13;
me13g=zeros(cc,1); me13bg=me13; me13Dg=me13; me13c1g=me13;
me13(1)=1; me13b(1)=ic1;me13c1(1)=1; 

%Use these for diagnostics
pdif=zeros(cc,1); pdif2=pdif;
s1=pdif; s2=pdif; s3=pdif; s4=pdif; s5=pdif; s6=pdif; s7=pdif; s8=pdif;

nid=zeros(cc,1);
mnd=nid; nrd=nid;
% for i=1:cc
% nid(i)=max(edges_all(pop_all==i));
% nrd(i)=mean(edgesphi_all(pop_all==i));
% mnd(i)=max(mneigh_all(pop_all==i));
% end;

phi=9;
x_mid=4;
a=(rn+sqrt(pi)/(2*sqrt(x_mid)));
b=(rn*(rn-1)+sqrt(pi)/(2*sqrt(x_mid)));
c=(1/2+sqrt(pi)/(2*sqrt(x_mid)));
p=pi*(rp*(rp+1)/2)+sqrt(pi*x_mid)-sqrt(pi)*x_mid/(2*sqrt(x_mid));


px=0;
nw1=(rn-1)*rn/np^2;
Bnpx=mu*(psn(:,1)*(nw1*px+rn^2).*(a*px+p))./(psn(:,1).*(nw1*px^2+rn^2*px)+psn(:,2).*(-nw1*px^2+(rn^3-rn^2)*px+p*rn^2));
Dnpx=mu*(psn(:,2)*(nw1*px+rn^2))./(psn(:,2).*(nw1*px+rn^2)-nw1.*psn(:,1)*px);
me13px=mean(log(1+(Bnpx-Dnpx)));

Bnpx2=mu*(psn(:,1)*rn^2*(rn+p))./(psn(:,2)*(rn^3+p*rn^2-rn^2)+psn(:,1)*rn^2);
Dnpx2=mu;
me13px2=mean(log(1+(Bnpx2-Dnpx2)));

for (ic=1:cc)
% %np=64;for (ic=1:cc)
% %np=64;
% %Benchmark
% %With Data: 
% nnmd=mnd(ic);
% nnid=nnmd./nid(ic);
% nnrd=nnmd./ic;%nrd(ic);
% BnD=mu*nid(ic)*nnid.*psn(:,1)./(psn(:,1).*nnid+psn(:,2).*(rn^2-nnid));
% DnD=mu*ic*nnrd.*psn(:,2)./(psn(:,2).*nnrd+psn(:,1).*(rn^2-nnrd));
% %Bn1=mu*(p+a*ic)*nni.*psn(:,1)./(psn(:,1).*nni+psn(:,2).*(rn^2-nni));
% %Dn1=mu*(p-pi+c*ic)*nnr.*psn(:,2)./(psn(:,2).*nnr+psn(:,1).*(rn^2-nnr));
% % Bn1=mu*(p+rn*ic).*r1./(r1+r0.*(rn^2/nni-1));
% % Dn1=mu*(ic).*r0./(r0+r1.*(rn^2/nnr-1));
% me13D(ic)=mean(log(1+(BnD-DnD)/ic));

c1f=(phi-1)-(pi*(sqrt(1/pi)+rp*(rp+1)/2)+((rn))-eq1*(rn/np)^2);
c2f=(phi-1)-(pi*(sqrt(1/pi)+rp*(rp+1)/2)+(rn*(rn-1))-(rn*(rn-1))/np^2);
n2=me13(ic); %(pi*(sqrt(me13(ic)/pi))+me13(ic)/2);
n1=c1f+pi*(sqrt(me13(ic)/pi)+rp*(rp+1)/2)+((rn)*me13(ic))-eq1*(rn/np)^2*me13(ic).^2;
mn=c2f+pi*(sqrt(me13(ic)/pi)+rp*(rp+1)/2)+(rn*(rn-1)*me13(ic))-(rn*(rn-1))/np^2*me13(ic).^2;

%if(mn>(rn^2-1)*me13(ic)) mn=(rn^2-1)*me13(ic); end;
kNi=(mn./n1);
%kNi2=(mn./n2);
kNi2=(mn./n2);
if(kNi2>(rn^2-1)) kNi2=(rn^2-1); end;
%So this only works if n1 is > the number used to calculate p1
p1=mu*kNi.*psn(:,1)./(psn(:,1).*kNi+psn(:,2).*(rn^2-kNi));
p2=mu*kNi2.*psn(:,2)./(psn(:,2)*kNi2+psn(:,1)*(rn^2-kNi2));

%n2=me13(ic);
Bn=n1*p1;
Dn=n2*p2;
pdif(ic)=mean(Bn-Dn);
me13(ic+1)=mean(me13(ic)+(Bn-Dn));

%Drop quadratics in Ni and Nr
a=(rn+sqrt(pi)./(2*sqrt(me13b(ic))));
b=(rn*(rn-1)+sqrt(pi)./(2*sqrt(me13b(ic))));
%c=(1/2+sqrt(pi)./(2*sqrt(me13b(ic))));
p=pi*(rp*(rp+1)/2)+sqrt(pi*me13b(ic))-sqrt(pi)*me13b(ic)/(2*sqrt(me13b(ic)));
c2=(phi-1)-b-p;
c1=(phi-1)-a-p;

%nnm=6*me13b(ic)+2;
nnm=(c2+p+b.*me13b(ic));
%nni=2*me13b(ic)+6;
nni=(c1+p+a.*me13b(ic));
nnr=me13b(ic); 
%nnr=(p-pi+c.*ic);
anni=nnm/nni;
annr=nnm/nnr;
%if(nnr >(rn^2-1)*ic) nnr=(rn^2-1)*ic; end; 
%if(nni >(rn^2-1)*ic) nni=(rn^2-1)*ic; end; 
%nnr=kNi2(ic); 
%pdif(ic)=mean((anni.*psn(:,1)./(psn(:,1).*anni+psn(:,2).*(rn^2-anni)))-(annr.*psn(:,2)./(psn(:,2).*annr+psn(:,1).*(rn^2-annr))));

Bn1=mu*anni*nni.*psn(:,1)./(psn(:,1).*anni+psn(:,2).*(rn^2-anni));
Dn1=mu*annr*nnr.*psn(:,2)./(psn(:,2).*annr+psn(:,1).*(rn^2-annr));
me13b(ic+1)=(me13b(ic)+(mean(Bn1)-mean(Dn1)));

%Checkinuc3a=mean(log(1+mu*((r1./r0).*(((f*rcc+h)*(g*rcc+h))./(rcc*(phi*(f*rcc+h)+(g*rcc+h)*((r1./r0)-1))))-(g*rcc+h)./((g*rcc+h)*(1-(r1./r0))+phi*rcc*(r1./r0)))));ng some algebra. This should be (roughly) the same as the
%non-quadratics 
nw1=(rn-1)*rn/np^2;
Bnc1=mu*(psn(:,1)*me13c1(ic).*(nw1*me13c1(ic)+rn^2).*(a*me13c1(ic)+p))./(psn(:,1).*(nw1*me13c1(ic)^2+rn^2*me13c1(ic))+psn(:,2).*(-nw1*me13c1(ic)^2+(rn^3-rn^2)*me13c1(ic)+p*rn^2));
Dnc1=mu*(psn(:,2)*me13c1(ic).*(nw1*me13c1(ic)+rn^2))./(psn(:,2).*(nw1*me13c1(ic)+rn^2)-nw1.*psn(:,1)*me13c1(ic));
me13c1(ic+1)=mean(me13c1(ic)+(Bn1-Dn1));


end;
mind=1;
for (ic=1:cc)
% %np=64;for (ic=1:cc)
% %np=64;
% %Benchmark
% %With Data: 
% nnmd=mnd(ic);
% nnid=nnmd./nid(ic);
% nnrd=nnmd./ic;%nrd(ic);
% BnD=mu*nid(ic)*nnid.*psn(:,1)./(psn(:,1).*nnid+psn(:,2).*(rn^2-nnid));
% DnD=mu*ic*nnrd.*psn(:,2)./(psn(:,2).*nnrd+psn(:,1).*(rn^2-nnrd));
% %Bn1=mu*(p+a*ic)*nni.*psn(:,1)./(psn(:,1).*nni+psn(:,2).*(rn^2-nni));
% %Dn1=mu*(p-pi+c*ic)*nnr.*psn(:,2)./(psn(:,2).*nnr+psn(:,1).*(rn^2-nnr));
% % Bn1=mu*(p+rn*ic).*r1./(r1+r0.*(rn^2/nni-1));
% % Dn1=mu*(ic).*r0./(r0+r1.*(rn^2/nnr-1));
% me13D(ic)=mean(log(1+(BnD-DnD)/ic));
ic=ic*mind;
c1f=(phi-1)-(pi*(sqrt(1/pi)+rp*(rp+1)/2)+((rn))-eq1*(rn/np)^2);
c2f=(phi-1)-(pi*(sqrt(1/pi)+rp*(rp+1)/2)+(rn*(rn-1))-(rn*(rn-1))/np^2);
n2=ic;%(pi*(sqrt(ic/pi))+ic/2);
n1=c1f+pi*(sqrt(ic/pi)+rp*(rp+1)/2)+((rn-1)*ic)-eq1*(rn/np)^2*ic.^2;
mn=c2f+pi*(sqrt(ic/pi)+rp*(rp+1)/2)+(rn*(rn-1)*ic)-(rn*(rn-1))/np^2*ic.^2;
%if(mn>(rn^2-1)*ic) mn=(rn^2-1)*ic; end;
kNi=(mn./n1);
%kNi2=(mn./n2);
kNi2=(mn./n2);
if(kNi2>(rn^2-1)) kNi2=(rn^2-1); end;
%So this only works if n1 is > the number used to calculate p1
p1=mu*kNi.*psn(:,1)./(psn(:,1).*kNi+psn(:,2).*(rn^2-kNi));
p2=mu*kNi2.*psn(:,2)./(psn(:,2)*kNi2+psn(:,1)*(rn^2-kNi2));

%n2=ic;
Bn=n1*p1;
Dn=n2*p2;
s1(ic)=mean(p1);
s2(ic)=mean(p2);
pdif(ic)=mean(Bn-Dn);
%me13(ceil(ic/mind+1))=mean(log(1+(Bn-Dn)/ic));
me13g(ceil(ic/mind+1))=(log(1+(mean(Bn/ic)-mean(Dn/ic))));


%Drop quadratics in Ni and Nr
a=(rn+sqrt(pi)./(2*sqrt(ic)));
b=(rn*(rn-1)+sqrt(pi)./(2*sqrt(ic)));
%c=(1/2+sqrt(pi)./(2*sqrt(ic)));mean(Bn-Dn)
p=pi*(rp*(rp+1)/2)+sqrt(pi*ic)-sqrt(pi)*ic/(2*sqrt(ic));
c2=(phi-1)-b-p;
c1=(phi-1)-a-p;

%nnm=6*ic+2;
nnm=(c2+p+b.*ic);
%nni=2*ic+6;
nni=(c1+p+a.*ic);
nnr=ic; 
%nnr=(p-pi+c.*ic);
anni=nnm/nni;
annr=nnm/nnr;
%if(nnr >(rn^2-1)*ic) nnr=(rn^2-1)*ic; end; 
%if(nni >(rn^2-1)*ic) nni=(rn^2-1)*ic; end; 
%nnr=kNi2(ic); 
%p1=anni.*psn(:,1)./(psn(:,1).*anni+psn(:,2).*(rn^2-anni));
%p2=annr.*psn(:,2)./(psn(:,2).*annr+psn(:,1).*(rn^2-annr));
p1=mu*anni./((psn(:,1)./psn(:,2)).*anni+(rn^2-anni));
p2=mu.*annr./(annr+(psn(:,1)./psn(:,2)).*(rn^2-annr));
%pdif(ic)=mean(p1-p2);
Bn1=mu*anni*nni.*psn(:,1)./(psn(:,1).*anni+psn(:,2).*(rn^2-anni));
Dn1=mu*annr*nnr.*psn(:,2)./(psn(:,2).*annr+psn(:,1).*(rn^2-annr));
s3(ic)=mean(p1);
s4(ic)=mean(p2);
pdif2(ic)=mean(Bn1-Dn1);
me13bg(ceil(ic/mind+1))=mean(log(1+(Bn1-Dn1)/ic));

%Checking some algebra. This should be (roughly) the same as the
%non-quadratics 
nw1=(rn-1)*rn/np^2;
Bnc1=mu*(psn(:,1)*ic.*(nw1*ic+rn^2).*(a*ic+p))./(psn(:,1).*(nw1*ic^2+rn^2*ic)+psn(:,2).*(-nw1*ic^2+(rn^3-rn^2)*ic+p*rn^2));
Dnc1=mu*(psn(:,2)*ic.*(nw1*ic+rn^2))./(psn(:,2).*(nw1*ic+rn^2)-nw1.*psn(:,1)*ic);
me13c1g(ceil(ic/mind+1))=mean(log(1+(Bnc1-Dnc1)/ic));

end;

%bndn=cc^(-1/2)/sqrt(3.14)*(psn(:,1)./(psn(:,1)+psn(:,2)*2)-psn(:,2)./(psn(:,2)+psn(:,1)*2))+psn(:,1)./(cc*(psn(:,1)+psn(:,2)*2));
ic=rcc;

a=(rn+sqrt(pi)./(2*sqrt((ic))));
b=(rn*(rn-1)+sqrt(pi)./(2*sqrt((ic))));
p=pi*(rp*(rp+1)/2)+sqrt(pi*ic)-sqrt(pi)*ic/(2*sqrt(ic));
f=a;
g=b;
h=p;
c2=(phi-1)-g-h;
c1=(phi-1)-f-h;

phi=rn^2;
ld1=mean(log(1+mu*((psn(:,1)./psn(:,2)).*(phi./( (phi-1)+(psn(:,1)./psn(:,2))))-1)));
pw1=mean(log(1+mu*(psn(:,1)./psn(:,2).*((phi-1)./(phi-2+psn(:,1)./psn(:,2)))-1)))+mean(log(1+mu*(psn(:,1)./psn(:,2).*((phi-1)./(phi-2+psn(:,1)./psn(:,2)))-1)))^2;
%nuc1=mean(log(1+mu*(psn(:,1)./psn(:,2).*((rn+b)./(rn+b-1+psn(:,1)./psn(:,2)))-1)));
nuc1=mean(log(1+(mu*r1*(rcc*rn+p).*(rn^2+nw1*rcc))./(r0.*(rcc*(rn^3-rn^2)+p.*rn^2-nw1*rcc^2)+r1.*(rcc*rn^2+nw1*rcc^2))-(mu*r0.*(rn^2+nw1*rcc))./(r0.*(rn^2+nw1*rcc)-nw1*r1.*rcc)));
nuc1a=mean(log(1+mu*((r1./r0).*((rcc*rn+p).*(phi+nw1*rcc)./(phi*(rcc*rn+p)+rcc*(nw1*rcc+phi)*((r1./r0)-1)))-((rcc*nw1+phi)./(phi+nw1*rcc*(1-(r1./r0)))))));
nuc2=mean(log(1+mu*(((r1./r0)*(rcc*rn+p))./((rcc*rn+p)+rcc*(r1./r0-1))-1)));
nuc3=mean(log(1+mu*((r1.*(f*rcc+h)*(g*rcc+h))./(rcc*(r0.*((f*phi-g)*rcc+h*phi-h)+r1.*(g*rcc+h)))-(r0.*(g*rcc+h))./(r1.*((phi-g)*rcc-h)+r0.*(g*rcc+h)))));
nuc3a=mean(log(1+mu*((r1./r0).*(((f*rcc+h)*(g*rcc+h))./(rcc*(phi*(f*rcc+h)+(g*rcc+h)*((r1./r0)-1))))-(g*rcc+h)./((g*rcc+h)*(1-(r1./r0))+phi*rcc*(r1./r0)))));
Dn=((f*rcc+h)*(g*rcc+h))./(rcc*(phi*(f*rcc+h)+(g*rcc+h)*((r1./r0)-1)))-r0./r1.*(g*rcc+h)./((g*rcc+h)*(1-(r1./r0))+phi*rcc*(r1./r0))+r0./r1;
nuc3b=mean(log(1+mu*((r1./r0).*Dn-1)));
nuc3_new=mean(log(1+mu*((r1./r0).*(((f*rcc+h+c1)*(g*rcc+h+c2))./(rcc*(phi*(f*rcc+h+c1)+(g*rcc+h+c2)*((r1./r0)-1))))-(g*rcc+h+c2)./((g*rcc+h+c2)*(1-(r1./r0))+phi*rcc*(r1./r0)))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For the critical cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1=psn(:,1);
r0=psn(:,2);

ic=50;
a=(rn+sqrt(pi)./(2*sqrt((ic))));
b=(rn*(rn-1)+sqrt(pi)./(2*sqrt((ic))));
p=pi*(rp*(rp+1)/2)+sqrt(pi*ic)-sqrt(pi)*ic/(2*sqrt(ic));
f=a;
g=b;
h=p;
c1=(phi-1)-f-h;
c2=(phi-1)-g-h;
%h adjusted 
h1=h+c2;
h=h+c1;
%Survival
s1=1;%(mu);

%Cubic root

%First round, but I think that these have been messed up. 
% a=((f*g*phi-f*g^2)*r1.^2+(f-1)*g^2.*r0.*r1+(g^2-f*g*phi)*r0.^2);
% %b=((f*h1+g*h-mu*g)*phi-2*f*g*h1-g^2*h+mu*g^2)*r1.^2+(-mu*f*phi^2+(mu*f+mu)*g*phi+(2*f-2)*g*h1+g^2*h-2*mu*g^2).*r0.*r1...
% %    +((-f*h1-g*h-mu*f*g)*phi+2*g*h1+mu*g^2)*r0.^2;
% b=(((f*h1+g*h-g*s1)*phi-2*f*g*h1-g^2*h+s1*g^2)*r1.^2+(-f*phi^2+(f+1)*s1*g*phi+(2*f-2)*g*h1+g^2*h-2*s1*g^2).*r0.*r1+((-f*h1-g*h-s1*f*g)*phi+2*g*h1+s1*g^2)*r0.^2);
% %c= (((h-mu)*h1*phi-f*h1^2+(2*mu*g-2*g*h)*h1)*r1.^2+(-mu*h*phi^2+((mu*f+mu)*h1+mu*g*h)*phi+(f-1)*h1^2+(2*g*h-4*mu*g)*h1).*r0.*r1...
% %    +(((-h-mu*f)*h1-mu*g*h)*phi+h1^2+2*mu*g*h1)*r0.^2);
% c=(((h-s1)*h1*phi-f*h1^2+(2*s1*g-2*g*h)*h1)*r1.^2+(-h*phi^2+((f+1)*h1+g*h)*s1*phi+(f-1)*h1^2+(2*g*h-4*s1*g)*h1).*r0.*r1+(((-h-s1*f)*h1-g*s1*h)*phi+h1^2+2*s1*g*h1)*r0.^2);
% %d=(mu-h)*h1^2*r1.^2+(mu*h*h1*phi+(h-2*mu)*h1^2).*r0.*r1+(mu*h1^2-mu*h*h1*phi)*r0.^2;
% d=(s1-h)*h1^2*r1.^2+(s1*h*h1*phi+(h-2*s1)*h1^2).*r0.*r1+(h1^2-h*h1*phi)*s1*r0.^2;
%After 1/x

%This is with the term s1, after dividing Bn/x and Dn/x. This gives
%imaginary results, so maybe not correct? 
% a=(((g^2-g*phi)*r1.^2+(-f*phi^2+(f+1)*g*phi-2*g^2).*r0.*r1+(g^2-f*g*phi)*r0.^2)*s1+(f*g*phi-f*g^2)*r1.^2+(f-1)*g^2.*r0.*r1+(g^2-f*g*phi)*r0.^2);
%   a=((g^2-g*phi)*r1.^2+(-f*phi^2+(f+1)*g*phi-2*g^2).*r0.*r1+(g^2-f*g*phi)*r0.^2)*s1+(f*g*phi-f*g^2)*r1.^2+(f-1)*g^2.*r0.*r1+(g^2-f*g*phi)*r0.^2;
% % b=(((2*g*h1-h1*phi).*r1.^2+(-h*phi^2+((f+1)*h1+g*h)*phi-4*g*h1).*r0.*r1+((-f*h1-g*h)*phi+2*g*h1).*r0.^2)*s1...
% %     +((f*h1+g*h)*phi-2*f*g*h1-g^2*h).*r1.^2+((2*f-2)*g*h1+g^2*h).*r0.*r1+((-f*h1-g*h)*phi+2*g*h1).*r0.^2);
%   b=(((2*g*h1-h1*phi)*r1.^2+(-h*phi^2+((f+1)*h1+g*h)*phi-4*g*h1).*r0.*r1+((-f*h1-g*h)*phi+2*g*h1)*r0.^2)*s1...
%       +((f*h1+g*h)*phi-2*f*g*h1-g^2*h)*r1.^2+((2*f-2)*g*h1+g^2*h).*r0.*r1+((-f*h1-g*h)*phi+2*g*h1)*r0.^2);
% %c=(h1^2*r1.^2+(h*h1*phi-2*h1^2).*r0.*r1+(h1^2-h*h1*phi).*r0.^2)*s1+(h*h1*phi-f*h1^2-2*g*h*h1).*r1.^2+((f-1)*h1^2+2*g*h*h1).*r0.*r1+(h1^2-h*h1*phi).*r0.^2;
%   c=((h1^2*r1.^2+(h*h1*phi-2*h1^2).*r0.*r1+(h1^2-h*h1*phi)*r0.^2)*s1+(h*h1*phi-f*h1^2-2*g*h*h1)*r1.^2+((f-1)*h1^2+2*g*h*h1).*r0.*r1+(h1^2-h*h1*phi)*r0.^2);
% % d=h*h1^2*r1.^2+h*h1^2.*r0.*r1;
%   d=h*h1^2*r1.^2+h*h1^2.*r0.*r1;

%This is with the term s1, so solving Bn-Dn=s1
a=((f*g*phi-f*g^2).*r1.^2+(f-1)*g^2.*r0.*r1+(g^2-f*g*phi).*r0.^2);
b=((g^2-g*phi)*r1.^2+(-f*phi^2+(f+1)*g*phi-2*g^2).*r0.*r1+(g^2-f*g*phi).*r0.^2)*s1+((f*h1+g*h)*phi-2*f*g*h1-g^2*h).*r1.^2+((2*f-2)*g*h1+g^2*h).*r0.*r1+((-f*h1-g*h)*phi+2*g*h1).*r0.^2;
c=(((2*g*h1-h1*phi).*r1.^2+(-h*phi^2+((f+1)*h1+g*h)*phi-4*g*h1).*r0.*r1+((-f*h1-g*h)*phi+2*g*h1).*r0.^2)*s1+(h*h1*phi-f*h1^2-2*g*h*h1).*r1.^2+((f-1)*h1^2+2*g*h*h1).*r0.*r1+(h1^2-h*h1*phi).*r0.^2);
d=(h1^2.*r1.^2+(h*h1*phi-2*h1^2).*r0.*r1+(h1^2-h*h1*phi).*r0.^2)*s1-h*h1^2.*r1.^2+h*h1^2.*r0.*r1;


r3= ((-b.^3./(27*a.^3)+(b.*c)./(6*a.^2)-d./(2*a))+sqrt((-b.^3./(27*a.^3)+(b.*c)./(6*a.^2)-d./(2*a)).^2+(c./(3*a)-b.^2./(9*a.^2)).^3)).^(1/3)...
    +((-b.^3./(27*a.^3)+(b.*c)./(6*a.^2)-d./(2*a))-sqrt((-b.^3./(27*a.^3)+(b.*c)./(6*a.^2)-d./(2*a)).^2+(c./(3*a)-b.^2./(9*a.^2)).^3)).^(1/3)...
    -b./(3*a);

p=-b./(3.*a); q=p.^3+(b.*c-3.*a.*d)./(6.*a.^2); r=c./(3.*a);
r3b=(q+(q.^2+(r-p.^2).^3).^(1/2)).^(1/3)+(q-(q.^2+(r-p.^2).^3).^(1/2)).^(1/3)+p;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%For a good general check on the behavior of Dn:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rcc=mean(r3);
ic=rcc;
a=(rn+sqrt(pi)./(2*sqrt((ic))));
b=(rn*(rn-1)+sqrt(pi)./(2*sqrt((ic))));
p=pi*(rp*(rp+1)/2)+sqrt(pi*ic)-sqrt(pi)*ic/(2*sqrt(ic));
f=a;
g=b;
h=p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These constants are a new feature. 
%They should be implemented for everything.
%When I have more time to blow. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c2=(phi-1)-g-h;
c1=(phi-1)-f-h;

a=0.01:0.01:10;
DnA=((f*rcc+h+c1)*(g*rcc+h+c2))./(rcc*(phi*(f*rcc+h+c1)+(g*rcc+h+c2)*((a)-1)))...
    -(1./a).*((g*rcc+h+c2)./((g*rcc+h+c2)*(1-(a))+phi*rcc*(a))-1);

%Versus the full version 
n2=ic; %(pi*(sqrt(ic/pi))+ic/2);
c1f=(phi-1)-(pi*(sqrt(1/pi)+rp*(rp+1)/2)+((rn))-eq1*(rn/np)^2);
c2f=(phi-1)-(pi*(sqrt(1/pi)+rp*(rp+1)/2)+(rn*(rn-1))-(rn*(rn-1))/np^2);
n1=c1f+pi*(sqrt(ic/pi)+rp*(rp+1)/2)+((rn)*ic)-eq1*(rn/np)^2*ic.^2;
mn=c2f+pi*(sqrt(ic/pi)+rp*(rp+1)/2)+(rn*(rn-1)*ic)-(rn*(rn-1))/np^2*ic.^2;
kNi=(mn./n1);
kNi2=(mn./n2);
%if(kNi2>(rn^2-1)) kNi2=(rn^2-1); end;
%So this only works if n1 is > the number used to calculate p1
p1=kNi./(a.*kNi+(rn^2-kNi));
p2=kNi2./(kNi2+a.*(rn^2-kNi2));
% p1a=1./(a+(phi/kNi-1));
% p2a=(1./a)./(1+a.*(phi/kNi2-1));
p1a=1./(a+(phi/((g*rcc+(phi-1))/(f*rcc+(phi-1)))-1));
p2a=(1./a)./(1+a.*(phi/((g*rcc+h+c2)/(rcc))-1));
%n2=ic; (1./a)
Bn=n1*p1/ic;
Dn=(1./a).*(n2*p2/ic-1);
DnB=(Bn-Dn);

%Versus LD 
ldn=phi./((phi-1)+a);

%Verus PW
pwn=(phi-1)./((phi-2)+a);

rcc=mean(r3);
ic=rcc;
a=(rn+sqrt(pi)./(2*sqrt((ic))));
b=(rn*(rn-1)+sqrt(pi)./(2*sqrt((ic))));
p=pi*(rp*(rp+1)/2)+sqrt(pi*ic)-sqrt(pi)*ic/(2*sqrt(ic));
f=a;
g=b;
h=p;

c2=(phi-1)-g-h;
c1=(phi-1)-f-h;

a=0.01:0.01:10;

%This is an approximate value for rp=1, based on c=1, f=rn, g=rn(rn-1),h=rn
h1=h+c2;
h=h+c1;

DnC=(((a-1)*h+(a-1)*rn+a-1)*h1^2+(((a-1)*rn^2+(2-2*a)*rn)*h+(a-1)*rn^3+(1-a)*rn^2+(2-2*a)*rn)*h1+(-rn^4+(1-a)*rn^3+(a-1)*rn^2)*h-rn^5+(1-a)*rn^4+(a-1)*rn^2)./((a.^2-2*a+1)*h1^2+((a-1)*rn^2*h+(a-1)*rn^3+(a.^2-3*a+2)*rn^2+(-2*a.^2+4*a-2)*rn)*h1+((1-a)*rn^3-rn^4)*h-rn^5+(2-2*a)*rn^4+(-a.^2+3*a-2)*rn^3+(a.^2-2*a+1)*rn^2);
