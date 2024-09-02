function [pos_final,z_final] = DLGWO(param , model ,fobj)
it=1;
np=param.nPop;
nx=model.nVar;
maxit=param.MaxIt;
pp_pbest=zeros(np,nx); 
pv=zeros(np,nx);
optimal_pos=zeros(1,nx);
z=zeros(np);
z_pbest=zeros(np);
z_iter=zeros(maxit);
z_alpha=inf;
z_beta=inf;
z_delta=inf;
elitism=1; 
k_max=0.9; 
k_min=0.4; 
c=0;
d=1;
x=rand;
lb= c.*ones( 1,nx);    
ub= d.*ones( 1,nx);    
varmax=ub.*ones(1,nx);
varmin=lb.*ones(1,nx);
limvel=0.1;  
velmax=limvel.*(varmax(1,1:nx)-varmin(1,1:nx)); 
velmin=-velmax;
pp=ones(np,1)*(ub-lb).*rand(np,nx)+ones(np,1)*lb;
pv=ones(np,1)*(velmax-velmin).*rand(np,nx)+ones(np,1)*velmin;
alpha=zeros(1,nx);
beta=zeros(1,nx);
delta=zeros(1,nx);
gv_alpha=zeros(np,nx);
gv_beta=zeros(np,nx);
gv_delta=zeros(np,nx);
SolbestX=struct(); 
for j=1:np
%     pp(j,1:nx)=(varmax-varmin).*rand(1,nx)+varmin;
    gv_alpha(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
    gv_beta(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
    gv_delta(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
end
% Start the optimization process

% Objective function evaluations and determine the personal best solutions
% and objectives, if elitism is applied
for j=1:np
    [z(j),Sol]=fobj(pp(j,:));
    
    % Elitism
    if elitism==1
        z_pbest(j)=z(j);
        pp_pbest(j,1:nx)=pp(j,1:nx);
    end
end
for j=1:np
    if z(j)<=z_alpha
        z_alpha=z(j);
        alpha(1,1:nx)=pp(j,1:nx);
        SolbestX=Sol;
    elseif z(j)>z_alpha && z(j)<=z_beta
        z_beta=z(j);
        beta(1,1:nx)=pp(j,1:nx);
    elseif z(j)>z_beta && z(j)<=z_delta
        z_delta=z(j);
        delta(1,1:nx)=pp(j,1:nx);
    end
end
z_final(it)=z_alpha;
pos_final.BestCost = z_final(it);
% The Main Loop
for it=2:maxit
   x = 3.9 * x * (1 - x);  % Logistic Eq.(27)
    aa_max=0.9;
    aa_min=0.4;
    aa=2-2*(x*aa_min+(aa_max-aa_min)*it/maxit);%Eq.(28)

    k=k_max-(k_max-k_min)*(it-1)/(maxit-1);%Eq.(23) 
    Flag4ub=pp(j,:)>ub;
    Flag4lb=pp(j,:)<lb;
    pp(j,:)=(pp(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    for j=1:np        
        a_beta(1,1:nx)=(2*rand(1,nx)-ones(1,nx))*(aa); 
        c_beta(1,1:nx)=1;
        gv_beta(j,1:nx)=k*(sign(a_beta(1,1:nx)).*abs(gv_beta(j,1:nx)))+...
            a_beta(1,1:nx).*abs(c_beta(1,1:nx).*beta(1,1:nx)-pp(j,1:nx)); % Eq.(21)
        
        a_delta(1,1:nx)=((1/3)*rand(1,nx)-ones(1,nx))*(aa); 
        c_delta(1,1:nx)=1;
        gv_delta(j,1:nx)=k*(sign(a_delta(1,1:nx)).*abs(gv_delta(j,1:nx)))+...
            a_delta(1,1:nx).*abs(c_delta(1,1:nx).*delta(1,1:nx)-pp(j,1:nx)); % Eq.(22)
       
        a_alpha(1,1:nx)=(2*rand(1,nx)-ones(1,nx))*(aa); 
        c_alpha(1,1:nx)=1;
        gv_alpha(j,1:nx)=k*(sign(a_alpha(1,1:nx)).*abs(gv_alpha(j,1:nx)))+...
            a_alpha(1,1:nx).*abs(c_alpha(1,1:nx).*alpha(1,1:nx)-pp(j,1:nx)); % Eq.(20)
       C2= (1/3)*exp(-5*(it/maxit));%Eq.(31)
       C3= (1/3) * (1 - (it / maxit).^2);%Eq.(30)
       C1= 1 - (C2+ C3);%Eq.(29)
        pv(j,1:nx)=C1*((alpha-gv_alpha(j,1:nx))-pp(j,1:nx))+C2*((beta-gv_beta(j,1:nx))-pp(j,1:nx))+C3*((delta-gv_delta(j,1:nx))-pp(j,1:nx));
        % Return back the velocity of the particles if going beyond the velocity boundaries
        flag4lbv=pv(j,:)<velmin(1,:);
        flag4ubv=pv(j,:)>velmax(1,:);
        pv(j,:)=(pv(j,:)).*(~(flag4lbv+flag4ubv))+velmin.*flag4lbv+velmax.*flag4ubv;
        pp(j,1:nx)=pp(j,1:nx)+pv(j,1:nx); % Eq.(35)
            
        % Return back the position of the particles if going beyond the position boundaries
        flag4lbp=pp(j,:)<varmin(1,:);
        flag4ubp=pp(j,:)>varmax(1,:);
        pp(j,:)=(pp(j,:)).*(~(flag4lbp+flag4ubp))+varmin.*flag4lbp+varmax.*flag4ubp; 
    
    % Objective function evaluations and determining of the personal best solutions and objectives
       [z(j),Sol(j)] =fobj(pp(j,:));
    end
    if elitism==1
        for j=1:np
            if z_pbest(j)<z(j)
                z(j)=z_pbest(j);
                pp(j,:)=pp_pbest(j,:);
            else
                z_pbest(j)=z(j);
                pp_pbest(j,:)=pp(j,:);
            end
        end
    end
    for j=1:np
        if z(j)<=z_alpha
            z_alpha=z(j);
            alpha(1,1:nx)=pp(j,1:nx);
            SolbestX=Sol(j);
        elseif z(j)>z_alpha && z(j)<=z_beta
            z_beta=z(j);
            beta(1,1:nx)=pp(j,1:nx);
        elseif z(j)>z_beta && z(j)<=z_delta
            z_delta=z(j);
            delta(1,1:nx)=pp(j,1:nx);
        end
    end
    z_optimal(it)=z_alpha;
    optimal_pos(it,:)=alpha(1,:);

    z_iter(it)=z_optimal(it);
    z_final(it)=z_alpha;

logic=pp>ones(np,1)*ub;
    pp=logic.*(ones(np,1)*ub)+(1-logic).*pp;
    logic=pp<ones(np,1)*lb;
    pp=logic.*(ones(np,1)*lb)+(1-logic).*pp;
pos_final.Sol=SolbestX;
disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(z_final(it))]);
pos_final.BestCost = z_final(2:it);
pos_final.MaxIt = maxit;
end