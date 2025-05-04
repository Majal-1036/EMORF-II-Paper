function xp =robust_vbks_ind_self_modular_ind_sens_nsensors_xpp(y,x,x_0,dt,Q,R)
e_0=0.9;
    f_0=0.1;

    
    xe=x_0;

    m = x_0;
    P = Q;
    
    n = size(xe,1);
    alpha = 1;
    beta = 2;
    kappa = 0;
    lambda = alpha^2 * (n + kappa) - n;
    
    [WM,WC]=weights(n,alpha,beta,lambda);
    
    Mm = zeros(size(m,1), size(y,2) );
    Pm = zeros(size(P,1),size(P,2),size(y,2) );
    
    Mp = zeros(size(m,1), size(y,2) );
    Pp = zeros(size(P,1),size(P,2),size(y,2) );
    
    Ms = zeros(size(m,1), size(y,2) );
    Ps = zeros(size(P,1),size(P,2),size(y,2) );
    
    Dkp=zeros(size(m,1),size(m,1),size(y,2)-1);
%     N=length(y)
    N=size(y,2);
    
    Sig_zt_v=zeros(size(y,1),size(y,1),N);

    for jjj=1:N
    Sig_zt_v(:,:,jjj)=R;
    end
    
%    Forward Pass
    tau=1;

%     xpp=[];
    dim=size(y,1);
%     Rinv=inv(R);
    j=0;
    
     z_t_v=ones(dim,N);
     e_t_v = e_0*ones(dim,N);
     f_t_v = f_0*ones(dim,N);

while(tau>10^-4)
j=j+1;  
m = x_0;
P = Q;

for k=1:N

        SXp=sigma_gen(m,P,n,lambda);        
        % Propagate through the dynamic model 
        HX=dynamic_model(SXp,dt);        
        
        % Compute the predicted mean and covariance
        [mpr,Ppr]=mean_cov(WM,WC,HX);   
        Ppr = Ppr + Q;
        
        Mm(:,k) = mpr;
        Pm(:,:,k) = Ppr;
             
        SX=sigma_gen(mpr,Ppr,n,lambda);        
        
        if(k>1)
            Dkp(:,:,k)=cross_cov(SXp,HX,WC,m,mpr);
        end
        
        HY=meas_model(SX,dim);   
        
        [mu,Uk]=mean_cov(WM,WC,HY);    
        C=cross_cov(SX,HY,WC,mpr,mu);   
       
        Rkdi=inver_sig_test(Sig_zt_v(:,:,k));

        invS=Rkdi-Rkdi*((eye(dim)+Uk*Rkdi)^-1)*Uk*Rkdi;
        K=C*invS;
              
        m = mpr + K*(y(:,k) - mu);
        P = Ppr - C*K';        
    
        Mp(:,k) = m;
        Pp(:,:,k) = P;
end  

%     Backward Pass
Ms(:,N)=m;
Ps(:,:,N)=P;

for k=N-1:-1:1
G=Dkp(:,:,k+1)*inv(Pm(:,:,k+1));    
Ms(:,k)=Mp(:,k)+G*(Ms(:,k+1)-Mm(:,k+1));
Ps(:,:,k)=Pp(:,:,k)+G*(Ps(:,:,k+1)-Pm(:,:,k+1))*G';
end

xp(:,:,j)=Ms(:,:);

%Update    
for k=1:N

    SX=sigma_gen(Ms(:,k),Ps(:,:,k),n,lambda);
    HY=meas_model(SX,dim);        
    [my,Yk]=mean_cov(WM,WC,HY);  
            
    Bk=(my-y(:,k))*(my-y(:,k))'+Yk;
    Sig_zt=Sig_zt_v(:,:,k);
    z_t=z_t_v(:,k);
    e_t=e_t_v(:,k);
    f_t=f_t_v(:,k);

    
    for iz=1:dim
       
    z_t(iz)=1;
    pz1=exp( -.5*Bk(iz,iz)/R(iz,iz)+( psi(e_t(iz))-psi(e_t(iz)+f_t(iz))));
    z_t(iz)=0;
    pz0=exp( ( psi(f_t(iz))-psi(e_t(iz)+f_t(iz))) );
    rat=pz1/(pz1+pz0);
    
        if(rat>.5)
        z_t(iz)=1;
        else
        z_t(iz)=0;
        end
        
        if((z_t(iz))==0)
        Sig_zt(iz,iz)=inf;
        else
        Sig_zt(iz,iz)=R(iz,iz);
        end
    end 
    
    e_t = e_0 + z_t;
    f_t = f_0 + 1 - z_t;
    
    
    z_t_v(:,k)=z_t;
    e_t_v(:,k)=e_t;
    f_t_v(:,k)=f_t;
    Sig_zt_v(:,:,k)=Sig_zt;
    
    %%%%%%%%%%%%%%%%%%
    
    
    
    
end

if(j==1)
tau=1;
else
tau=norm(xp(:,:,j)-xp(:,:,j-1))/norm(xp(:,:,j-1));
end

end



time_mse_sp=mean(mean((x(:,:)-Ms(:,:)).^2));



function SS=sigma_gen(x,P,n,lambda)  
    P=nearestSPD(P);
    A = chol(P,'lower');
    SS = [zeros(size(x)) A -A];
    SS = sqrt(n + lambda)*SS + repmat(x,1,size(SS,2));
end

function HX=dynamic_model(SX,dt)
     HX=zeros(size(SX,1),size(SX,2));
        for i=1:size(SX,2)
        w_t=SX(5,i);
        F_t=[1 sin(w_t*dt)/w_t 0 (cos(w_t*dt)-1)/w_t 0;...
        0 cos(w_t*dt) 0 -sin(w_t*dt) 0;...
        0 (1-cos(w_t*dt))/w_t 1 sin(w_t*dt)/w_t 0;...
        0 sin(w_t*dt) 0 cos(w_t*dt) 0;...
        0 0 0 0 1];
        HX(:,i) = F_t*SX(:,i);
        end
end

function [x,P]=mean_cov(WM,WC,HX)
%         x= zeros(size(xin));
%         P = zeros(size(Pin));
x=0;P=0;
        for i=1:size(HX,2)
            x = x + WM(i) * HX(:,i);
        end
        for i=1:size(HX,2)
            P = P + WC(i) * (HX(:,i) - x) * (HX(:,i) - x)';
        end
end

function [C]=cross_cov(SX,HY,WC,xm,mu)
    C  = zeros(size(SX,1),size(HY,1));
    for i=1:size(SX,2)
            C = C + WC(i) * (SX(:,i) - xm) * (HY(:,i) - mu)';
    end
end

% function HY=meas_model(SX,m)
%            for ii=1:size(SX,2)
%            for jj=1:m
%               HY(jj,ii)=sqrt((SX(1,ii)-(jj-1)*350)^2+(SX(3,ii)-(350*mod(jj+1,2)))^2);
%            end
%            end
%  
% end

function HY=meas_model(SX,m)
           for jj=1:m
           for ii=1:size(SX,2)
              HY(jj,ii)=sqrt((SX(1,ii)-(1-1)*350)^2+(SX(3,ii)-(350*mod(1+1,2)))^2)-sqrt((SX(1,ii)-(jj+1-1)*350)^2+(SX(3,ii)-(350*mod(jj+1+1,2)))^2) ;
           end
           end
 
end


function [WM,WC]=weights(n,alpha,beta,lambda)
    WM = zeros(2*n+1,1);
    WC = zeros(2*n+1,1);
    
    for j=1:2*n+1
        if j==1
            wm = lambda / (n + lambda);
            wc = lambda / (n + lambda) + (1 - alpha^2 + beta);
        else
            wm = 1 / (2 * (n + lambda));
            wc = wm;
        end
        WM(j) = wm;
        WC(j) = wc;
    end
end



function [inv_sig]=inver_sig_test(Sig)
    inv_sig=zeros(size(Sig));
    indc_inf=find(diag(Sig)~=inf);
    inv_sig(indc_inf,indc_inf)=Sig(indc_inf,indc_inf)^-1;
end




end