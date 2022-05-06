% Stresses in a Cooling Elastic Shell
% Jakob Kintzele, Princeton Geosciences
% Last Update: 5/6/22

%%  ==== Parameters ==== %%

rho_i=920; %shell density
rho_w=1000; %ocean density
g=1.315; %surface gravity
E=5*10^9; %Young's Modulus, Manga & Wang (2009) 
nu=0.33; %Poisson Ratio
m=E/(2*(1+nu)); %elastic parameter
l=E*nu/((1+nu)*(1-2*nu));%elastic parameter 
beta=4*10^-10; %ocean compressibility, Manga & Wang (2009) 
R=1560*10^3; %satellite radius, Manga & Wang (2009) 
rc=1440*10^3; %core radius, Manga & Wang (2009) 
Ts=100;  % surface temperature, Nimmo (2002)
Te=180; % elastic layer basal temperature, Nimmo (2003) 
Tc=270; % conductive layer basal temperature 
Tb=270; % viscous layer basal temperature, Nimmo (2003) 
kappa=10^-6;% ice thermal diffusivity, Nimmo (2003) 
alpha=(10^-4)/3;% ice thermal expansivity, Nimmo (2003) 
c=2100; %specific heat capacity [J/(kg K)]
L=334*10^3; %Latent heat of fusion [J/kg]
lambda=0.636;
lambdae=0.27;

%% ==== Mesh ==== %%
N=100; % grid cells
Nt=12000; % timesteps
tconv=86400*365.25*10^6; % s->Ma
tmax=(25*10^3/(2*lambda))^2/kappa; % final time
t=linspace((2.4*10^3)^2/(4*kappa*lambda^2),tmax,Nt); % time vector
dt=(t(2)-t(1)); % timestep

rv=R-2*lambda.*sqrt(kappa.*(t)); % inner radius
rvp=-lambda.*sqrt(kappa./(t)); % inner radius velocity
re=R-2*lambdae.*sqrt(kappa.*(t)); %elastic radius
rep=-lambdae.*sqrt(kappa./(t)); % elastic radius velocity
dz=2*lambdae.*sqrt(kappa.*(t-t(1))); % freezing amt.
dzdt=lambdae.*sqrt(kappa./(t)); % freezing rate
rfix=zeros(4,N); % radial position 
for j=1:N
%rfix(1,j)=rv(end)+(j)*(R-rv(end))/(N);%element tops
rfix(2,j)=rv(end)+(j-0.5)*(R-rv(end))/(N);%element mids
%rfix(3,j)=rv(end)+(j-1)*(R-rv(end))/(N);%element bottoms
end
rfix(4,1)=rv(end); rfix(4,end)=R;
rfix(4,2)=rv(end)+10;
%rfix(4,3)=rv(end)+100;
rfix(4,3:end-2)=linspace(rv(end)+100,R-10,N-4);
%rfix(4,end-2)=R-10;
rfix(4,end-1)=R-1; %initial surface fracture depth

%% ==== Stefan Solution (Turcotte & Schubert, 2002) ==== %%
Tcart=zeros(N,Nt); % temperature 
for i=1:Nt
    for j=1:N
        if rfix(4,j)>=rv(i)
            Tcart(j,i)=(Te-Ts)*erf((R-rfix(4,j))/(2*sqrt(kappa*t(i))))/erf(lambdae)+Ts;
        elseif rfix(4,j)<rv(i) 
            Tcart(j,i)=Tb;    
        end
    end 
end

num=5;
f50=figure(1);
p50 = uipanel('Parent',f50,'BorderType','none'); 
p50.Title = 'Temperature Evolution'; 
p50.TitlePosition = 'centertop'; 
p50.FontSize = 12;
p50.FontWeight = 'bold';
subplot(1,2,1,'Parent',p50)
hold on 
grid on
plot(Tcart(:,1),(R-rfix(4,:))./10^3,'k-','LineWidth',2)
for i=1:num
plot(Tcart(:,i*length(t)/(num)),(R-rfix(4,:))./10^3,'LineWidth',2)
end
title('Stefan Solution')
ylim([0 25])
set(gca, 'YDir', 'Reverse')
ylabel('Depth [km]')
xlabel('Temperature [K]')
text(255, 1.5, '0.1 Ma')
text(255, 11.6, '2.5 Ma')
text(255, 19, '7.4 Ma')
text(220, 24.5, '12.2 Ma')

subplot(1,2,2,'Parent',p50)
hold on 
grid on
plot(t./tconv, (R-re)./10^3,'c-','LineWidth',2)
plot(t./tconv, (R-rv)./10^3,'b-','LineWidth',2)
legend('Elastic Layer', 'Full Shell')
xlabel('time [Ma]')
title('Shell Structure')
ylabel('Thickness [km]')
ylim([0 25])
set(gca, 'YDir', 'Reverse')
xlim([t(1)./tconv, t(end)./tconv])
%% ==== Thermal Stress (Nimmo, 2004) ==== %%
sigmaT=zeros(N,Nt);
sigmaTb=zeros(1,Nt);
duTBdt=zeros(1,Nt);
for i=2:Nt
    for j=1:N

        if Tcart(j,i)<=Te
sigmaT(j,i)=sigmaT(j,i-1)-dt*... 
alpha*E/(1-nu)*((kappa*(Te - Ts)*((-rfix(4,j) + R)/exp((rfix(4,j) - R)^2 ...
/(4*kappa*t(i))) + (2*kappa*t(i)*(2*rfix(4,j)^3 + (R - 2*lambdae*sqrt(kappa*t(i)))^3)* ...
     (R^2 + 4*kappa*t(i) + 4*kappa*lambdae^2*t(i) - 4*lambdae*R*sqrt(kappa*t(i)) ...
     - exp(lambdae^2)*(R^2 + 4*kappa*t(i)) +  ...
      2*exp(lambdae^2)*sqrt(pi)*R*sqrt(kappa*t(i))*erf(lambdae)))/ ...
    (exp(lambdae^2)*rfix(4,j)^3*(R^3 - (R - 2*lambdae*sqrt(kappa*t(i)))^3)) -  ...
   (2*kappa*t(i)*((rfix(4,j)^2 + 4*kappa*t(i))/exp((rfix(4,j) - R)^2/(4*kappa*t(i))) -  ...
      (R^2 + 4*kappa*(1 + lambdae^2)*t(i) - 4*lambdae*R*sqrt(kappa*t(i)))/exp(lambdae^2) -  ...
      2*sqrt(pi)*R*(sqrt(kappa*t(i))*erf(lambdae) + sqrt(kappa)*sqrt(t(i))*erf((rfix(4,j) -...
      R)/(2*sqrt(kappa)*sqrt(t(i)))))))/rfix(4,j)^3))/ ...
 (2*sqrt(pi)*(kappa*t(i))^(3/2)*erf(lambdae)));

sigmaTb(i)=-dt*... %basal stress
alpha*E/(1-nu)*((kappa*(Te - Ts)*((-re(i-1) + R)/exp((re(i-1) - R)^2 ...
/(4*kappa*t(i))) + (2*kappa*t(i)*(2*re(i-1)^3 + (R - 2*lambdae*sqrt(kappa*t(i)))^3)* ...
     (R^2 + 4*kappa*t(i) + 4*kappa*lambdae^2*t(i) - 4*lambdae*R*sqrt(kappa*t(i)) ...
     - exp(lambdae^2)*(R^2 + 4*kappa*t(i)) +  ...
      2*exp(lambdae^2)*sqrt(pi)*R*sqrt(kappa*t(i))*erf(lambdae)))/ ...
    (exp(lambdae^2)*re(i-1)^3*(R^3 - (R - 2*lambdae*sqrt(kappa*t(i)))^3)) -  ...
   (2*kappa*t(i)*((re(i-1)^2 + 4*kappa*t(i))/exp((re(i-1) - R)^2/(4*kappa*t(i))) -  ...
      (R^2 + 4*kappa*(1 + lambdae^2)*t(i) - 4*lambdae*R*sqrt(kappa*t(i)))/exp(lambdae^2) -  ...
      2*sqrt(pi)*R*(sqrt(kappa*t(i))*erf(lambdae) + sqrt(kappa)*sqrt(t(i))*erf((re(i-1) -...
      R)/(2*sqrt(kappa)*sqrt(t(i)))))))/re(i-1)^3))/ ...
 (2*sqrt(pi)*(kappa*t(i))^(3/2)*erf(lambdae)));

duTBdt(i)=(1-nu)/E*re(i)*... % basal displacement rate 
    (alpha*E/(1-nu)*((kappa*(Te - Ts)*((-re(i) + R)/exp((re(i) - R)^2 ...
/(4*kappa*t(i))) + (2*kappa*t(i)*(2*re(i)^3 + (R - 2*lambdae*sqrt(kappa*t(i)))^3)* ...
     (R^2 + 4*kappa*t(i) + 4*kappa*lambdae^2*t(i) - 4*lambdae*R*sqrt(kappa*t(i)) ...
     - exp(lambdae^2)*(R^2 + 4*kappa*t(i)) +  ...
      2*exp(lambdae^2)*sqrt(pi)*R*sqrt(kappa*t(i))*erf(lambdae)))/ ...
    (exp(lambdae^2)*re(i)^3*(R^3 - (R - 2*lambdae*sqrt(kappa*t(i)))^3)) -  ...
   (2*kappa*t(i)*((re(i)^2 + 4*kappa*t(i))/exp((re(i) - R)^2/(4*kappa*t(i))) -  ...
      (R^2 + 4*kappa*(1 + lambdae^2)*t(i) - 4*lambdae*R*sqrt(kappa*t(i)))/exp(lambdae^2) -  ...
      2*sqrt(pi)*R*(sqrt(kappa*t(i))*erf(lambdae) + sqrt(kappa)*sqrt(t(i))*erf((re(i) -...
      R)/(2*sqrt(kappa)*sqrt(t(i)))))))/re(i)^3))/ ...
 (2*sqrt(pi)*(kappa*t(i))^(3/2)*erf(lambdae)))) ...
-alpha*(Te-Ts)/erf(lambdae)*exp(-lambdae^2)*lambdae/t(i)*re(i); %thermal contraction term (Hillier & Squyres)
        end
    end
end

num=5;
figure(2)
grid on 
hold on
for i=1:num
plot(sigmaT(:,i*length(t)/(num))./10^6,(R-rfix(4,:))./10^3);
%t(i*length(t)/(num))./tconv %time slice values [Ma]
end
plot(sigmaT(:,3)./10^6,(R-rfix(4,:))./10^3,'k-','LineWidth',2)
plot(sigmaT(:,end)./10^6,(R-rfix(4,:))./10^3,'r-','LineWidth',2)
text(0.5, 0.3, '0.1 Ma')
text(0.5, 5.1, '2.5 Ma')
text(0.5, 8.05, '7.4 Ma')
text(0.5, 10.4, '12.2 Ma')
text(-6.5, 1, 'tension')
text(15, 1, 'compression')
set(gca, 'YDir', 'Reverse')
ylabel('Depth [km]')
xlabel('Tangential Stress [MPa]')
ylim([0 11])
xlim([-7 24])
title('Thermal Stress')

figure(22)
grid on 
hold on
for i=1:num
plot((1-nu)/E*sigmaT(:,i*length(t)/(num))*100,(R-rfix(4,:))./10^3);
%t(i*length(t)/(num))./tconv %time slice values [Ma]
end
plot((1-nu)/E*sigmaT(:,3)*100,(R-rfix(4,:))./10^3,'k-','LineWidth',2)
plot((1-nu)/E*sigmaT(:,end)*100,(R-rfix(4,:))./10^3,'r-','LineWidth',2)
title('Tangential Strain [%]') 
set(gca, 'YDir', 'Reverse')
%%  ==== Volume Change Effect (Manga & Wang, 2007) ==== %%
sigmaV=zeros(N,Nt);
dPexdt=zeros(1,Nt);
Pex=zeros(1,Nt);
A=zeros(1,Nt);
dAdt=zeros(1,Nt);
uP=zeros(1,Nt);
uP(1)=dz(1)*(rho_w-rho_i)/rho_w;
for i=2:Nt

A(i)=3*rv(i)^2/(beta*(rv(i)^3-rc^3));
dAdt(i)=3/beta*rvp(i)*(2*rv(i)*(rv(i)^3-rc^3)-rv(i)^2*3*rv(i)^2)/((rv(i)^3-rc^3)^2); 

    dPexdt(i)=(-dAdt(i)* ... %pressurization rate
        ( 0 ... % thermal effect (set to zero for purely volume change effects)
        -dz(i)*(rho_w-rho_i)/rho_w)-A(i)* ...
        ( 0 ...% thermal effect (set to zero for purely volume change effects)
         -dzdt(i)*(rho_w-rho_i)/rho_w)-dAdt(i)*uP(i-1))/ ...
        (1+re(i)/E*((1+0.5*(R/re(i))^3)/((R/re(i))^3-1)*(1-nu)+nu)*(dAdt(i)*dt+A(i)));
    
    uP(i)=uP(i-1)+dt*dPexdt(i)*re(i)/E*((1+0.5*(R/re(i))^3)/((R/re(i))^3-1)*(1-nu)+nu); %Basal disp.
    
    Pex(i)=Pex(i-1)+dt*dPexdt(i); %pressure
    for j=1:N 
        if Tcart(j,i)<=Te
sigmaV(j,i)=sigmaV(j,i-1)-dt*(1+1/2*(R/rfix(4,j))^3)/((R/re(i))^3-1)*dPexdt(i);
        elseif Tcart(j,i)>Te
sigmaV(j,i)=sigmaV(j,i-1)+dt*dPexdt(i);
        end
    end
end


num=5;
f7=figure(7);
p = uipanel('Parent',f7,'BorderType','none'); 
p.Title = 'Volume Change Effect'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';
subplot(1,2,1,'Parent',p)
grid on 
hold on
plot(sigmaV(:,3)./10^6,(R-rfix(4,:))./10^3,'k-','LineWidth',2);
for i=1:num
plot(sigmaV(:,i*length(t)/(num))./10^6,(R-rfix(4,:))./10^3);
%t(i*length(t)/(num))./tconv %timeslice values [Ma]
end
plot(sigmaV(:,end)./10^6,(R-rfix(4,:))./10^3,'r-','LineWidth',2);
text(-0.8, 0.5, '0.1 Ma')
text(-1.8, 0.5, '2.5 Ma')
text(-3.1, 0.5, '7.4 Ma')
text(-4, 3, '12.2 Ma')
set(gca, 'YDir', 'Reverse')
xlabel('Tangential Stress [MPa]')
ylabel('Depth [km]')
ylim([0 (R-re(end))/10^3])
%xlim([-10 0.5])
title('Elastic Layer')
subplot(2,2,4,'Parent',p)
grid on 
hold on
plot(t(1:100:end)./tconv, (-dz(1:100:end)*(rho_w-rho_i)/rho_w+uP(1:100:end)),'c-','LineWidth',3);
ylabel('Radial Displacement [m]')
title('Net Basal Displacement')
xlabel('time [Ma]')
xlim([t(1)./tconv, t(end)./tconv])
subplot(2,2,2,'Parent',p)
grid on 
hold on
plot(t(1:100:end)./tconv,Pex(1:100:end)./10^3,'b-','LineWidth',3)
ylabel('Excess Pressure [kPa]')
xlim([t(1)./tconv, t(end)./tconv])
xlabel('time [Ma]')
title('Viscous Ice and Ocean')


%% ==== Overpressure from Thermal Contraction ==== %%
% == Basal Displacement from Thermal Contraction ==%
uTB=zeros(1,Nt);
duTBdt=zeros(1,Nt);
for i=1:Nt
uTB(i)=re(i)*(1/2)*alpha*(Ts-Te);
if i>1
duTBdt(i)=(uTB(i)-uTB(i-1))/dt;
end
end

% == overpressurization and response == %

uPT=zeros(1,Nt);
dPexTdt=zeros(1,Nt);
PexT=zeros(1,Nt);
%uPT(1)=-uTB(1);
for i=2:Nt
    dPexTdt(i)=(-dAdt(i)*(uTB(i))-A(i)*(duTBdt(i))-dAdt(i)*uPT(i-1))/ ...
        (1+re(i)/E*((1+0.5*(R/re(i))^3)/((R/re(i))^3-1)*(1-nu)+nu)*(dAdt(i)*dt+A(i)));
    
    uPT(i)=uPT(i-1)+dt*dPexTdt(i)*re(i)/E*((1+0.5*(R/re(i))^3)/((R/re(i))^3-1)*(1-nu)+nu); %Hillier & Squyres
    PexT(i)=PexT(i-1)+dt*dPexTdt(i);
end

% == induced stress field ==%
sigmaVT=zeros(N,Nt);
for i=2:Nt
    for j=1:N 
        if Tcart(j,i)<=Te
sigmaVT(j,i)=sigmaVT(j,i-1)-dt*(1+1/2*(R/rfix(4,j))^3)/((R/re(i))^3-1)*dPexTdt(i);
        elseif Tcart(j,i)>Te
sigmaVT(j,i)=sigmaVT(j,i-1)+dt*dPexTdt(i);
        end
    end
end

num=5;
f9=figure(99);
p9 = uipanel('Parent',f9,'BorderType','none'); 
p9.Title = 'Thermal Contraction'; 
p9.TitlePosition = 'centertop'; 
p9.FontSize = 12;
p9.FontWeight = 'bold';
subplot(1,2,1,'Parent',p9)
grid on 
hold on
for i=1:num
plot(sigmaVT(:,i*length(t)/(num))./10^6,(R-rfix(4,:))./10^3);
%t(i*length(t)/(num))./tconv %time slice values [Ma]
end
plot(sigmaVT(:,3)./10^6,(R-rfix(4,:))./10^3,'k-','LineWidth',2);
plot(sigmaVT(:,end)./10^6,(R-rfix(4,:))./10^3,'r-','LineWidth',2);

set(gca, 'YDir', 'Reverse')
xlabel('Tangential Stress [MPa]')
ylabel('Depth [km]')
ylim([0 (R-re(end))/10^3])
title('Elastic Layer')
text(-.35, 0.5, '0.1 Ma')
text(-.8, 0.5, '2.5 Ma')
text(-1.4, 0.5, '7.4 Ma')
text(-1.9, 3, '12.2 Ma')

subplot(2,2,2,'Parent',p9)
title('Viscous Ice and Ocean') 
grid on 
hold on
plot(t(1:10:end)./tconv,(PexT(1:10:end))./10^3,'b-','LineWidth',3)
ylabel('Excess Pressure [kPa]')
xlabel('Time [Ma]')
xlim([t(1)./tconv, t(end)./tconv])

subplot(2,2,4,'Parent',p9)
grid on 
hold on
%plot(t./tconv,(uTB)./10^3,'r-','LineWidth',3)
%plot(t./tconv,(-dz)./10^3,'c-','LineWidth',3)
%plot(t./tconv,(uPT)./10^3,'b-','LineWidth',3)
%plot(t./tconv,(uP)./10^3,'b--','LineWidth',3)
plot(t(1:10:end)./tconv,(uPT(1:10:end)+uTB(1:10:end)),'c-','LineWidth',3)
title('Net Basal Displacement')
ylabel('Radial Displacement [m]')
xlabel('Time [Ma]')
xlim([t(1)./tconv, t(end)./tconv])

%% ==== Overburden ==== %%
sigmaW=zeros(N,Nt);
for i=1:Nt
    for j=1:N
        if Tcart(j,i)>=Tb
            sigmaW(j,i)=g*(rho_i*(R-rv(i))+rho_w*(rv(i)-rfix(4,j)));
        elseif Tcart(j,i)<Tb
            sigmaW(j,i)=rho_i*g*(R-rfix(4,j));
        end
    end
end

figure(7)
grid on 
hold on
for i=1:num
plot(sigmaW(:,i*length(t)/(num))./10^6,(R-rfix(4,:))./10^3);
%t(i*length(t)/(num))./tconv %time slice values [Ma]
end
plot(sigmaW(:,end)./10^6,(R-rfix(4,:))./10^3,'r-','LineWidth',2);
plot(sigmaW(:,1)./10^6,(R-rfix(4,:))./10^3,'k-','LineWidth',2);
set(gca, 'YDir', 'Reverse')
xlabel('Overburden [MPa]')
ylabel('Depth [km]')
ylim([0 (R-rv(end))/10^3])

%% ==== Net Stress ==== %%
sigma=sigmaV+sigmaW+sigmaVT+sigmaT;

f22=figure(10);
p2 = uipanel('Parent',f22,'BorderType','none'); 
p2.Title = 'Cumulative State'; 
p2.TitlePosition = 'centertop'; 
p2.FontSize = 12;
p2.FontWeight = 'bold';
subplot(1,2,1,'Parent',p2)
grid on 
hold on
num=5;
for i=1:num
plot(sigma(:,i*length(t)/(num))./10^6,(R-rfix(4,:))./10^3);
%t(i*length(t)/(num))./tconv %time slice values [Ma]
end
plot(sigma(:,3)./10^6,(R-rfix(4,:))./10^3,'k-','LineWidth',2);
plot(sigma(:,end)./10^6,(R-rfix(4,:))./10^3,'r-','LineWidth',2);
xlim([-10 33]);
set(gca, 'YDir', 'Reverse')
xlabel('Tangential Stress [MPa]')
%text(1, 0.5, '0.1 Ma')
text(4, 2.5, '0.1 Ma')
text(7,5, '2.5 Ma')
text(12, 8,  '7.4 Ma')
text(15, (R-re(end))./10^3, '12.2 Ma')

ylabel('Depth [km]')
ylim([0 25])
title('Total Ice Shell Stress')

subplot(1,2,2,'Parent',p2)
grid on 
hold on
%plot(t(1:100:end)./tconv,(Pex(1:100:end)+PexT(1:100:end))./10^3,'b-','LineWidth',3)
plot(t(1:100:end)./tconv,(PexT(1:100:end)+Pex(1:100:end))./10^3,'c-','LineWidth',3)
%plot(t(1:100:end)./tconv,(PexT(1:100:end))./10^3,'k--','LineWidth',3)
%legend('Cumulative Pressure','Volume Change','Thermal Contraction','Location','NW')
ylabel('Excess Pressure [kPa]')
xlabel('Time [Ma]')
xlim([t(1)./tconv, t(end)./tconv])
title('Ocean Pressurization')
