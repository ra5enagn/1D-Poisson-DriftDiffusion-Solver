%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                         %%
%%             Equilibrium Poisson Equation Solver for Diodes              %%
%%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc ;
clear;
close all;

% Defining the Fundamental and Material Constants %

q   = 1.602E-19;
kb  = 1.38E-23;
eps = 1.05E-12; %F/cm
T   = 300;
ni  = 1.5E10;
Vt  = kb*T/q;
RNc = 2.8E19;
dEc = Vt*log(RNc/ni);

u_n = 1650;   % Electron mobility in cm^2/(V.s)
u_p = 650;    % Hole mobility in cm^2/(V.s)
Dn = Vt * u_n;  % Diffusion constant for electrons
Dp = Vt * u_p;  % Diffusion constant for hole

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define input parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read Doping Values %

Na = 1E17;
Nd = 1E17;

% Read Device dimensions: extent of the n and p-regions

x_n = 0.5e-4;
x_p = 0.5e-4;

delta_acc = 1E-9;       % Preset the Tolerance for convergence

%biasing voltage
v_bias = 0.75;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate relevant parameters for the simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ldn = sqrt(eps*Vt/(q*Nd));
Ldp = sqrt(eps*Vt/(q*Na));
Ldi = sqrt(eps*Vt/(q*ni));

x_max = x_n + x_p;

% Setting the grid size based on the extrinsic Debye lengths %
% minimum Debye length wins

dx = Ldn;
if(dx > Ldp)
    dx=Ldp;
end
dx = dx/5;

% Calculate the required number of grid points and renormalize dx %

n_max = x_max/dx;
n_max = round(n_max);

dx = dx/Ldi;  % Renormalize lengths with Ldi

% set up size of variables used

dop = zeros(n_max,1);
fi  = zeros(n_max,1);
fip = zeros(n_max,1);
n   = zeros(n_max,1);
p   = zeros(n_max,1);
a   = zeros(n_max,1);
b   = zeros(n_max,1);
c   = zeros(n_max,1);
f   = zeros(n_max,1);
v   = zeros(n_max,1);
Ec   = zeros(n_max,1);
el_field   = zeros(n_max,1);
ro   = zeros(n_max,1);
xx1 = zeros(n_max,1);
delta_vec = zeros(n_max,1);


an   = zeros(n_max,1);
bn   = zeros(n_max,1);
cn   = zeros(n_max,1);
fn   = zeros(n_max,1);

ap   = zeros(n_max,1);
bp   = zeros(n_max,1);
cp   = zeros(n_max,1);
fp   = zeros(n_max,1);

jp   = zeros(n_max-1,1);
jn   = zeros(n_max-1,1);
jt   = zeros(n_max-1,1);

jp_zero   = zeros(n_max-1,1);
jn_zero   = zeros(n_max-1,1);
jt_zero   = zeros(n_max-1,1);


% Set up the doping vector dop(x)=Nd(x)-Na(x) that is normalized with ni

for i = 1:n_max
    xxx = i*dx*Ldi;
    if(xxx <= x_n)
        dop(i) = Nd/ni;
    elseif ((xxx > x_n))
        dop(i) = -Na/ni;
    end
end

% Initialize the potential based on the requirement of charge
% neutrality throughout the whole structure

for i = 1: n_max
    zz = 0.5*dop(i);
    if(zz > 0)
        xx = zz*(1 + sqrt(1+1/(zz*zz)));
    elseif(zz < 0)
        xx = zz*(1 - sqrt(1+1/(zz*zz)));
    end
    fi(i) = log(xx);    % initial potential
    n(i) = xx;          % initial electron density
    p(i) = 1/xx;        % initial hole density
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %%
%        Solving for the Equillibirium Case         %%
%                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (A) Define the elements of the coefficient matrix for the internal nodes and
%     initialize the forcing function

dx2 = dx*dx;
for i = 2: n_max-1
    a(i) = 1/dx2;
    c(i) = 1/dx2;
    b(i) = -(2/dx2+exp(fi(i))+exp(-fi(i)));
    f(i) = exp(fi(i))-exp(-fi(i))-dop(i) - fi(i)*(exp(fi(i))+exp(-fi(i)));
end

% (B) Define the elements of the coefficient matrix and initialize the forcing
%     function at the ohmic contacts

a(1) = 0;
c(1) = 0;
b(1) = 1;
f(1) = fi(1);
a(n_max) = 0;
c(n_max) = 0;
b(n_max) = 1;
f(n_max) = fi(n_max);

% (C) Start the iterative procedure for the solution of the linearized Poisson
%     equation using LU decomposition method:

flag_conv = 0;  % set convergence flag of the Poisson loop to zero
k_iter= 0;

while(~flag_conv && k_iter <= 2)

    k_iter = k_iter + 1;

    % call LU Decomposition function to solve for fip

    [fip] = LU_Decomposition(a,b,c,f,n_max);

    for i = 1:n_max
        delta_vec(i) = fip(i) - fi(i);
    end

    delta_max = max(abs(delta_vec));

    % Test convergence. If convergence not achieved:
    %   (1) update central coefficient b(i)
    %   (2) update forcing function f(i)

    if(delta_max < delta_acc)
        flag_conv = 1;
    else
        for i = 2: n_max-1
            b(i) = -(2/dx2+exp(fip(i))+exp(-fip(i)));
            f(i) = exp(fip(i))-exp(-fip(i))-dop(i)-fip(i)*(exp(fip(i))+exp(-fip(i)));
        end
    end
    fi = fip;
end
if(k_iter >1000)
    disp('Failed operation')
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Non Equilibrium solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%new electrons and holes from poison solver

for i=1:n_max
    Ec(i) = dEc - Vt*fi(i);
    n(i) = exp(fi(i)); %n to be normalised...
    p(i) = exp(-fi(i));
    ro(i) = p(i) - n(i) + ni*dop(i);
end
    for i=1:n_max-1
        jn_zero(i) = (((q*n(i+1)*ni*Dn)/(dx * Ldi))* berfun(fi(i+1)-fi(i))) - (((q*n(i)*ni*Dn)/(dx * Ldi))* berfun(fi(i)-fi(i+1)));
        jp_zero(i) = (((q*p(i)*ni*Dn)/(dx * Ldi))* berfun(fi(i+1)-fi(i))) - (((q*p(i+1)*ni*Dn)/(dx * Ldi))* berfun(fi(i)-fi(i+1)));
        jt_zero(i) = jn(i)+jp(i);
    end
fi_zero = fi;

Ec_zero = Ec;
n_zero = n;
p_zero = p;


function Ber = berfun(x)
flag_conv = 0;
if(x>0.01)
    Ber = x*exp(-x)/(1-exp(-x));
elseif(x<0 && abs(x)> 0.01)
    Ber = x/(exp(x)-1);
elseif(x == 0)
    Ber = 1;
else
    temp_term = 1;
    sum = temp_term;
    i = 0;
    while(~flag_conv)
        i = i + 1;
        temp_term = temp_term*x/(i+1);
        if( sum + temp_term == sum)
            flag_conv = 1;
        end
        sum = sum + temp_term;
    end
    Ber = 1/sum;
end
end


%current densities from new n &p

an(1) = 0;
bn(1) = 1;
cn(1) = 0;
fn(1) = n(1);

an(n_max) = 0;
bn(n_max) = 1;
cn(n_max) = 0;
fn(n_max) = n(n_max);

ap(1) = 0;
bp(1) = 1;
cp(1) = 0;
fp(1) = p(1);

ap(n_max) = 0;
bp(n_max) = 1;
cp(n_max) = 0;
fp(n_max) = p(n_max);


f(1)=fi(1);
iteration=0;

for counter = 0:Vt:abs(v_bias)
    if v_bias < 0 
        volt_bias = -1;
    else
        volt_bias = +1;
    end
    fi(n_max) =fi(n_max) + volt_bias; %(vt/vt  = cmpounding step size)
    f(n_max) = fi(n_max);

    TAUN0 = 10e-6;
    TAUP0 = 10e-6;
    iteration =iteration++1;

    k_iter= 0;

    while(1)

        k_iter = k_iter + 1;


        for i = 2:n_max-1
            tp = TAUP0;
            tn = TAUN0;

            an(i) = (Dn) * berfun(fi(i-1) - fi(i));
            bn(i) = -((Dn) * berfun(fi(i) - fi(i+1)) + (Dn) * berfun(fi(i) - fi(i-1)));
            cn(i) = (Dn) * berfun(fi(i+1) - fi(i));
            fn(i) = ((n(i)*p(i) - 1)/(tp * (n(i)+1) + tn * (p(i)+1)))*(dx2* Ldi^2);

            ap(i) = (Dp/(dx2* Ldi^2)) * berfun(fi(i) - fi(i-1));
            bp(i) = -((Dp/(dx2* Ldi^2)) * berfun(fi(i+1) - fi(i)) + (Dp/(dx2* Ldi^2)) * berfun(fi(i-1) - fi(i)));
            cp(i) = (Dp/(dx2* Ldi^2)) * berfun(fi(i) - fi(i+1));
            fp(i) = ((n(i)*p(i)) - 1)/(tp*(n(i)+1)+tn*(p(i)+1));
        end

        [n] = LU_Decomposition(an,bn,cn,fn,n_max);
        [p] = LU_Decomposition(ap,bp,cp,fp,n_max);

        for i = 2:n_max-1
            b(i) = -((2/dx2)+n(i)+p(i));
            f(i) = -(p(i)-n(i) + dop(i))-(fi(i)*(n(i)+p(i)));
        end


        % call LU Decomposition function to solve for fip

        [fip] = LU_Decomposition(a,b,c,f,n_max);

        for i = 1:n_max
            delta_vec(i) = fip(i) - fi(i);
        end

        delta_max = max(abs(delta_vec));

        % Test convergence. If convergence not achieved:
        %   (1) update central coefficient b(i)
        %   (2) update forcing function f(i)

        if(delta_max < delta_acc)
            break;
        else
            for i = 2: n_max-1
                b(i) = -(2/dx2+exp(fip(i))+exp(-fip(i)));
                f(i) = exp(fip(i))-exp(-fip(i))-dop(i)-fip(i)*(exp(fip(i))+exp(-fip(i)));
            end
        end
        fi = fip; 
    end
    %fprintf('iteration %d\n',k_iter)
    for i=1:n_max-1
        jn(i) = ((q*Dn)/(dx*Ldi))*(n(i+1)*berfun(fi(i+1)-fi(i))-n(i)*berfun(fi(i)-fi(i+1)))*ni;
        jp(i) = ((q*Dp)/(dx*Ldi))*(p(i)*berfun(fi(i+1)-fi(i))-p(i+1)*berfun(fi(i)-fi(i+1)))*ni;
        jt(i) = jn(i)+jp(i);
    end
    jt_avg(iteration) = -1*mean(jt);
    voltage (iteration)= counter;
end
for i=1:n_max
    Ec(i) = dEc - Vt*fi(i);
    ro(i) = p(i) - n(i) + ni*dop(i);
end

xx1(1) = dx*Ldi*1e4; % convert distance in um
for i = 2:n_max-1
    xx1(i) = xx1(i-1) + dx * Ldi * 1e4;
    el_field(i) = -(fi(i+1) - fi(i-1)) * Vt /(2 * dx * Ldi);
    el_field_zero(i) = -(fi_zero(i+1) - fi_zero(i-1)) * Vt /(2 * dx * Ldi);
end
xx1(n_max) = xx1(n_max-1) + dx*Ldi*1e4;
el_field(1) = el_field(2);
el_field(n_max) = el_field(n_max-1);

el_field_zero(1) = el_field_zero(2);
el_field_zero(n_max) = el_field_zero(n_max-1);


figure(1)
plot(xx1, Ec.','DisplayName',['E_{c} @ V_{bias} = ' num2str(v_bias) 'V'],'LineWidth',2,'LineStyle','-')
hold on;
 grid on;
 plot(xx1, Ec.' - 1.12,'DisplayName',['E_{v}  @ V_{bias} = ' num2str(v_bias) 'V'],'LineWidth',2,'LineStyle','-');
 plot(xx1,0,'HandleVisibility','off'); 
 plot(xx1, Ec_zero.','DisplayName',['E_{c}  @ V_{bias} = 0V'],'LineWidth',2,'LineStyle','--')
 plot(xx1, Ec_zero.'-1.12,'DisplayName','E_{v}  @ V_{bias} = 0V','LineWidth',2,'LineStyle','--')
xlabel('x [um]');
ylabel('Conduction and Valence Bands [eV]');
title('Conduction and Valence bands vs Position');
hold off;
legend('NumColumns',2)

figure(2)
semilogy(xx1, n*ni,'DisplayName',['n  @ V_{bias} = ' num2str(v_bias) 'V'],'LineWidth',2,'LineStyle','-')
hold on;
grid on
semilogy(xx1, p*ni,'DisplayName',['P  @ V_{bias} = ' num2str(v_bias) 'V'],'LineWidth',2,'LineStyle','-')
semilogy(xx1, n_zero*ni,'DisplayName','n  @ V_{bias} = 0V','LineWidth',2,'LineStyle','--')
semilogy(xx1, p_zero*ni,'DisplayName','p  @ V_{bias} = 0V','LineWidth',2,'LineStyle','--')
xlabel('x [um]');
ylabel('Electron & Hole Densities [1/cm^3]');
title('Electron & Hole Densities vs Position');
legend('Location','east');


figure(3)
plot(voltage,jt_avg ,'DisplayName',['Hole current @ V_{bias} = ' num2str(v_bias) 'V']','LineWidth',2,'LineStyle','-')
grid on;
xlabel('Voltage');
ylabel('Current A/um');
title('Current Vs Voltage');

figure(4)
plot(xx1(1:n_max-1), -jp,'DisplayName',['Hole current @ V_{bias} = ' num2str(v_bias) 'V'],'LineWidth',2,'LineStyle','-')
hold on;
grid on;
plot(xx1(1:n_max-1), -jn,'DisplayName',['Electron current @ V_{bias} = ' num2str(v_bias) 'V']','LineWidth',2,'LineStyle','-')
plot(xx1(1:n_max-1), -jt,'DisplayName',['Total current @ V_{bias} = ' num2str(v_bias) 'V']','LineWidth',2,'LineStyle','-')
plot(xx1(1:n_max-1), -jp_zero,'DisplayName','Hole current @ V_{bias} = 0.V','LineWidth',2,'LineStyle','--')
plot(xx1(1:n_max-1), -jn_zero,'DisplayName','Electron current @ V_{bias} = 0.V','LineWidth',2,'LineStyle',':')
plot(xx1(1:n_max-1), -jt_zero,'DisplayName','Total current @ V_{bias} = 0.V','LineWidth',2,'LineStyle','-.')
xlabel('x [um]');
ylabel('Current density A/um^{2}');
title('Current density Vs Position');
hold off;
legend ('NumColumns',1,'Location','east')

max_field= max(el_field);
ind = find(el_field == max_field);
max_field_zero= max(el_field_zero);
ind_zero = find(el_field_zero == max_field_zero,1,"last");

% Calculate quasi-Fermi levels
for i = 1:n_max
    qfn(i) = (fi(i)-log(n(i)))*Vt;    % Electron quasi-Fermi level
    qfp(i) = (fi(i)+log(p(i)))*Vt;    % Hole quasi-Fermi level
end

figure(5)
plot(xx1,qfn,'DisplayName','Quasi-Fermi for electrons','LineWidth',2,'LineStyle','-')
hold on;
grid on;
plot(xx1,qfp,'DisplayName','Quasi-Fermi for holes','LineWidth',2,'LineStyle','-')
hold off
title ('Quasi-Fermi Energy level Vs position')
ylabel('Quasi Fermi energy level in eV')
xlabel('poition in um')
legend ('NumColumns',1,'Location','east')


figure(6)
plot(xx1, el_field,'DisplayName',['Electron Field @ V_{bias} = ' num2str(v_bias) 'V'],'LineWidth',2,'LineStyle','-')
hold on;
grid on;
plot(xx1, el_field_zero,'DisplayName','Electron Field (@ V_{bias} = 0V)','LineWidth',2,'LineStyle','--')
plot (xx1(ind),max_field,'Marker','*','DisplayName',['max Electron Field (@ V_{bias} = 0 V) = '  num2str(max_field_zero)]);
plot (xx1(ind_zero),max_field_zero,'Marker','o','DisplayName',['max Electron Field (@ V_{bias} = ' num2str(v_bias) 'V) = ' num2str(max_field)]);
hold off
xlabel('x [um]');
ylabel('Electric Field [V/cm]');
title('Field Profile vs Position');
legend Show