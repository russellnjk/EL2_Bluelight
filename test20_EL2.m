function result=test20_EL2()
close all;
clear all;
global default apply tspan1 tspan2 Init1 Init2 M1 M2 lb ub fval S1 S2;

dt = 1; %min
tspan1 = 0:dt:6; %min
tspan2 = [0 0.75:1:6.75];

default=[
1;           %hill coeff    
0.152;       %k_a 1/hr k_d 1/hr
4.32e-05;    %synmRFP M/hr
2.305e-1;    %km M
8.316;       %degmRFP 1/hr
0.312;       %synRFP 1/hr
0.1124;      %degRFP 1/hr Ctrl
0.2418/0.1124;      %degRFP 1/hr DAS
0.2888/0.2418;      %degRFP 1/hr AAV
0.4469/0.2888;];    %degRFP 1/hr YbaQ

ub=[1.2;2; 0.001;5;0.8;0.3;3.5;2;2;];
lb=[1;0.0001;0.000001;0.0001;0.1;0.1;1;1;1;];
apply=[1;1;1;1;0;1;1;1;1;1;];

Init1  = [0          0          0           0;  %EL2a
         1.0*10^-6  1.0*10^-6  0.75*10^-6  0.43*10^-6;  %RFP
         0          0          0           0;];  %mRFP

Init2  = [0          0           0           0;  %EL2a
         0.25*10^-6  0.25*10^-6  0.1*10^-6  0.1*10^-6;  %RFP
         0           0           0           0;];  %mRFP

  
file1="Inducible_2hrON4hrOFF.csv";
M1=dlmread(file1, ',', [2 0 8 4]);
disp(M1);
%Insert M2 definition here
file2="Inducible_3hrON3hrOFF.csv";
M2=dlmread(file2, ';', [2 0 8 4]);
disp(M2);

stdev1="";
if stdev1==""
    S1=[];
else
    S1=dlmread(stdev1, ',', [1 0 7 4]);
end
disp(S1);
stdev2="";
if stdev2==""
    S2=[];
else
    S2=dlmread(stdev2, ',', [1 0 8 4]);
end
disp(S2);


i=[default(1:6) default(1:6) default(1:6) default(1:6);default(7) default(8) default(9) default(10);];
initial=get_sse2(i);
disp(initial);

trials=3;
outputs=zeros(7,length(default)-6,trials);
fval=[];
for i=1:trials
    disp("Trial number: "+num2str(i));
    outputs(:,:,i)=get_min(); 
    disp(outputs(:,:,i));
end

[best_sse,index]=min(fval);
disp(index);
disp(best_sse);
result=outputs(:,:,index);

final  =get_sse2(result);
disp(final);
end

function result=get_min()
global apply default lb ub fval;
i=0;
for ii=1:length(apply)
    if apply(ii)==1
        i=i+1;
    else
    end
end
x0=ones(i,length(default(1,:)));
i=1;
for ii=1:length(default(:,1))
    if apply(ii)==1
        x0(i,:)=default(ii,:);
        i=i+1;
    else
    end
end
disp('x0');
disp(x0);
%Optimize parameters x0
options = optimset('MaxFunEvals',2000);
[x,sse] = simulannealbnd(@get_sse,x0, lb, ub);
[x,sse] = fminsearch(@get_sse,x);
fval=[fval;sse];
%Replaces default values with new ones
result=get_param(x);
end

function sse = get_sse(x) %Takes in guess vector as input 
%Construct paramter vector
p=get_param(x);
disp(p);
%Construct model
result1=plot_model(p,0,1);
residual1=get_residual(result1,1);
residual1=sum(sum(residual1));

result2=plot_model(p,0,2);
residual2=get_residual(result2,2);
residual2=sum(sum(residual2));

sse=residual1+residual2;
end

function sse = get_sse2(x) %Takes in array as input


result1=plot_model(x,1,1);
residual1=get_residual(result1,1);

residual1=sum(sum(residual1));

result2=plot_model(x,1,2);
residual2=get_residual(result2,2);

residual2=sum(sum(residual2));

sse=residual1+residual2;
end

function residual=get_residual(r,light)
global M1 M2;
if light==1
    M=M1;
else
    M=M2;
end

temp_M=M(2:length(M(:,1)),:);
temp_r=r(2:length(M(:,1)),:);

residual=(temp_M-temp_r).^2;
end

function p=get_param(x)
global apply default;
p=zeros(7,length(default)-6);
ii=1;
for i=1:length(default(:,1))
    if apply(i)==0
        if i==7
            p(7,1+i-7)=default(i);
        elseif i>7
            p(7,1+i-7)=default(i)*default(7);
        else
            for iii=1:length(p(i,:))
            p(i,:)=default(i,:);
            end
        end
    elseif i==7
        p(7,1+i-7)=x(ii,1);
        ii=ii+1;
    elseif i>7
        p(7,1+i-7)=x(ii)*p(7,1+i-8);
        ii=ii+1;
    else
        p(i,:)=x(ii,1)*ones(1,length(default(1,:)));
        ii=ii+1;
    end    
end
end

function result=plot_model(x,make_plot,light)
h        = x(1,:);
k_a      = x(2,:);
k_d      = x(2,:);
syn_mRFP = x(3,:);
km       = x(4,:);
deg_mRFP = x(5,:);
syn_RFP  = x(6,:);
deg_RFP  = x(7,:);

global Init1 Init2 tspan1 tspan2;
if light==1
    Init=Init1;
    tspan=tspan1;
else
    Init=Init2;
    tspan=tspan2;
end
numtspan = length(tspan);

numi = length(Init(1,:));
result = zeros(numtspan,numi+1);
result2= zeros(numtspan,numi+1);

for i = 1:numi
    [t,y] = ode15s(@(t,y) SolveODE_param(t,y,h(i),k_a(i),k_d(i),syn_mRFP(i),km(i),deg_mRFP(i),syn_RFP(i),deg_RFP(i),light), tspan, Init(:,i));
    result(:,1)=t;
    result(:,1+i)=y(:,2); 
    if i==1
        result2(:,1)=t;
        result2(:,2)=y(:,1);
    else
    end
end

global M1 M2 S1 S2;
if make_plot==1
    if light==1
        M=M1;
        S=S1;
    else
        M=M2;
        S=S2;
    end
    figure;
    for ii=1:numi
    plot(result(:,1),result(:,ii+1),'linewidth',2); hold on;
    end
    set(gca,"ColorOrderIndex",1);  
    for iii=1:numi
    if size(M)==size(S)
        errorbar(M(:,1),M(:,iii+1),S(:,iii+1),'o','MarkerSize',5,'linewidth',2); hold on
    else
        plot(M(:,1),M(:,iii+1),'+','MarkerSize',10,'linewidth',2); hold on
    end
    legend('Cm','Dm','Am','Ym','Ce','De','Ae','Ye','Location','Northwest');
    end
    figure;
    plot(result2(:,1),result2(:,2),'linewidth',2); 
end
end

function dydt=SolveODE_param(t,y,h,k_a,k_d,syn_mRFP,km,deg_mRFP,syn_RFP,deg_RFP,light)
blue=is_blue(t,light);

EL2a= y(1);
RFP=  y(2);

mRFP= y(3);

dEL2a= k_a*blue-k_d*EL2a;
dRFP = syn_RFP*mRFP-deg_RFP*RFP;

dmRFP= syn_mRFP*(EL2a^h)/(km^h+EL2a^h)-deg_mRFP*mRFP;

dydt=[dEL2a;dRFP;dmRFP;];
end

function blue=is_blue(t,light)
if light==1
    if t<2
        blue=1;
    else
        blue=0;
    end
else
    if t<3
        blue=1;
    else
        blue=0;
    end
end
end