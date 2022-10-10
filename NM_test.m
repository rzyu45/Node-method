%% Calculate thermal dynamics using node method 
clear
clc
load("test_data.mat")

%parameter
parameter.density=960; %m^3/s
parameter.area=1.54; %m^2
parameter.Cp=4182; 
parameter.Ta=-10; % celsius degree
parameter.L=9250; %m
parameter.lambda=1/0.35;

%calculation parameter
T=20*3600;% total calculation time
opt.t=180;% temporal step size
N=T/t;% total steps

%variable
index.pipe=1;
num.pipe=1;
m=volume_per_second*parameter.density; %kg/s

%initial pipe instance
pipe_luhua=NMpipe(parameter,index,num,opt,m(1),boundary_condition(1));
Tout=zeros(N+1,1);
Tout(1)=initial_condition(end);

for t_count=1:N
    % push mass flow into the pipe
    pipe_luhua=pipe_luhua.push(m(t_count+1));
    % generate outlet temperature
    Tout(t_count+1)=pipe_luhua.generate_Tout();
    % update inlet temperature
    pipe_luhua=pipe_luhua.update_T(boundary_condition(t_count+1));
end

%Compare the results with the measured data
plot(Tout_measured,'k')
hold on
plot(Tout,'Color','r')
hold off