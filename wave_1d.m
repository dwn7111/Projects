clear;
%%Defining Space
Lx = 20;
dx = 0.1;
nx = fix(Lx/dx);
x = linspace (0, Lx, nx);

%%Defining Time
T = 20;

%%Defining Variable
wn = zeros(nx, 1);
wnm = wn; %%Wn at time n-1
wnp = wn; %%Wn at time n+1

%%Parameter
CFL = 1; %%CFL = v*dt/dx
v = 1;
dt = (CFL*dx)/v;

%%Initial Condition
wn(49:51)=[0.1 0.2 0.1];
wnp1(:) = wn(:);

%%Looping
t = 0;
while t<T
    %Either the boundary will reflect the wave or absorb, we have both
    %conditions
    
    %Refleting Boundary condition
    wn([1 end]) = 0;
    
    %Absorbing Boundary Condition
    %wnp1(1) = wn(2) + ((CFL -1)/(CFL +1))*(wnp1(2)-wn(1));
    %wnp1(end) = wn(end-1) + ((CFL -1)/(CFL +1))*(wnp1(end-1)-wn(end));
    
    %solution
    t = t+dt;
    wnm1 = wn;
    wn = wnp1; %Saving our arrays
    
    %Source
    %wn(50) = dt^2 + 20*sin(20*pi*t/T);
    
    for i = 2:nx-1
        wnp1(i) = 2*wn(i) - wnm1(i) + CFL^2*(wn(i+1)) - 2*wn(i)+ wn(i-1);
    end
  
    clf;
    plot(x,wn);
    title (sprintf('t = %.2f',t));
    axis([0 Lx -0.5 0.5]);
    shg;pause(0.01);
end
wn(end)
x(end)