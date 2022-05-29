clear;
% This program solves the shallow water equation for a singular drop
% falling on surface of water
%% defining space and other parameters
n = 50;                  % grid
g = 9.8;                 % gravitational accleration
dx = 1.0;
dy = dx;
dt = 0.01;   % you can change the speed of simulaiton by 
              % changing the dt to a smaller value

T = 25;                  % time limit
d = drop(1.5,21);        % drop's dimensions

b = 0.5;

[surfplot,top] = initgraphics(n); 

while 1==1 %This loop restarts the solving, after time limit is reached at t = T
%     pause(0.5) 
    % defining Variable (Zeroes)
    h = ones(n+2,n+2); 
    hx = zeros(n+1,n+1); 
    hy = zeros(n+1,n+1);
    
    U = zeros(n+2,n+2);  
    Ux = zeros(n+1,n+1); 
    Uy = zeros(n+1,n+1); 
   
    V = zeros(n+2,n+2);
    Vx = zeros(n+1,n+1);
    Vy = zeros(n+1,n+1);

    %% Solution Loop
    
    %Looping conditionals
    t = 0;
    s = 0;
    ss = 0;
    %loop
    while t<T
       t = t+dt;                 % Loop counter
       ss = ss+1;
       % defining water drop position
       % This loop conditions the drop to drop only oncce per cycle of
       % solving
       if  s <2
           q1 = rand;q2 = rand;
           a = 21;                          % taken from width of the drop
           i = ceil(q1*25)+(1:a);           % x coordinate of the center of the drop
           j = ceil(q2*29)+(1:a);           % y coordinate of the center of the drop
           h(i,j) = h(i,j) + b*d;           % max height of the center of the drop, it can be changed by changing value of b above
           s = s+1;
       end
       % position of the drop on the grid can be changed by changing q1 and
       % q2 between 0 & 1
       

       % Refleting Boundary conditions
       h(:,1) = h(:,2);      U(:,1) = U(:,2);       V(:,1) = -V(:,2);
       h(:,n+2) = h(:,n+1);  U(:,n+2) = U(:,n+1);   V(:,n+2) = -V(:,n+1);
       h(1,:) = h(2,:);      U(1,:) = -U(2,:);      V(1,:) = V(2,:);
       %barrier
       h(n+2,:) = h(n+1,:);  U(n+2,:) = -U(n+1,:);  V(n+2,:) = V(n+1,:); 

       %% 
       % Solving first half step
       % In x
       i = 1:n+1;
       j = 1:n;

       hx(i,j) = (h(i+1,j+1)+h(i,j+1))/2 - dt/(2*dx)*(U(i+1,j+1)-U(i,j+1));
       
       Ux(i,j) = (U(i+1,j+1)+U(i,j+1))/2 -  ...
                 dt/(2*dx)*((U(i+1,j+1).^2./h(i+1,j+1) + g/2*h(i+1,j+1).^2) - ...
                            (U(i,j+1).^2./h(i,j+1) + g/2*h(i,j+1).^2));
                        
       Vx(i,j) = (V(i+1,j+1)+V(i,j+1))/2 - ...
                 dt/(2*dx)*((U(i+1,j+1).*V(i+1,j+1)./h(i+1,j+1)) - ...
                            (U(i,j+1).*V(i,j+1)./h(i,j+1)));
       %Same In y
       i = 1:n;
       j = 1:n+1;
       
       hy(i,j) = (h(i+1,j+1)+h(i+1,j))/2 - dt/(2*dy)*(V(i+1,j+1)-V(i+1,j));

       Uy(i,j) = (U(i+1,j+1)+U(i+1,j))/2 - ...
                 dt/(2*dy)*((V(i+1,j+1).*U(i+1,j+1)./h(i+1,j+1)) - ...
                            (V(i+1,j).*U(i+1,j)./h(i+1,j)));
       Vy(i,j) = (V(i+1,j+1)+V(i+1,j))/2 - ...
                 dt/(2*dy)*((V(i+1,j+1).^2./h(i+1,j+1) + g/2*h(i+1,j+1).^2) - ...
                            (V(i+1,j).^2./h(i+1,j) + g/2*h(i+1,j).^2));
               
       % Second half step
       i = 2:n+1;
       j = 2:n+1;

      
       h(i,j) = h(i,j) - (dt/dx)*(Ux(i,j-1)-Ux(i-1,j-1)) - ...
                         (dt/dy)*(Vy(i-1,j)-Vy(i-1,j-1));
      
       U(i,j) = U(i,j) - (dt/dx)*((Ux(i,j-1).^2./hx(i,j-1) + g/2*hx(i,j-1).^2) - ...
                         (Ux(i-1,j-1).^2./hx(i-1,j-1) + g/2*hx(i-1,j-1).^2)) ...
                       - (dt/dy)*((Vy(i-1,j).*Uy(i-1,j)./hy(i-1,j)) - ...
                         (Vy(i-1,j-1).*Uy(i-1,j-1)./hy(i-1,j-1)));
      
       V(i,j) = V(i,j) - (dt/dx)*((Ux(i,j-1).*Vx(i,j-1)./hx(i,j-1)) - ...
                         (Ux(i-1,j-1).*Vx(i-1,j-1)./hx(i-1,j-1))) ...
                       - (dt/dy)*((Vy(i-1,j).^2./hy(i-1,j) + g/2*hy(i-1,j).^2) - ...
                         (Vy(i-1,j-1).^2./hy(i-1,j-1) + g/2*hy(i-1,j-1).^2));

      
       % Values of variables are stored and updated
       if ss == 10   
       % The value of ss can be used control the speed of dropping the drop and subsequent plot
       
          k = abs(U(i,j)) + abs(V(i,j));
          set(surfplot,'zdata',h(i,j),'cdata',k);
          set(top,'string',sprintf('t = %.2f',t))
          drawnow
          ss=0;
       end
    end
end
%% Functions
%drop Function
function d = drop(h,w) %h - height , w - width
   [x,y] = ndgrid(-1:(2/(w-1)):1);
   d = h*exp(-5*(x.^2+y.^2)); % drop's Shape is taken to be something like a Gaussian
end
  %%
% Surf plot of the function
function [surfplot,top] = initgraphics(n)
   clf
   x = (0:n-1)/(n-1);
   surfplot = surf(x,x,ones(n,n),zeros(n,n));
   axis([0 1 0 1 -1 3])
    caxis([0 3])
   shg
   colorbar;
   top = title('Shallow Water');
end