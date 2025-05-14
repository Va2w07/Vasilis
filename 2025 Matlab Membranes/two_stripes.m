% make the stripes longer
clear; clc; format long; % Class A -Lasers
Nx=256/2; Lx=80; Ny=256/2; Ly=80; non=1;
x=(Lx/Nx)*(-Nx/2:Nx/2-1)'; dx=x(2)-x(1); kx=(2*pi/Lx)*[0:Nx/2-1 -Nx/2:-1]';
y=(Ly/Ny)*(-Ny/2:Ny/2-1)'; dy=y(2)-y(1); [X,Y]=meshgrid(x,y);
ky=(2*pi/Ly)*[0:Ny/2-1 -Ny/2:-1]'; [Kx,Ky]=meshgrid(kx,ky);
L= -Kx.^2 - Ky.^2; field=zeros(100,Nx,Ny); mu=10; D=[]; pow=[];
% Define the square boundary limits
x_min = 30;  % Minimum x-coordinate for the square boundary
x_max = 80;  % Maximum x-coordinate for the square boundary
y_min = 30;  % Minimum y-coordinate for the square boundary
y_max = 80;  % Maximum y-coordinate for the square boundary
% Dynamic Distance Between the Two Potentials
distance = 0;  % Define distance between potentials (could vary dynamically)

rhoc = 4 * (1 + 0 * X); % a
rhb = 1 * (1 + 0 * X); % b
rhx1 = 30 * (1 + 0 * X); % b
rhx2 = -30 * (1 + 0 * X); % b
rh = 5 * (1 + 0 * X); % c

Dpump=10;
for tx=1:Nx
   for ty=1:Ny
       % Calculate radial distance for x and y axis
       rr_x = abs(x(tx));  % Distance in x direction
       rr_y = abs(y(ty));  % Distance in y direction
       rr_x1 = x(tx);  % Distance in x direction
       rr_y1 = y(ty);  % Distance in y direction
       rr = x(tx);  % third strip
       rr2=y(ty);
       % Set absorbing square boundary
      % if ( rr_x <= x_max && rr_y >= y_min && rr_y <= y_max)
       if ( rr_y <= y_max && rr_x >= x_min && rr_x <= x_max )
            po(tx, ty) = 1*i;    % Potential at the square boundary
            D0(tx, ty) = 0;    % Pump is zero inside the square boundary
       elseif ( rr_x <= x_max && rr_y >= y_min && rr_y <= y_max )
           po(tx, ty) = 1*i;    % Potential at the square boundary
           D0(tx, ty) = 0;    % Pump is zero inside the square boundary
        % Set strips
       elseif (rr_y1 > (rhb(tx, ty) + distance) && rr_y1 < ( rhoc(tx, ty) + distance) && rr_x1 < rhx1(tx, ty) && rr_x1 > rhx2(tx, ty))
           po(tx, ty) = 3.5 + 0.00001i*0.1;  % Potential inside the waveguide
           D0(tx, ty) = Dpump;        % Pump inside the waveguide
       elseif (rr_y1 < -(rhb(tx, ty) + distance) && rr_y1 > -(rhoc(tx, ty) + distance) && rr_x1 < rhx1(tx, ty) && rr_x1 > rhx2(tx, ty))
           po(tx, ty) = 3.5 + 0.00001i*0.1;  % Potential inside the waveguide
           D0(tx, ty) = Dpump;        % Pump inside the waveguide
       else
           po(tx, ty) = 3.5 - 100 * 1i;  % Loss outside the waveguide
           D0(tx, ty) = 0;       % No pump outside the waveguide
       end
   end
end


figure(1); surface(X,Y,real(po)); shading interp; axis square; axis tight; hold on;
title('potential (strips) Real part +pump +u0')
figure(2); surface(X,Y,imag(po)); shading interp; axis square; axis tight;
title('potential (strips) Imaginary part')
figure(1); surface(X,Y,D0); shading interp; axis square; axis tight; hold on;
%title('pump')

%%%%%%%%%%%%%%%%% dynamics
u = (exp(-((X+2.5).^2 + (Y-0).^2)/15));  % closer_strips A Gaussian-shaped initial field

%u = (exp(-((X+13).^2 + (Y-0).^2)/15));  % far_apart_strips A Gaussian-shaped initial field
%u = randn(Nx, Ny); % some input noise
v=fft2(u); S=1000; zmax=100; h=zmax/S; z=0; dim=0; count=0;
dim=0; savestep=10;

figure(1); surface(X,Y,10*abs(u).^2); shading interp; axis square; axis tight; hold on;
% figure(10); plot([-15,15], [12,12],'-m','LineWidth',3); hold on;
% plot([15,-15],[-12,-12],'-m','LineWidth',3); hold on;

for m=1:S
   % Reflective boundary conditions
%     u(1, :) = u(2, :);      % Left boundary (reflective)
%     u(Nx, :) = u(Nx-1, :);  % Right boundary (reflective)
%     u(:, 1) = u(:, 2);      % Bottom boundary (reflective)
%     u(:, Ny) = u(:, Ny-1);  % Top boundary (reflective)
   vm=v;  zm=z; % First step of Runge-Kutta
    par=exp(1i*z*L); g1=1i*exp(-1i*z*L);
  g2=fft2((po-non*1i*D0./(1+mu*abs(ifft2(v.*par)).^2)).*(ifft2(v.*par)));
   g=g1.*g2; z1=h*g;

  z=zm+0.5*h; v=vm+0.5*z1; % Second step of Runge-Kutta
  par=exp(1i*z*L); g1=1i*exp(-1i*z*L);
 g2=fft2((po-non*1i*D0./(1+mu*abs(ifft2(v.*par)).^2)).*(ifft2(v.*par)));
  g=g1.*g2; z2=h*g;

  z=zm+0.5*h; v=vm+0.5*z2; %  Third step of Runge-Kutta
  par=exp(1i*z*L); g1=1i*exp(-1i*z*L);
 g2=fft2((po-non*1i*D0./(1+mu*abs(ifft2(v.*par)).^2)).*(ifft2(v.*par)));
  g=g1.*g2; z3=h*g;

  z=zm+h; v=vm+z3; % Fourth step of Runge-Kutta
  par=exp(1i*z*L); g1=1i*exp(-1i*z*L);
 g2=fft2((po-non*1i*D0./(1+mu*abs(ifft2(v.*par)).^2)).*(ifft2(v.*par)));
  g=g1.*g2; z4=h*g;

 v=vm+(z1+2*z2+2*z3+z4)/6;  z=zm+h; % Final solution

  count=count+1; if (count==savestep)
  u=ifft2(v.*exp(1i*z*L)); dim=dim+1
  field(dim,:,:)=u; power(dim)=sum(sum(abs(u).^2))*dx*dy; dis(dim)=z;
  crosspower(dim)=sum(sum(abs(field(dim,:,127)).^2))*dx*dy;
  count=0; end;
end; u=ifft2(v.*exp(1i*z*L));

figure(4); surface(X,Y,abs(u).^2); shading interp; axis square; axis tight;hold on;
plot([-30,30], [30,30],'-m','LineWidth',3); hold on;
plot([30,-30],[-30,-30],'-m','LineWidth',3); hold on;
title('abs(u).^2 after Runge-Kutta') ;


% Create video writer
video = VideoWriter('an_stepbystep.mp4', 'MPEG-4');
video.FrameRate = 10; open(video);

fie=zeros(Nx,Ny); for tr=1:60
fie(:,:)=field(tr,:,:);
figure(5); surf(X,Y,abs(fie).^2); hold on;
plot([-32,32], [30,30],'-m','LineWidth',3); hold on;
plot([32,-32],[-30,-30],'-m','LineWidth',3); hold on;
title(['Time is ',num2str(tr),'sec']);
shading interp; view([-90 90]); axis tight; axis square;
view([0 90]); drawnow; pause(0.5);
frame = getframe(gcf);
   writeVideo(video, frame);end;
% Finalize video
close(video);
disp('Video saved as dotdot.mp4');