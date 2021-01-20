%    PROGRAM wtrrkt
%
%    this code solves for the thrust and trajectory
%    histories for a water rocket.  code written by
%    s. heister during week of 4/1-4/5/96.
%
      global V AA cd At po p pa rho eta h Vuo f ve cdrag rhoa m n me

%
%    constants
%
%     print*, "input init p_ull, poly. const,n & noz diam"
%     read(*,*) po,n,dthroat

	  format compact
	  fract=[.125,.195,.2,.21,.22.23,.24,.25,.26,.27,.28.29,.30,.31,.32,.33,.34,.35,.39,.52]
	  po=input('Input init p_ull  ')
%	  n=input('Input polynomial constant, n  ')
	  dthroat=11/32
%	  dthroat=input('Input nozzle diameter  ')
      dt=0.00025                  %time step, sec
      g=32.174                    %gravity, f/s^2
      rho = 62.4/1728.            %density of h2o, lb/in^3
      me=6/16.                    %rocket empty mass, lb
      n=1.4                       %polytropic constant   
%     pi=4.*atan(1.);             %pi
      AA=(pi*3.5^2)/4             %projected area, in^2
      cd=0.9                      %discharge coeff.
      cdrag=0.4                   %drag coefficient
      rhoa=0.0765                 %air density, lb/f^3
      pa=14.7                     %ambient pressure, psi
      At=0.25*pi*dthroat^2        %nozzle throat area, in^2
      V=2.*0.03532*1728.          %bottle volume, in^3
	  itmax=10000				  %number of time points in simulation
	  iprint=50					  %print interval (print every iprint)
	  iff=length(fract)			  % number of fuel fractions
	  fuel=zeros(iff,3);
	  ipflag=0					  %ipflag=1 to print time history
%
%    set initial conditions
%

% >>>  Start MAIN LOOP on fuel fraction
for jj=1:iff
	  ff=fract(jj)
%     ff=0.09+0.01*jj	  %fuel fraction
      t=-dt;
      Vuo=V*(1-ff)		  %initial air fraction
      vu=Vuo;			  %vu is air volume, in^3
      m=me+rho*(V-vu);	  %weight of fueled rocket, lb
      alt=0.0;            %rocket altitude, ft.
      h=(V-vu)/AA;        %fluid height assuming cyl. tank
      eta=5.;

		for ii=1:10
		      p=po*(Vuo/vu)^n;
		      ve=sqrt(2.*g*(p-pa+rho*eta*h)/(12.*rho*(1.-(At/AA)^2))) ;              
		      f=rho*cd*ve^2*At*12./g;        %thrust, lb
		      eta=f/m;                       %vehicle acceleration, g's
		      vr=0.0;                        %rocket velocity, f/s
		      alt=0.0;                       %rocket altitude, f
		      timp=0.0;                      %total impulse, lb-sec
		end

		if ipflag==1
%      		Print header
			disp('time,sec,   ull.vol   p_ullage   thrust,lb   m_rkt,lb   vr,f/s  alt,ft')
		end

%    begin time stepping.....
%
      	fsave=f;
		out=zeros((round(itmax/iprint+1)),7);
		kk=0;
		for ii=1:itmax
		      t=t+dt;
		      y1=dvudt(vu);
		      z1=accel(vr,vu);
		      timp=timp+0.5*dt*(fsave+f);
		      fsave=f;
		%     if(ii.eq.1.or.mod(ii,iprint).eq.0) write(*,120) t,vu/V,p,f,m,vr,alt
			  if ii==1 | rem(ii,iprint)==0
			    	kk=kk+1;
			  		out(kk,:)=[t,vu/V,p,f,m,vr,alt];
					if ipflag==1; disp(out(kk,:)); end
			  end
			
		      vu2=vu+0.5*dt*y1;
		      vr2=vr+0.5*dt*z1;
		      y2=dvudt(vu2);
		      z2=accel(vr2,vu2);                 
		      vu3=vu+0.5*dt*y2;
		      vr3=vr+0.5*dt*z2;
		      y3=dvudt(vu3);
		      z3=accel(vr3,vu3);                 
		      vu4=vu+dt*y3;
		      vr4=vr+dt*z3;
		      y4=dvudt(vu4);
		      z4=accel(vr4,vu4);               
		      vu=vu+dt*(y1+2.*y2+2.*y3+y4)/6.;
		      alt=alt+dt*(vr+2.*vr2+2.*vr3+vr4)/6.;
		      vr=vr+dt*(z1+2.*z2+2.*z3+z4)/6.;
		      if vr<0.0; break; end
		end
		% end time stepping
		%save fuel fraction results
		fuel(jj,:)=[ff,alt,timp];
		disp(' ')
		disp(fuel(jj,:))
		disp(' ')
		plot(out(:,1),out(:,6))
		xlabel('time, sec')
		ylabel('velocity, ft/sec')
		str=['Velocity time history for a fuel fraction of ',num2str(ff)];
		title(str)
%
end
% >>>  End MAIN LOOP
fuel
plot(fuel(:,1),fuel(:,2))
xlabel('initial fuel fraction')
ylabel('final altitude, ft')
title('Effect of fuel fraction on altitude')
save third fuel -ascii