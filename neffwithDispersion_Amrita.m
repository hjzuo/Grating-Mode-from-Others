% This example computes all of the field components of a
% Ridge waveguide on SiO2 sub and coated by IPG

% %                 IPG cladding 
% % 
% %         ^       -------------      ^
% %         |       |           |       |  rh 
% %         h2      |           |       |
 %          | -------             ---------
% %         |         TeO2 film
% %         -----------------------------
% %         ^
% %         |
% %         h1
% %         |          SiO2 sub
% %         |
% %         -------------------------------
% % 
clear all;
clf;
P1 = [7.49717782E-20 -5.08979494E-16 1.42386851E-12 -2.10794779E-9 1.75069041E-6 -7.91674295E-4 1.61107043E+0]; %SiO2
%P2 = [2.77866191E-13 -1.63422779E-9 3.67920538E-6 -3.82404897E-3 4.22689944E+0]; % Ge11
%P2 = [-2.251936E-16 1.690474E-12 -5.152364E-9 8.040037E-6 -6.531151E-3 4.957051E+0]; % Ge11.5 for air clad
P2 = [1.974996E-13 -1.219626E-9 2.895649E-6 -3.189796E-3 4.095021E+0];  %annealed
%P2 = [2.631e-19 -2.147e-15 7.340e-12 -1.353e-8 1.431e-5 -8.370e-3 4.554] %%%As2S3 not annealled Michael data.
P3 = [7.80953561E-14 -3.74069408E-10 6.69117335E-7 -5.39302445E-4 1.68215884E+0]; %IPG
stepsize = 20;
%startWavelength = 1400;
%endWavelength = 1500;
%lambda = 1540;
lambda=[1000:stepsize:1600];
nmodes = 2;
c = 3e-7; % c value for D in ps.nm/km
c1 = 3e-4; % c value in m/ps
omega = 2*pi*c1./(lambda.*10^-9);
%omega = linspace(1,1000,length(lambda));
effectiveIndex = zeros (length(lambda), nmodes);
width = [1000:1000:4000]; 
etchdepth = [0.3:0.1:0.7];
height=[500:100:800];
cladding_index = [1:0.25:2.5];
  %for nii = 1:length(cladding_index);
      %n3 = cladding_index(nii);
%  for zii = 1:length(height);
%       h2 = height(zii);
 %for yii = 1:length(width);
      %w = width(yii);
%for xii= 1:length(etchdepth);
     %ed = etchdepth(xii);
for ii=1:length(lambda)
    wavelength = lambda(ii);
    n3 = 1.0008; % refractive index of air cladding
    n1 = polyval(P1,wavelength);          % SiO2 lower cladding
    n2 = polyval(P2,wavelength);           % Ge11 core
    %n2 = 2;
    %n3 = polyval(P3,wavelength);          % IPG upper cladding
    
    h1 = 1462.4;           % SiO2 thickness(nm)
    h2 = 686.2;           % Ge11 film thick (nm)
    h3 = 2000;           % IPG thickess cladding (nm)
    rh = 0.6*h2;            % Ridge heigh

    dx = 50;             % grid size (x)
    dy = dx;             % grid size (y)

    %lambda = 1400;       % wavelength (nm)
    %nmodes = 8;          % number of modes to compute
    w = 2000;            % waveguide full-width (nm)
    side = 2000;         % space on side of waveguide (nm)

    fprintf (1,'generating index mesh...\n');

    [x,y,xc,yc,nx,ny,eps] = waveguidemeshfull([n1,n2,n3],[h1,h2,h3],rh,w/2,side,dx,dy);

    [x1,y1,xc1,yc1,nx1,ny1,eps1] = waveguidemeshfull([n1,n2,n3],[h1,h2,h3],rh,w/2,side,dx,dy);


    [x2,y2,xc2,yc2,nx2,ny2,eps2] = waveguidemeshfull([n1,n2,n3],[h1,h2,h3],rh,w/2,side,dx,dy);
    % % % Now we stretch out the mesh at the boundaries: N S E W
    [x,y,xc,yc,dx,dy] = stretchmesh(x,y,[20,20,40,40],[1.5,1.5,1.5,1.5]);
    fprintf (1,'solving for eigenmodes...'); t = cputime;

    [Hx,Hy,neff] = wgmodes (wavelength, n2, nmodes, dx, dy, eps, '0000');
    fprintf (1,'done (cputime = %7.3f)\n', cputime-t);
    wavelength
    neff
    effectiveIndex(ii,:)=neff;
end
    for i = 1:nmodes
        for j = 1:length(lambda)
            effectiveIndex1 (j,:) = -sort(-effectiveIndex(j,:));
        end
    end
   % effectiveIndex2 = effectiveIndex1;
%         for j = 1:nmodes
%             effectiveIndex2 (:,j) = -sort(-effectiveIndex1(:,j));
%     end
effectiveIndex2=effectiveIndex1;
    figure(1); clf; plot (lambda,effectiveIndex2, 'x', lambda, effectiveIndex2); 
    %figure(rh); plot (lambda, abs(effectiveIndex1), 'x', lambda, abs(effectiveIndex1));
     xlabel ('wavelength');
     ylabel('Neff');
     title1=sprintf('neff of Ge11 waveguide %dnm width, %dnm height, %dnm ridge %dnm clad', w, h2, rh, n3);
     title(title1);
     figureindex = sprintf('Ge11_5IPGClad/Neff %dnm height %dnm ridge % dnm width %dnm clad',h2, rh, w, n3);
     saveas(gcf, figureindex,'jpg');
     clf;
    

%         P = polyfit(lambda',effectiveIndex2(:,1),4)
%         disp = (12*P(1).*lambda.^2 + 6*P(2).*lambda+ 2*P(3));
%         disp = disp.*lambda/c;
           
   for j = 1:nmodes;
        for i = 1:length(lambda)-1;
            dndl(i,j) = (effectiveIndex2(i+1,j)- (effectiveIndex2(i,j)));
        end
    end

    for j = 1:nmodes;
        for i = 1:length(lambda)-2;
            d2ndl(i,j) = (dndl(i+1,j) - dndl(i,j));
        end
    end
    lambdanew = lambda(2:length(lambda)-1); 

    for b = 1: nmodes
        for a = 1:length(lambdanew)
            D(a,b) = -lambdanew(1,a)*d2ndl(a,b)/(stepsize^2);
        end
    end
    D = D/c;
%     for aii = 1:nmodes
%     disp(xii,aii) =D(floor(length(lambdanew)/2),aii);
%     end 
        
 figure(100); plot (lambdanew', D, 'x', lambdanew', D);
 xlabel('Wavelength'); ylabel('D');
 title2=sprintf('GVD of Ge11 waveguide %dnm width, %dnm height, %dnm ridge, %dnm clad', w, h2, rh, n3);
 title(title2);
 figureindex1 = sprintf('Ge11_5IPGClad/GVD %dnm height %dnm ridge %dnm width %dnm clad',h2, rh, w, n3);
 saveas(gcf, figureindex1,'jpg');
 clf;
%  figure(3); plot(lambda',disp);
%  temp1 = diff(abs(effectiveIndex2));
%  temp2 = diff(temp1);
%  for yii = 1:nmodes
%  GVD(:,yii) = (1/c)*lambdanew'.*temp2(:,yii)/(stepsize^2);
%  figure(3); plot(lambdanew',GVD(:,yii),'x',lambdanew',GVD(:,yii));
%  title(rh);
%  hold on;
%  end
%  hold off;
% end 
% end
%  end
for bii = 1;
beta(:,bii) = effectiveIndex2(:,bii).*(omega'/c1);
dbdo = diff(beta(:,bii))./diff(omega');
omega1 = diff(omega)/2+omega(1:length(omega')-1);
beta2 = diff(dbdo)./diff(omega1');
omega2 = diff(omega1)/2 + omega1(1:length(omega1)-1);
beta3 = diff(beta2)./diff(omega2');
omega3 = diff(omega2)/2 + omega2(1:length(omega2)-1);
figure(200);plot(omega2,beta2, 'x', omega2, beta2);
end
