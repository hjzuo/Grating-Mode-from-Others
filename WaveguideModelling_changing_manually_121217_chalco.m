% This example computes all of the field components of a
% Ridge waveguide
% %         ^        coat (n3)
% %         |
% %         | h3
% %         ^         -------------       ^
% %         |         |           |       |  rh 
% %         h2        |           |       |
% %         |   -------            ---------
% %         |         core film (n2)
% %         -   ----------------------------
% %         ^
% %         |
% %         h1
% %         |          sub (n1)
% %         |
% %         -   ------------------------------
% % 
clear all;
close all
% % starttime=datestr(now)
starttime=clock;
timestamp=datestr(now, 30);


% % SPECIFY THE RUNNING RANGE HERE
% % -----------------------------------
dx0=25;                 %Grid size nm
dy0=25;                %grid size nm
side = 2000;            % space on side of waveguide (nm) matters
h1 = 5000;              % buffer thickness(nm)
h3 = 1000;              % coating  thickess - cladding (nm)

nmodes=2;
lambda=3900:100:4900;      % Wavelength range
w_array=4000;         % waveguide width array
h2_array=3000;         % film thickness 
rh_por_array=0.5;       % etch depth portion to thickness (move in side film thickness)
% % -------------------------------------

iw_length=length(w_array);
il_length=length(lambda);
ih2_length=length(h2_array);
irh_length=length(rh_por_array);

% % % Initialising matrix to store results
All_eff_index=zeros(il_length, nmodes, iw_length, ih2_length, irh_length);
All_diff_index=zeros(il_length-1, nmodes, iw_length, ih2_length, irh_length);
All_GVD=zeros(il_length-2, nmodes, iw_length, ih2_length, irh_length);


for i=1:length(lambda)
    
% mgf2_sellmeier;
 n_mgf2(i)=sqrt(1+0.48755108*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.04338408^2)+0.39875031*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.09461442^2)...
   +2.3120353*(lambda(i)/1000)^2/((lambda(i)/1000)^2-23.793604^2));%200-7000nm
 n1(i)=n_mgf2(i);
%  
% baf2_sellmeier;
%  n_baf2(i) =sqrt( 1 + 0.643356*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.057789^2) + 0.506762*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.10968^2)...
%     + 3.8261*(lambda(i)/1000)^2/((lambda(i)/1000)^2-46.3864^2));%270-10300nm
%  n1(i)=n_baf2(i);
 
 %  n_sio2(i)=sqrt(1+0.663044*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.060^2)+0.517852*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.106^2)...
%          +0.175912*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.119^2)+0.565380*(lambda(i)/1000)^2/((lambda(i)/1000)^2-8.844^2) ...
%          + 1.675299*(lambda(i)/1000)^2/((lambda(i)/1000)^2-20.742^2));%180-3000nm
 %n1(i)=n_sio2(i);
 
 %as2s3_sellmeier;
%   as2s3(i)=1+1.898368*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.0225)+1.922298*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.0625)...
%            +0.8765138*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.1225)+0.118878*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.2025)...
%            +0.956998*(lambda(i)/1000)^2/((lambda(i)/1000)^2-749.99847);
%     n_as2s3(i)=sqrt(as2s3(i));   %some paper
%  as2s3(i)=1+4.72553*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.25920^2)+0.46765*(lambda(i)/1000)^2/((lambda(i)/1000)^2-22.33309^2);
%  n_as2s3(i)=sqrt(as2s3(i)); 
 %n2(i)=n_as2s3(i);  %c zhiyong yang
   %ge11 (GeAsS) sellmeier
%  ge11(i)=1+4.18011*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.31679^2)+0.35895*(lambda(i)/1000)^2/((lambda(i)/1000)^2-22.77018^2);
%    n_ge11(i)=sqrt(ge11(i));
%    %n4(i)=n_ge11(i);
%     n1(i)=n_ge11(i); 
    %ge12.5_sellmeier
 ge12(i)=1+5.78525*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.28795^2)+0.39705*(lambda(i)/1000)^2/((lambda(i)/1000)^2-30.39338^2);
 n_ge12(i)=sqrt(ge12(i));
 n2(i)=n_ge12(i);
 
    %silicon_sellmeier
%     n_silicon(i) = sqrt( 1 + 10.6684293*(lambda(i)/1000)^2/((lambda(i)/1000)^2-0.301516485^2) + 0.003043475*(lambda(i)/1000)^2/((lambda(i)/1000)^2-1.13475115^2)...
%         + 1.54133408*(lambda(i)/1000)^2/((lambda(i)/1000)^2-1104.0^2) );
%     n2(i)=n_silicon(i);%1360-11000nm;
    


%teflon_sellmeier c duk;
n_teflon(i)=sqrt(1+0.8628*(lambda(i)/1000)^2)/((lambda(i)/1000)^2-0.1173^2)+0.0958*(lambda(i)/1000)^2/((lambda(i)/1000)^2-7.3196^2);%1700-5600nm
n3(i)=n_teflon(i);

% n4(i)=1;% air cladding;
end;   
    
    % %   CHECKING MATERIAL DISPERSION
    last=length(lambda);
    lambda2=lambda;
    lambda2(last)=[];
    lambda2(1)=[];
    
     mat_dis1=diff(n2);
     GVD=-(10^7/3)*lambda2.*diff(mat_dis1)/(lambda(2)-lambda(1))^2; 

     figure(9)
     clf;
     set(gca,'FontSize',20);
     plot(lambda2, GVD,'LineWidth',2);
     hold
     x=get(gca,'xlim');
     y=0;
     hold on
     plot(x,[y y],'r');
      
    xlabel('Wavelength(nm)','FontSize',20);
    ylabel('Dispersion parameter (ps/km.nm)','FontSize',20);
    title1=sprintf('Dispersion parameter of the film');
    title(title1,'FontSize',20);
    name=sprintf('Dispersion parameter of the films_%s', timestamp);
    saveas(gcf,name, 'png')
 

% %     START THE LOOP TO CALCULATE neff FOR DIFFERENT WAVELENGTH, RIDGE
% % % % AND WIDTH
count=0;
effectiveIndex = zeros(length(lambda), nmodes);

for iw = 1:iw_length;            % waveguide full-width (nm)

    w=floor(w_array(iw))
for ih2=1:ih2_length;           %% Looping thickness
    h2=floor(h2_array(ih2))
    
    for irh =1:irh_length;            % Ridge heigh
      rh=floor(h2*rh_por_array(irh))
          for ii=1:il_length
            count=count+1;
            dx = dy0;             % grid size (x)
            dy = dx0;             % grid size (y) 
    %         fprintf (1,'generating index mesh...\n');

            [x,y,xc,yc,nx,ny,eps] = wgmesh([n1(ii),n2(ii),n3(ii)],[h1,h2/2,h3*3],rh,w/2,side,dx,dy);
%             waveguidemeshfull([n1(ii),n2(ii),n3(ii)],[h1,h2,h3],rh,w/2,side,dx,dy);
             

            % % % Now we stretch out the mesh at the boundaries: N S E W
            [x,y,xc,yc,dx,dy] = stretch(x,y,[20,20,40,40],[1.5,1.5,1.5,1.5]);
    %         fprintf (1,'solving for eigenmodes...'); t = cputime;

            [Hx,Hy,neff] = wgmodes (lambda(ii), n2(ii), nmodes, dx, dy, eps, '0000');
             effectiveIndex(ii,:)=sort(neff,'descend');
          end;

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % 
    % % Ordering the right neff base on 2nd order derivatives or GVD

% %     effectiveIndex2=smooth(effectiveIndex);
    effectiveIndex2=effectiveIndex;
    
    All_eff_index(:,:,iw,ih2, irh)=effectiveIndex2;
    % % Ploting output
    figure(1)
    clf;

% %     figure(2)
% %     clf;

    figure(3)
    clf;
    
%     figure(4)
%     clf;

    % % Wavelength for 1st and 2nd derivative
    lambda1=lambda;
    last=length(lambda);
    lambda1(last)=[];
    lambda1=lambda1+(lambda(2)-lambda(1))/2;

    lambda2=lambda;
    lambda2(last)=[];
    lambda2(1)=[];
   
        for ii=1:nmodes
         figure(1)
         hold all;
         set(gca,'FontSize',20);

           plot(lambda, effectiveIndex2(:,ii), '*-','LineWidth',2);

           % [p, s, mu] = polyfit(lambda',effectiveIndex2(:,ii), 7);
%             [effIndex_fit, delta]=polyval(p, lambda', s, mu);
%             delta=effIndex_fit-effectiveIndex2(:,ii);
% %             plot(lambda, effIndex_fit);       
            
% % Calculating GVD=(lambda/c)*d^2n/dl^2ex_fit);    
% %             ng=diff(effIndex_fit);
            ng=diff(effectiveIndex2(:,ii));
            GVD=-(10^7/3)*lambda2'.*diff(ng)/(lambda(2)-lambda(1))^2;
            All_diff_index(: ,ii, iw, ih2, irh)=ng;
            All_GVD(:,ii, iw,ih2, irh)=GVD;
    
% %         figure(2)
% %         hold all
% %         plot(lambda1,ng,'-*')

        figure(3)
        hold all
        plot(lambda2, GVD,'*-','LineWidth',2); 
        
%         figure(4)
%         hold all 
%         plot(lambda, delta);

        end

    figure(1)
    title1=sprintf('neff  %dnm wide, %dnm thick, %dnm etched', w, h2, rh);
    title(title1);
    xlabel('Wavelegnth (nm)');
    ylabel('Effective indices');
    legend('mode 1', 'mode 2', 'mode 3', 'mode 4')
    legend('boxoff')
    
    name=sprintf('WG_%dnm_%dnm_%dnm_neff_%s',w, h2, floor(rh), timestamp);
    saveas(gcf,name, 'png')

% % %     figure(2)
% % %     title1=sprintf('div of neff of TeO2 waveguide %dnm wide, %dnm thick, %dnm etched', w, h2, rh);
% % %     title(title1);
% % %     xlabel('Wavelegnth (nm)');
% % %     ylabel('del neff');
% % %     legend('mode 1','mode 2','mode 3','mode 4')
% % %     
% % %     name=sprintf('TeO2_%dnm_%dnm_%dnm_neff_diff_%s',w, h2, floor(rh), timestamp);
% % %     saveas(gcf,name, 'png')

    figure(3)
    set(gca,'FontSize',20);

    plot([lambda(1) lambda(last)], [0, 0]);
%     plot([1550, 1550], [-200, 200]);
    title1=sprintf('Dispersion parameter %dnm wide, %dnm thick, %dnm etched', w, h2, rh);
    title(title1);
    xlabel('Wavelegnth (nm)');
    ylabel('Dispersion parameter (ps/km.nm)');
  ylim([-150 400]);
    legend('mode 1', 'mode 2','Location','NorthWest')
     legend('boxoff')
  
    name=sprintf('WG_%dnm_%dnm_%dnm_GVD_%s',w, h2, floor(rh), timestamp);
    saveas(gcf,name, 'png')

 
% %     SAVING MATLAB DATA FILE
    name=sprintf('WG_%dnm_%dnm_%dnm_%s',w, h2, floor(rh), timestamp);
    save(name)
    
    
% % % Estimate time to finish     
    esendtime=(1/60)*etime(clock,starttime)/count*(il_length*ih2_length*irh_length*iw_length-count);

    fprintf (1,'%d wavelengths of %d done, time to finish %.2f mins (%.2f hours)\n',count, il_length*ih2_length*irh_length*iw_length, esendtime, esendtime/60);

    
    end;        %ridge
    end;        %thickness
end;            %width

% % Saving final result arrays
    name=sprintf('WG_All_summary_%s', timestamp);
    save(name,'nmodes','lambda', 'w_array','rh_por_array', 'All_eff_index', 'All_diff_index', 'All_GVD')
    
    
% % %     Plotting 1 variable at the time, keeping other at centre
% %  GVD
% % figure(300)
% % clf;
% % for ii=1:nmodes
% %     d=size(All_GVD);
% %     s=floor(1/2*d);
% %     plot_val=[];
% %     
% %     
% %      s(1)=s(1)+floor((1550-lambda(s(1)))/(lambda(2)-lambda(1)))-1;            %wavelength
% %      s(3)=d(3)-1;        %width
% %     s(4)=d(4)-2;      %film thick
% %     
% %     for i_index=1:d(5)  %etched depth
% %         
% %         j_index=(i_index-1)*d(4)*d(3)*d(2)*d(1) + (s(4)-1)*d(3)*d(2)*d(1)...
% %                 +(s(3)-1)*d(2)*d(1) + (ii-1)*d(1) + s(1)+2;
% %             
% % % %         j_index=(s(5)-1)*d(4)*d(3)*d(2)*d(1)+(s(4)-1)*d(3)*d(2)*d(1)...
% % % %                 +(s(3)-1)*d(2)*d(1) + (s(2)-1)*d(1) + s(1)+2
% % 
% % % %         j_index=d2(1)+2+(ii-1)*d(1)+(d2(3)-1)*d(2)*d(1)+(i_index-1)*d(3)*d(2)*d(1);
% % 
% %         plot_val(i_index)=All_GVD(j_index);
% %         
% %     end;
% %     figure(300)
% %     hold all
% %     plot(rh_por_array*100, plot_val,'*-')
% % 
% % end
% %     figure(300)
% %     xlabel('Etched depth (%)');
% %     ylabel('Dispersion parameter (ps/km.nm)');
% % %     ylim([-200 200])
% %     title1=sprintf('Dispersion parameter %dnm wl %dnm wide, %dnm thick', lambda(s(1)+2), w_array(s(3)), h2_array(s(4)));
% %     title(title1);
% %     legend('mode 1', 'mode 2', 'mode 3', 'mode 4')
% %     name=sprintf('WG_GVD_ridge_%dnm_%dnm_%dnm_ridge_%s', lambda(s(1)+2), w_array(s(3)), h2_array(s(4)), timestamp);
% %     saveas(gcf,name, 'png')
% %     
% % 
% % figure(400)
% % clf;
% % for ii=1:nmodes
% %     d=size(All_GVD);
% %     s=floor(1/2*d);
% %     plot_val=[];
% % %     
% %      s(1)=s(1)+floor((1550-lambda(s(1)))/(lambda(2)-lambda(1))) -1;            %wavelength
% %     s(3)=d(3)-1;      %width
% % %     s(5)=3;         %etched depth
% % 
% % 
% %     for i_index=1:d(4)    %thick
% %         j_index=(s(5)-1)*d(4)*d(3)*d(2)*d(1)+(i_index-1)*d(3)*d(2)*d(1)...
% %                 +(s(3)-1)*d(2)*d(1) + (ii-1)*d(1) + s(1)+2;
% % % %         j_index=d2(1)+2+(ii-1)*d(1)+(i_index-1)*d(2)*d(1)+(d2(4)+1-1)*d(3)*d(2)*d(1);
% %         plot_val(i_index)=All_GVD(j_index);
% %     end;
% %     figure(400)
% %     hold all
% %     plot(h2_array, plot_val,'*-')
% % 
% % end
% %     figure(400)
% %     xlabel('Film thick (nm)');
% %     ylabel('Dispersion parameter ps/km.nm');
% %      ylim([-200 200])
% %     title1=sprintf('Dispersion parameter %dnm wl, %dnm wide, %d percent etched', lambda(s(1)+2), w_array(s(3)), floor(rh_por_array(s(5))*100));
% %     title(title1);
% %     legend('mode 1', 'mode 2', 'mode 3', 'mode 4')
% %     
% %     name=sprintf('WG_GVD_thick_%dm_%dnm_%d_%s',  lambda(s(1)+2), w_array(s(3)), floor(rh_por_array(s(5))*100), timestamp);
% %     saveas(gcf,name, 'png')
% % 
% %     
% % % % %     Plotting 1 variable at the time, keeping other at centre
% % % % INDEX
% % figure(500)
% % clf;
% % for ii=1:nmodes
% %     d=size(All_eff_index);
% %     s=floor(1/2*d);
% %     plot_val=[];
% %     
% %      s(1)=s(1)+floor((1550-lambda(s(1)))/(lambda(2)-lambda(1))) -1;            %wavelength
% %     s(3)=d(3)-1;        %width
% %     s(4)=d(4)-2;      %film thick
% % %     
% %     for i_index=1:d(5)  %etched depth       
% %         j_index=(i_index-1)*d(4)*d(3)*d(2)*d(1) + (s(4)-1)*d(3)*d(2)*d(1)...
% %                 +(s(3)-1)*d(2)*d(1) + (ii-1)*d(1) + s(1)+2;
% %             
% % % %         j_index=(s(5)-1)*d(4)*d(3)*d(2)*d(1)+(s(4)-1)*d(3)*d(2)*d(1)...
% % % %                 +(s(3)-1)*d(2)*d(1) + (s(2)-1)*d(1) + s(1)+2
% % 
% % % %         j_index=d2(1)+2+(ii-1)*d(1)+(d2(3)-1)*d(2)*d(1)+(i_index-1)*d(3)*d(2)*d(1);
% % 
% %         plot_val(i_index)=All_eff_index(j_index);
% %         
% %     end;
% %     figure(500)
% %     hold all
% %     plot(rh_por_array*100, plot_val,'*-')
% % 
% % end
% %     figure(500)
% %     xlabel('Etched depth (%)');
% %     ylabel('Index');
% %     title1=sprintf('Index %dnm wl %dnm wide, %dnm thick', lambda(s(1)+2), w_array(s(3)), h2_array(s(4)));
% %     title(title1);
% %     legend('mode 1', 'mode 2', 'mode 3', 'mode 4')
% %     name=sprintf('WG_index_ridge_%dnm_%dnm_%dnm_ridge_%s', lambda(s(1)+2), w_array(s(3)), h2_array(s(4)), timestamp);
% %     saveas(gcf,name, 'png')
% %     
% % 
% % figure(600)
% % clf;
% % for ii=1:nmodes
% %     d=size(All_eff_index);
% %     s=floor(1/2*d);
% %     plot_val=[];
% %         
% %      s(1)=s(1)+floor((1550-lambda(s(1)))/(lambda(2)-lambda(1)))-1;            %wavelength
% %     s(3)=d(3)-1;      %width
% % %     s(5)=3;         %etched depth
% % 
% % 
% %     for i_index=1:d(4)    %thick
% %         j_index=(s(5)+1-1)*d(4)*d(3)*d(2)*d(1)+(i_index-1)*d(3)*d(2)*d(1)...
% %                 +(s(3)-1)*d(2)*d(1) + (ii-1)*d(1) + s(1)+2;
% % % %         j_index=d2(1)+2+(ii-1)*d(1)+(i_index-1)*d(2)*d(1)+(d2(4)+1-1)*d(3)*d(2)*d(1);
% %         plot_val(i_index)=All_eff_index(j_index);
% %     end;
% %     figure(600)
% %     hold all
% %     plot(h2_array, plot_val,'*-')
% % 
% % end
% %     figure(600)
% %     xlabel('Film thick (nm)');
% %     ylabel('Index');
% %     title1=sprintf('Index %dnm wl, %dnm wide, %d percent etched', lambda(s(1)+2), w_array(s(3)), floor(rh_por_array(s(5))*100));
% %     title(title1);
% %     legend('mode 1', 'mode 2', 'mode 3', 'mode 4')
% %     
% %     name=sprintf('WG_index_thick_%dm_%dnm_%dnm_GVD_%s',  lambda(s(1)+2), w_array(s(3)), floor(rh_por_array(s(5))*100), timestamp);
% %     saveas(gcf,name, 'png')
    
    
endtime=datestr(now)