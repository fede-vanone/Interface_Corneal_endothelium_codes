function [Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,var)
[a c k NR]=size(Y);
%Y=Y(:,time_points,:,:,:);
if nargin<4
    disp(['ERROR = Choose one or more outputs from var 1 and variable ',num2str(c),' of the model'])
    error('eFAST: the output components for the sensitivity is missing. Not enough input arguments.')
end
 %%% done by using trapz
for u=var
        ttest_sti=zeros(k,k);
        ttest_si=ttest_sti;
        for i=1:k %loop through parameters
            % Initialize AV,AVi,AVci to zero.
            AV = 0;
            AVi = 0;
            AVci = 0;
            for L=1:NR%length(Y(1,u,i,:))
                ind=isnan(Y(:,u,i,L));
                rr=1:length(Y(:,u,i,L));
                Yn=interp_nonan(Y,u,i,L,77);
                shift=50;
                while isnan(Yn(:,u,i,L))~=zeros(size(Yn(:,u,i,L)))
                    shift=shift+10;
                    Yn=interp_nonan(Y,u,i,L,150+shift);
                end

                    
%                 figure;
%                 plot(rr,Yn(:,u,i,L))
%                 hold on 
%                 plot(rr,Y(:,u,i,L),'o')
                %Y(ind,u,i,L)=0;
                Z(:,u,i,L) = (Yn(:,u,i,L)-mean(Y(~ind,u,i,L)))';
                % Fourier coeff. at [1:OMi/2].
                N=length(Yn(:,u,i,L));
                NQ = (N-1)/2;
                N0 = NQ+1;
                COMPL = 0;
                Y_k = [Z(:,u,i,L);Z(1,u,i,L)];
                F=trapz(pi*(2*(0:N)-N)/N,Y_k)/2/pi;
                for j=1:OMi/2
                    AC(j)=0; BC(j)=0;
                    ANGLE = pi*(2*(0:N)-N)/N;
                    C_VEC = cos(ANGLE*j);
                    S_VEC = sin(ANGLE*j);
                    
                        AC(j) = trapz(ANGLE,Y_k.*C_VEC')/2/pi;
                        BC(j) = trapz(ANGLE,Y_k.*S_VEC')/2/pi;
                    F=F+AC(j)*cos(j*ANGLE)+BC(j)*sin(j*ANGLE);
                    COMPL = COMPL+AC(j)^2+BC(j)^2;
                end
                Vci(L) = 2*COMPL;
%                 for j=OMi/2+1:OMi*MI*2
%                     ANGLE = pi*(2*(0:N)-N)/N;
%                     C_VEC = cos(ANGLE*j);
%                     S_VEC = sin(ANGLE*j);
%                     
%                         AC(j) = trapz(ANGLE,Y_k.*C_VEC')/2/pi;
%                         BC(j) = trapz(ANGLE,Y_k.*S_VEC')/2/pi;
%                     F=F+AC(j)*cos(j*ANGLE)+BC(j)*sin(j*ANGLE);
%                 end
%                 figure
%                 plot(ANGLE,F)
%                 hold on 
%                 plot(ANGLE,Y_k,'--')
                % Computation of V_{(ci)}.
                Vci(L) = 2*COMPL;
                %AVci = AVci+Vci;
                % Fourier coeff. at [P*OMi, for P=1:MI].
                COMPL = 0;
                for j=OMi:OMi:OMi*MI 
                    AC(j)=0; BC(j)=0;
                    ANGLE = pi*(2*(0:N)-N)/N;
                    C_VEC = cos(ANGLE*j);
                    S_VEC = sin(ANGLE*j);
                    
                        AC(j) = trapz(ANGLE,Y_k.*C_VEC')/2/pi;
                        BC(j) = trapz(ANGLE,Y_k.*S_VEC')/2/pi;
                    
                    COMPL = COMPL+AC(j)^2+BC(j)^2;
                end% Computation of V_i.
                Vi(L) = 2*COMPL;
                %AVi = AVi+Vi;
                % Computation of the total variance
                % in the time domain.
                V(L) = Z(:,u,i,L)'*Z(:,u,i,L)/length(Z(:,u,i,L));
                clear Z Yn
            end %L
            
            % Computation of sensitivity indexes.
            %AV = AV/length(Y(1,t,1,i,:)); %AV=average total variance
            %AVi = AVi/length(Y(1,t,1,i,:));
            %AVci = AVci/length(Y(1,t,1,i,:));
            Si(i,u) = mean(Vi)/mean(V);
            Sti(i,u) = 1-mean(Vci)/mean(V);
            rangeSi(i,:,u) = Vi./V;
            rangeSti(i,:,u) = 1-(Vci./V);
        end %i
end