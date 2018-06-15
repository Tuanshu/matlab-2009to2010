         clear all;
		 %���W�{��
   w=0;
         
%===========���ܰѼ�===============
          Material='LT';              % LT or LN
          Scanning_method='wavelength';         %wavelength or power
          Scanning_method_2='ratio';     %ratio or random pitch
          QPM_wavelength=1.55;   
          Fundamental_wavelength=1.53:0.01:1.57;        %micron
          Ratio=1;
          Averaged_power=0.12;                      %�����J�g���\�v(W)
          Repetition_rate=4000;         %Hz
          Pulsewidth=27E-9;             %second
          Averaging_times=1;
%=================================

%===========���`���ʪ��Ѽ�=========
            q_p=6000;                %�N�C��QPM pitch�����o��h��pitch (spitially) (QPM pitch/q_p=dz)
            random_pitch=1;        %���dz, �̦n�]��q_p/2 (coherent length)���]��, ���M�|�J�촫poling��V�ɤ@��pitch�Q���_�����p(���z�Ϲ�����)
            number_of_period=250;      %�g���ƶq, ���ƭȭ��WQPM�g��(�o�O��X�Ӫ�)�Ydevice����
%====================================

if strcmp(Scanning_method,'wavelength')
    Scanning_parameter=Fundamental_wavelength';
else if strcmp(Scanning_method,'power')
        Scanning_parameter=Averaged_power';
    end
end
if strcmp(Scanning_method_2,'ratio')
    Scanning_parameter_2=Ratio';
else if strcmp(Scanning_method_2,'random_pitch')
        Scanning_parameter_2=random_pitch';
    end
end


%===========�w���ŧi�}�C����==========
Z01(1:length(Scanning_parameter))=0;
W3(1:length(Scanning_parameter))=0;
Z03(1:length(Scanning_parameter))=0;
area(1:length(Scanning_parameter))=0;
SHG_power_temp(1:length(Scanning_parameter))=0;
SHG_power(1:length(Scanning_parameter),1:length(Scanning_parameter_2))=0;
SHG_wavelength_t(1:length(Scanning_parameter))=0;
Averaging_temp(1:length(Scanning_parameter))=0;
A(1:q_p)=0;
B(1:q_p)=0;
delta_k_SHG(1:q_p)=0;
position(1:q_p)=0;
M_SHG(1:q_p)=0;
power_run_A(1:q_p)=0;
power_run_B(1:q_p)=0;
%SHG_power_temp(1:Averaging_times,1:(Scanning_range/Scanning_pitch)+1)=0;
%==================================== 
         

		 W01=0.07;            %��J����e�������y�b�|(mm)
		 fL=50;                %�z��J�Z(mm)
		 n=1;                   %�u�Ť���g�v
		 C=3E8;                 %���t(m/s)
		 pi=3.1415926;          %��P�v
		 d=0.6E-6;              %���P���첾�q(m)

		 sus=8.854E-12;         %�u�Ť����q�`��(F/m)   (C^2/N/m^2)
		 perm=pi*4E-7;          %�u�Ť��ɺϫY��(H/m)   (N/A^2)
         
         T=25;                  %�ū�

if strcmp(Material,'LN')
         
		    fe=(T-24.5)*(T+570.82);                     %Sellimeier equation �Ѽ� for LN
	        c1=5.35583;                                
		    c2=0.100473;
		    c3=0.20692;
		    c4=100;
		    c5=11.34927;
		    c6=-1.5334e-2;
		    d1=4.629e-7;
		    d2=3.826e-8;
		    d3=-8.9e-9;
		    d4=2.657e-5	;	 
%=============================�H�U�OLT��
else if strcmp(Material,'LT')
	        a=4.514261                      ;           
		    b=0.011901;
		    c=0.110744;
            d=-0.02323;
            e=0.076144;
            f=0.195596;
            bT=(1.82194*1E-8)*(T+273.15)^2;
            cT=(1.5662*1E-8)*(T+273.15)^2;
    end
end          
              
        
			
	      %power_p=q*10


		  %do r=27,30,1
          %write(*,*),r


		  %vi=0.1
          %do rv_round=1,14,1

	      %call random_seed()
		  %do rv=1,301,1
		  %call random_number(r1)
		  %r1_array(rv)=1-vi+2*vi*r1
		  %end do

          %=======================================!
		  %     �j��-1 (�̳̳̥~��) ==> ���P��Ratio or random pitch
		  %=======================================!           

          
for y=1:length(Scanning_parameter_2)         
          
          
          
          %=======================================!
		  %     �j��0 (�̳̥~��) ==> �h������SHG�W��
		  %=======================================! 
          
for z=1:Averaging_times
          
          
          %=======================================!
		  %     �j��1 (�̥~��) ==> �]�i�� or power�ܤ�
		  %=======================================! 

for j=1:length(Scanning_parameter)                                
        %write(*,*),j
% lam_p=QPM_wavelength;
        
if strcmp(Scanning_method,'wavelength')
            lam_p=Fundamental_wavelength(j);   %���W���i��(um) �Oj�����
            power_p=Averaged_power/(Repetition_rate*Pulsewidth);    %��peak power��
else if strcmp(Scanning_method,'power')
            lam_p=Fundamental_wavelength;
            power_p=Averaged_power(j)/(Repetition_rate*Pulsewidth);    %��peak power�� �Oj�����
    end
end
		    lam_c=lam_p/2;             %���W���i��(um)
		
            
			%write(1,*) lam_p*1E3
		
			 
			fre_p=2*pi*C/lam_p*1E6;   %���W���W�v(Hz)    �o�̪�fre���O�����W�v
		    
		    fre_c=2*pi*C/lam_c*1E6;   %���W���W�v(Hz)

            
      if strcmp(Material,'LN')
			np=(c1+d1*fe+(c2+d2*fe)/(lam_p^2-(c3+d3*fe)^2)+(c4+d4*fe)/(lam_p^2-(c5)^2)+c6*lam_p^2)^0.5;  %���W����g�v FOR LN  		    
		    nc=(c1+d1*fe+(c2+d2*fe)/(lam_c^2-(c3+d3*fe)^2)+(c4+d4*fe)/(lam_c^2-(c5)^2)+c6*lam_c^2)^0.5;  %���W����g�v FOR LN
            
      else if strcmp(Material,'LT')            
            np=(a+(b+bT)/(lam_p^2-(c+cT)^2)+e/(lam_p^2-f^2)+d*lam_p^2)^0.5 ;
		    nc=(a+(b+bT)/(lam_c^2-(c+cT)^2)+e/(lam_c^2-f^2)+d*lam_c^2)^0.5 ;
          end
      end
		    %write(*,*) np
			%write(*,*) nc

			Z01(j)=pi*(W01)^2*n/(lam_p*1E-3);
		
		    W3(j)=(fL/Z01(j))/((1+(fL/Z01(j))^2)^0.5)*W01; %�J�g�������y�b�|

		    Z03(j)=pi*(W3(j))^2*n/(lam_p*1E-3);             %rayleigh range

			area(j)=pi*(W3(j)*1E-3)^2;                       %���I���n(m^2)
		 
          A(1)=((2*power_p*(perm/sus)^(0.5))/(area(j)*fre_p))^(0.5);     %���W���q����l��
		  B(1)=0;                                                        %���W���q����l��
		  
          %�U�����Odesigned pitch, �ҥH�Mj�L��
		   %lamp_gpr=1.550                  %�p��_�l�g�����׹��������W���i��
		    lamp_gpr=QPM_wavelength;
		
            lamc_gpr=lamp_gpr/2;              %�p��_�l�g�����׹��������W���i��
		    
			
			%write(*,*) lamc_gpr
            
if strcmp(Material,'LN')
			np_gpr=(c1+d1*fe+(c2+d2*fe)/(lamp_gpr^2-(c3+d3*fe)^2)+(c4+d4*fe)/(lamp_gpr^2-(c5)^2)+c6*lamp_gpr^2)^0.5;  %�p��_�l�g�����׹��������W����g�v	FOR LN	    
		    nc_gpr=(c1+d1*fe+(c2+d2*fe)/(lamc_gpr^2-(c3+d3*fe)^2)+(c4+d4*fe)/(lamc_gpr^2-(c5)^2)+c6*lamc_gpr^2)^0.5;  %�p��_�l�g�����׹��������W����g�v
else if strcmp(Material,'LT')     		                        
            np_gpr=(a+(b+bT)/(lamp_gpr^2-(c+cT)^2)+e/(lamp_gpr^2-f^2)+d*lamp_gpr^2)^0.5 ;
		    nc_gpr=(a+(b+bT)/(lamc_gpr^2-(c+cT)^2)+e/(lamc_gpr^2-f^2)+d*lamc_gpr^2)^0.5 ;
    end
end
           	%gpr_1=(nc_gpr/lamc_gpr-np_gpr/lamp_gpr-np_gpr/lamp_gpr)^(-1) ;   	
%gpr_1=Scanning_pitch*(j-1);
        
		    %write(*,*) np_gpr,nc_gpr,nu_gpr
			
			gpr_1=(nc_gpr/lamc_gpr-np_gpr/lamp_gpr-np_gpr/lamp_gpr)^(-1);	%�@�����W������QPM�g��(um) 
			
            
		    dz=gpr_1/q_p*1E-6;              %�p�ⴹ����פ����Z(m)        
			%dz=(gpr_2+d*s+r*h)/600*1E-6 	%�p�ⴹ����פ����Z(m)		
                                        



		    %do d_change=1,10000,1
			%deff_1=0.0001*(d_change)*1.939E-22 
if strcmp(Material,'LN')
            d_BPM=sus*34.4E-12;                                                          %�����D�o�Oԣ��� to check FOR LN
            else if strcmp(Material,'LT') 
            d_BPM=sus*8.5E-12;              %    8.5E-12 m/V, �Ѧ�Hao 071220           
                end
end
            deff_1_p=d_BPM;                                                         %�D�u�ʫY��(AS/V^2) =d_BPM*2/(pi)
 %deff_1_p=0.8962E-22;
            deff_1_n=-deff_1_p;
            %write(*,*) 2*((perm/sus)**1.5)*((c/lam_p*1E6)**2)*(d_BPM**2)*(0.0032**2)/(np**3)*power_p/(area(j))*100
             %==============================================


		  %================================================================!
		  %     �j��2 (��2��) ==> ���w�g�����׻P�Ӽ�
		  %                       
		  %================================================================! 

        for h=1:number_of_period  %�g�����ƶq
		  
		  %===========================================================================!
		  %     �j��3 (��3��) ==> ���W�����ӻP���W��������                     
		  %===========================================================================! 
		
    for q=1:q_p                        %q��position��index   600���@�g��
             
if strcmp(Scanning_method_2,'ratio')
    
              if q<=q_p/2              %�g������,duty cycle=50%       �o��i�H�y�L�g����²�@�I, ���o�ˤ���M��
                    if (mod(q-1,random_pitch)==0)       
                        if rand<Ratio(y)
                            deff_1=deff_1_n;
                        else
                            deff_1=deff_1_p;
                        end
                    end
                else
                    if mod(q-1,random_pitch)==0
                        if rand<Ratio(y)
                            deff_1=deff_1_p;
                        else
                            deff_1=deff_1_n;
                        end
                    end                             %�S��check�ɤ��򳣤���, �]�N�O����
               end
 
    
else if strcmp(Scanning_method_2,'random_pitch')
        
           if q<=q_p/2              %�g������,duty cycle=50%       �o��i�H�y�L�g����²�@�I, ���o�ˤ���M��
                    if (mod(q-1,random_pitch(y))==0)       
                        if rand<Ratio
                            deff_1=deff_1_n;
                        else
                            deff_1=deff_1_p;
                        end
                    end
                else
                    if mod(q-1,random_pitch(y))==0
                        if rand<Ratio
                            deff_1=deff_1_p;
                        else
                            deff_1=deff_1_n;
                        end
                    end                             %�S��check�ɤ��򳣤���, �]�N�O����
            end
 
    end
end

		
  
                
			%write(*,*) i,deff_1
		   
			%write(*,*)np,nc

            sus_p=sus*np^2;                  %���W�������q�Y��
		   
		    sus_c=sus*nc^2;                  %���W�������q�Y��
            
			%write(*,*)  sus_p,sus_c

            lo_p=0;                           %���W���ǿ�l�� 
           
            lo_c=0;                           %���W���ǿ�l�� 
         

			abs_p=lo_p*(perm/sus_p)^(0.5);   %���W���l���Y��
            
            abs_c=lo_c*(perm/sus_c)^(0.5);   %���W���l���Y��
            
			%write(*,*) abs_p,abs_c
	
            delta_k_SHG(q)=2*pi*(nc/(lam_c*1E-6)-np/(lam_p*1E-6)-np/(lam_p*1E-6));   %���W�ۦ줣�ǰt�{�� �o�Ӧn�����άOarray
			
			%delta_k_1=0                                                    
			%write(*,*) delta_k_1

			position(q+1)=position(q)+dz;                   %�����ǻ���m
           
		    M_SHG(q)=i*delta_k_SHG(q)*position(q);
           

			k_1=deff_1*((perm/sus)*fre_p*fre_p*fre_c/(np*np*nc))^(0.5);                %���W�D�u�ʽ��X�Y��
           
           A(q+1)=A(q)+dz*(-0.5*abs_p*A(q)-0.5*i*k_1*B(q)*conj(A(q))*exp(-M_SHG(q)));  %�p����W��������
		   B(q+1)=B(q)+dz*(-0.5*abs_c*B(q)-0.5*i*k_1*A(q)^2*exp(M_SHG(q)));            %�p�⭿�W��������
          
		   power_run_A(q)=(0.5*(sus/perm)^(0.5))*fre_p*((abs(A(q)))^2)*area(j);     %�g�L�@�q���שҳѤU�����W���\�v(W) �Ǫ����{������2, �Q�ڮ����F
           power_run_B(q)=(0.5*(sus/perm)^(0.5))*fre_c*((abs(B(q)))^2)*area(j);       %�g�L�@�q���שҲ��ͪ����W���\�v(W)
                      
        
          % total_A((h-1)*q_p+(q+1))=power_run_A(q);
          % total_B((h-1)*q_p+(q+1))=power_run_B(q);
	
           
		 
    end                                %�j��3���� q
           
		   position(1)=position(q_p+1);           %���������
		   A(1)=A(q_p+1);                         %���W������ǻ�
           B(1)=B(q_p+1);                         %���W������ǻ�
	                                
           
        end                               %�j��2���� h

  

		% SHG_power_temp(z,j)=power_run_B(q_p)';                % (W)
         SHG_power_temp(j)=power_run_B(q_p)'*(Repetition_rate*Pulsewidth);                % (W), �b�o�̶}�l�Ҽ{�����\�v
         
		  position(1)=0;                        %���ܪi���e���N�����m�k�s           
         
         
         %write(6,*),real(power_run_A(600))+real(power_run_B(600))
          %write(*,*) real(power_run_B(600))
          %write(*,*) real(power_run_B(600))
	 	 
          %if( real(power_run_B(600)) > 120.098 ) then
     	  %write(*,*) "found"
	      %write(3,*) deff_1
	      %end if

	     %end do
		  
end                                    %�j��1���� j

%Averaging_temp=Averaging_temp+SHG_power_temp(z,:);
Averaging_temp=Averaging_temp+SHG_power_temp;
end

SHG_power(:,y)=Averaging_temp/Averaging_times;
Averaging_temp=0;

end
SHG_power_nanoW=SHG_power*1E9;
 plot(Scanning_parameter,SHG_power);
 %plot(power_run_B);		  
