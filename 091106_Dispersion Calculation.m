clear all;


%% ver. SHG-for NTU     &     wavelength

T=25;
fe_CLN=(T-24.5)*(T+570.82);                     %Sellimeier equation 參數 for LN
	        c1_CLN=5.35583;                                
		    c2_CLN=0.100473;
		    c3_CLN=0.20692;
		    c4_CLN=100;
		    c5_CLN=11.34927;
		    c6_CLN=-1.5334e-2;
		    d1_CLN=4.629e-7;
		    d2_CLN=3.826e-8;
		    d3_CLN=-8.9e-9;
		    d4_CLN=2.657e-5	;	 

            a_CLT=4.514261;           
		    b_CLT=0.011901;
		    c_CLT=0.110744;
            d_CLT=-0.02323;
            e_CLT=0.076144;
            f_CLT=0.195596;
            bT_CLT=(1.82194*1E-8)*(T+273.15)^2;
            cT_CLT=(1.5662*1E-8)*(T+273.15)^2;
            
            a_SLT=4.528254;           
		    b_SLT=0.012962;
		    c_SLT=0.242783;
            d_SLT=-0.02288;
            e_SLT=0.068131;
            f_SLT=0.177370;
            g_SLT=1.307470;
            h_SLT=7.061878;
            bT_SLT=(3.483933*1E-8)*(T+273.15)^2;
            cT_SLT=(1.607839*1E-8)*(T+273.15)^2;
            
            lam_f=0.56;
            
            ne_f_CLT=(a_CLT+(b_CLT+bT_CLT)./(lam_f.^2-(c_CLT+cT_CLT).^2)+e_CLT./(lam_f.^2-f_CLT^2)+d_CLT*lam_f.^2).^0.5 ;
            ne_f_SLT=(a_SLT+(b_SLT+bT_SLT)./(lam_f.^2-(c_SLT+cT_SLT).^2)+e_SLT./(lam_f.^2-f_SLT^2)+g_SLT./(lam_f.^2-h_SLT^2)+d_SLT*lam_f.^2).^0.5 ;
            ne_f_CLN=(c1_CLN+d1_CLN*fe_CLN+(c2_CLN+d2_CLN*fe_CLN)./(lam_f.^2-(c3_CLN+d3_CLN*fe_CLN)^2)+(c4_CLN+d4_CLN*fe_CLN)./(lam_f.^2-(c5_CLN)^2)+c6_CLN*lam_f.^2).^0.5;  %基頻光折射率 FOR LN 

                    
                    DS=diff(ne_f_CLN)/(0.0001);
                    DS=mean(DS);
                    COHERENT_LENGTH=0.098;
                    WALKING_LENGTH=24.7;
                    OCT_DISPERSION_BROADENED=((1+(COHERENT_LENGTH*WALKING_LENGTH*DS)^2)^0.5);