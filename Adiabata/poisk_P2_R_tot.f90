 program poisk_P2_R_tot
 implicit none 
 integer :: i
 integer, PARAMETER :: M=1000000
Double Precision, PARAMETER :: two_pii = 6.2831853071796, pii = 3.1415926536, c_cgs=2.99792458D+10, e_cgs=4.8032068D-10,&
                    four_pii= 12.566370614359D0,  &
                    eight_pii= 25.132741228718, & 
                    mp_cgs=1.6726231D-24, sp_ht_rel=5.0D0/3.0D0, &
                    kb_cgs=1.380658D-16
Double Precision, PARAMETER :: c2_cgs       = c_cgs**2.0D0   ! Speed of light squared, cm^2/s^2
Double Precision, PARAMETER :: mpc2_cgs     = mp_cgs*c_cgs**2.0D0    ! Rest mass of proton, ergs 
Double Precision, PARAMETER :: mpc_cgs      = mp_cgs*c_cgs 
Double Precision, PARAMETER :: Gam0=5.0D0/3.0D0 
Double Precision, PARAMETER :: delta_W=2.0D0, alpha_W=delta_W/(delta_W-1.0D0)  
Double Precision, PARAMETER :: InterstellarMagneticField=8.0D-3, &
    UpstreamNumberDensity=3.0D0, &
    gamma_sh=1.5D0, &
    !gamma_sh=2.0D0, &
    UpstreamFlowSpeed=c_cgs*(1.0D0-1.0D0/gamma_sh**2.0D0)**0.5D0, &
    UpstreamTemperature=1.0D+6
double precision, parameter :: P_w_Up=InterstellarMagneticField**2.0D0/eight_pii    
double precision, parameter :: P_th_Up=UpstreamNumberDensity*kb_cgs*UpstreamTemperature
  Double Precision :: Rtot_new, P2_Pxx, P2_En, T_2
  Double Precision :: Rtot_new_a, Rtot_new_b
  Double Precision :: betta0
  Double Precision :: Gam2,Qesc_Pxx, Qesc_En
  Double Precision :: P_w0, P_w2
  Double Precision :: P_th0
  Double Precision :: PP, EE
  Double Precision :: x1, x2
  Double Precision, PARAMETER :: x_m=100, x_h=x_m/M
  
  
  betta0=UpstreamFlowSpeed/c_cgs
  
  P_w0=P_w_Up/(mpc2_cgs*UpstreamNumberDensity)
  P_th0=P_th_Up/(mpc2_cgs*UpstreamNumberDensity)
  
  
  Gam2=1.56548D0
 ! Gam2=1.5147D0
  write(*,*) 'Gam2 0=', Gam2 
  Qesc_Pxx=0.0D0
  Qesc_En=0.0D0
        
  
   Rtot_new_a=2.0D0 
   Rtot_new_b=8.0D0     
   Rtot_new=(Rtot_new_a+Rtot_new_b)/2.0D0
   DO i=1,100
      P_w2=P_w0*gamma_sh**alpha_W*(Rtot_new**2.0D0-betta0**2.0D0)**(3.0D0/4.0D0)
      P2_Pxx=(1.0D0/(Gam2/((Gam2-1.0D0)*gamma_sh**2.0D0*(Rtot_new**2.0D0-betta0**2.0D0))+ &
              1.0D0/(gamma_sh**2.0D0*betta0**2.0D0)))* &
              (1.0D0-1.0D0/(gamma_sh*DSQRT(Rtot_new**2.0D0-betta0**2.0D0))+Qesc_Pxx+ &
              P_th0*(Gam0/(Gam0-1.0D0)+1.0D0/(gamma_sh**2.0D0*betta0**2.0D0))+ &
              P_w0*(delta_W+1.0D0/(gamma_sh**2.0D0*betta0**2.0D0))- &
              P_w2*(delta_W/(gamma_sh**2.0D0*(Rtot_new**2.0D0-betta0**2.0D0)) &
					+1.0D0/(gamma_sh**2.0D0*betta0**2.0D0))) 
      P2_En=((Gam2-1.0D0)*gamma_sh**2.0D0*(Rtot_new**2.0D0-betta0**2.0D0)/(Rtot_new*Gam2))* &
              (1.0D0-Rtot_new/(gamma_sh*DSQRT(Rtot_new**2.0D0-betta0**2.0D0))+Qesc_En+ &
              P_th0*Gam0/(Gam0-1.0D0)+ &
              P_w0*delta_W-P_w2*delta_W*Rtot_new/(gamma_sh**2.0D0*(Rtot_new**2.0D0-betta0**2.0D0)))
              
       IF ((P2_En-P2_Pxx).GT.0.0D0) then   
            Rtot_new_b=Rtot_new
            Rtot_new=(Rtot_new_a+Rtot_new_b)/2.0D0       
       else IF ((P2_En-P2_Pxx).LT.0.0D0) then 
            Rtot_new_a=Rtot_new
            Rtot_new=(Rtot_new_a+Rtot_new_b)/2.0D0 
       else
            exit
       end IF  
  
   end DO  
   
   T_2=P2_Pxx/(gamma_sh*(Rtot_new**2.0D0-betta0**2.0D0)**0.5D0)   
  
  PP=0.0D0
  EE=0.0D0
  DO i=1, M 
    x2=x_h*i
    x1=x_h*(i-1)
    PP=PP+0.5D0*x_h*(x1**2.0D0*Dexp(-(1.0D0+x1**2.0D0)**0.5D0/T_2) + &
        x2**2.0D0*Dexp(-(1.0D0+x2**2.0D0)**0.5D0/T_2)) 
    EE=EE+0.5D0*x_h*(((1.0D0+x1**2.0D0)**0.5D0-1.0D0)*x1**2.0D0*Dexp(-(1.0D0+x1**2.0D0)**0.5D0/T_2) + &
        ((1.0D0+x2**2.0D0)**0.5D0-1.0D0)*x2**2.0D0*Dexp(-(1.0D0+x2**2.0D0)**0.5D0/T_2))     
  end DO 
  
  Gam2=1.0D0+T_2*PP/EE
 
 write(*,*) 'P2_Pxx=', P2_Pxx
 write(*,*) 'P2_En=', P2_En
 write(*,*) 'T_2', T_2
 write(*,*) 'Rtot_new=', Rtot_new
 write(*,*) EE, T_2*PP
 write(*,*) 'Gam2=', Gam2
 
! PAUSE

 end program poisk_P2_R_tot

