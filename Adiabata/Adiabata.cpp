// Adiabata.cpp : Defines the entry point for the console application.
//
#include "stdio.h"
#include <stdlib.h>
#include "math.h"

double sqr(const double& a){
	return a*a;
}

int main()
{

	int  i;
	const int M=1000000;
	const double two_pii = 6.2831853071796;
	const double pii = 3.1415926536;
	const double c_cgs=2.99792458E+10;
	const double e_cgs=4.8032068E-10;
	const double four_pii= 12.566370614359;
	const double eight_pii= 25.132741228718; 
	const double mp_cgs=1.6726231E-24;
	const double sp_ht_rel=5.0/3.0;
	const double kb_cgs=1.380658E-16;
	const double c2_cgs       = c_cgs*c_cgs;   //! Speed of light squared, cm^2/s^2
	const double mpc2_cgs     = mp_cgs*c2_cgs; //    ! Rest mass of proton, ergs 
	const double mpc_cgs      = mp_cgs*c_cgs; 
	const double Gam0=5.0/3.0; 
	const double delta_W=2.0;
	const double alpha_W=delta_W/(delta_W-1.0);
	//const double InterstellarMagneticField=8.0E-3;
	const double InterstellarMagneticField=0.0064;
	//const double gamma_sh=1.5;
	const double gamma0 = 10.0;
	const double beta0 = sqrt(1.0 - 1.0/sqr(gamma0));
	const double beta_sh = 0.504;
	const double beta1 = (beta0 + beta_sh)/(1 + beta0*beta_sh);
	const double gamma_sh= 1.0/sqrt(1.0 - beta1*beta1);
	const double UpstreamNumberDensity=1.0/gamma_sh;
	const double UpstreamFlowSpeed=c_cgs*sqrt((1.0-sqr(1.0/gamma_sh)));
	const double UpstreamTemperature=1.0E-4*mpc2_cgs/kb_cgs;
	const double P_w_Up=sqr(InterstellarMagneticField)/eight_pii ;   
	const double P_th_Up=UpstreamNumberDensity*kb_cgs*UpstreamTemperature;
	double Rtot_new, P2_Pxx, P2_En, T_2;
	double Rtot_new_a, Rtot_new_b;
	double betta0;
	double Gam2,Qesc_Pxx, Qesc_En;
	double P_w0, P_w2;
	double P_th0;
	double PP, EE;
	double x1, x2;
	double x_m=100; 
	double x_h=x_m/M;
  
  
	betta0=UpstreamFlowSpeed/c_cgs;
  
	P_w0=P_w_Up/(mpc2_cgs*UpstreamNumberDensity);
	P_th0=P_th_Up/(mpc2_cgs*UpstreamNumberDensity);
  
    double Gam_a = 1.2;
	double Gam_b = 1.8;
	Gam2=1.56548;
	//printf("Gam2 0= %g \n", Gam2); 
	for(int j = 0; j < 10; ++j){
		double Gam_new = (Gam_a + Gam_b)/2;
		Gam2 = Gam_new;

		Qesc_Pxx=0.0;
		Qesc_En=0.0;
        
  
		Rtot_new_a=1.01;
		Rtot_new_b=8.0;     
		Rtot_new=(Rtot_new_a+Rtot_new_b)/2.0;
		for (i=1; i <= 100; ++i){
			P_w2=P_w0*pow(gamma_sh,alpha_W)*pow(sqr(Rtot_new)-sqr(betta0),1.0);
			P2_Pxx=(1.0/(Gam2/((Gam2-1.0)*sqr(gamma_sh)*(sqr(Rtot_new)-sqr(betta0)))+ 
              1.0/(sqr(gamma_sh)*sqr(betta0))))* 
              (1.0-1.0/(gamma_sh*sqrt(sqr(Rtot_new)-sqr(betta0)))+Qesc_Pxx+ 
              P_th0*(Gam0/(Gam0-1.0)+1.0/(sqr(gamma_sh)*sqr(betta0)))+ 
              P_w0*(delta_W+1.0/(sqr(gamma_sh)*sqr(betta0)))- 
              P_w2*(delta_W/(sqr(gamma_sh)*(sqr(Rtot_new)-sqr(betta0))) 
					+1.0/(sqr(gamma_sh)*sqr(betta0)))) ;
			P2_En=((Gam2-1.0)*sqr(gamma_sh)*(sqr(Rtot_new)-sqr(betta0))/(Rtot_new*Gam2))* 
              (1.0-Rtot_new/(gamma_sh*sqrt(sqr(Rtot_new)-sqr(betta0)))+Qesc_En+ 
              P_th0*Gam0/(Gam0-1.0)+ 
              P_w0*delta_W-P_w2*delta_W*Rtot_new/(sqr(gamma_sh)*(sqr(Rtot_new)-sqr(betta0))));
              
			if ((P2_En-P2_Pxx) > 0.0) {   
				Rtot_new_b=Rtot_new;
				Rtot_new=(Rtot_new_a+Rtot_new_b)/2.0;       
			} else if ((P2_En-P2_Pxx) <0.0) {
				Rtot_new_a=Rtot_new;
				Rtot_new=(Rtot_new_a+Rtot_new_b)/2.0;
			} else{
				break;
			}
			  
		} 
   
		T_2=P2_Pxx/(gamma_sh*sqrt(sqr(Rtot_new)-sqr(betta0)))  ; 
  
		PP=0.0;
		EE=0.0;
		for (i=1; i <= M ; ++i){
			x2=x_h*i;
			x1=x_h*(i-1);
			PP=PP+0.5*x_h*(sqr(x1)*exp(-sqrt((1.0+sqr(x1)))/T_2) + 
				sqr(x2)*exp(-sqrt((1.0 + sqr(x2)))/T_2)) ;
			EE=EE+0.5*x_h*((sqrt(1.0+sqr(x1))-1.0)*sqr(x1)*exp(-sqrt((1.0+sqr(x1)))/T_2) + 
				(sqrt((1.0+sqr(x2)))-1.0)*sqr(x2)*exp(-sqrt((1.0+sqr(x2)))/T_2)) ;    
		} 
  
		Gam_new=1.0+T_2*PP/EE;
		if(Gam_new > Gam2){
			Gam_a = Gam2;
		} else {
			Gam_b = Gam2;
		}
	}
	Gam2 = (Gam_a + Gam_b)/2;
 
	printf("P2_Pxx= %g \n", P2_Pxx);
	printf("P2_En=%g \n", P2_En);
	printf("T_2 = %g \n", T_2);
	printf("T_2 in K= %g \n", T_2*mpc2_cgs/kb_cgs);
	printf("Rtot_new= %g \n", Rtot_new);
	printf("EE = %g T_2*PP = %g \n", EE, T_2*PP);
	printf("Gam2= %g \n", Gam2);
	getchar();

	return 0;
}

