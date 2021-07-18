#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <sys/stat.h>
#include <omp.h>

typedef double prec;
#include "Ensemble_param.h"

class wave_eq_evolve
{ 
	Ensemble_param Y; //funcion de x
	Ensemble_param G; //funcion de x
	Ensemble_param Q; //funcion de x
	std::vector<prec> eta_np1; //funcion de x
	std::vector<prec> eta_n; //funcion de x 
	std::vector<prec> eta_nm1; //funcion de x
	Ensemble_param eta_0;  //funcion de t

	prec dt;
	prec dx;
	int N_L;
	int N_T;

	boost::mt19937 &rng;

	std::vector<std::vector<prec>> eta_save;
	void first_tic()
	{
		eta_n[0]=eta_0[0];
		for (int i = 1; i < eta_n.size()-1; ++i)
		{
			eta_n[i]=eta_nm1[i]+0.5*(1.0/Y[i])*pow(dt/dx,2)*(0.5*(Q[i]+Q[i+1])*(eta_nm1[i+1]-eta_nm1[i])-0.5*(Q[i]+Q[i-1])*(eta_nm1[i]-eta_nm1[i-1]));
		}
		eta_n[eta_n.size()-1]=0;
	}

	void evolve(prec t_save)
	{
		int idx_max=eta_np1.size()-1;
		int idt_max=eta_0.size()-1;
		int idt_save=t_save/dt;
		for (int i = 1; i < idt_max; ++i)
		{
	    	if((i+2)%10000==0)
	    	{
	    		std::cout << "time: " << (i+2)*dt << "/" << (idt_max+1)*dt <<'\r' << std::flush ;
	    	}
			eta_np1[0]=eta_0[i];
			#pragma omp parallel for
			for (int j = 1; j < idx_max; ++j)
			{
				eta_np1[j]=(1.0/(Y[j]+dt*0.5*G[j]))*(Y[j]*(2*eta_n[j]-eta_nm1[j])+pow(dt/dx,2)*(0.5*(Q[j]+Q[j+1])*(eta_n[j+1]-eta_n[j])-0.5*(Q[j]+Q[j-1])*(eta_n[j]-eta_n[j-1]))+0.5*G[j]*dt*eta_nm1[j]);
			}
			eta_np1[idx_max]=0;
			eta_nm1=eta_n;
			eta_n=eta_np1;
			if(i>idx_max-idt_save)
			{
				eta_save.push_back(eta_n);
			}
		}
		std::cout << std::endl;
	}

	void set_external_force()
	{
	    for (int i = 0; i < eta_0.size(); ++i)
	    {
	    	eta_0.assign(i,sin(i*dt));
	    }
	}

	public:

	    wave_eq_evolve(prec in_dt, prec in_dx, boost::mt19937 &in_rng, prec T, prec L, prec Y_med, prec G_med, prec Y_sigm=0, prec G_sigm=0) : N_L(L/in_dx), N_T(T/in_dt), Y(L/in_dx,"Y_i",Y_med,Y_sigm,in_rng,true,true), G(L/in_dx,"G_i",G_med,G_sigm,in_rng,true,true), Q(L/in_dx,"Q_i",1,0,in_rng,true,true), eta_n(L/in_dx), eta_nm1(L/in_dx,0), eta_np1(L/in_dx), eta_0(T/in_dt,"eta_i",0,0,in_rng,true,true), dt(in_dt), dx(in_dx), rng(in_rng)
	    {
	    	set_external_force();
	    }

	    void run(prec t_save=0)
	    {
	    	first_tic();
	    	evolve(t_save);
	    }

	    void save()
	    {
			mkdir("out", 0777);

			std::ofstream txtOut;
			txtOut.precision(std::numeric_limits< prec >::max_digits10);
			txtOut.open("out/X.txt");
			for (int j = 0; j < eta_n.size(); ++j)
			{
				txtOut << j*dx << " " << eta_n[j] << std::endl;
			}
			txtOut.close();
	    }
};

int main()
{
 	boost::mt19937 rng(static_cast<unsigned int>(std::time(0)));
	std::cout << "Hello World" << std::endl;

 	int nProcessors=omp_get_max_threads();
 	nProcessors=nProcessors-2;
    omp_set_num_threads(nProcessors);

    prec Y=1000;
    prec G=2500;
    prec dt=0.001;
    prec dx=0.001;
    prec T=1000;
    prec L=1;
	wave_eq_evolve a(dt,dx,rng,T,L,Y,G);
	a.run(100);
	a.save();

	return 0;
}