#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>

typedef double prec;
#include "Ensemble_param.h"

class wave_eq_evolve
{ 
	Ensemble_param Y; //funcion de x
	Ensemble_param G; //funcion de x
	std::vector<prec> eta_next; //funcion de x
	std::vector<prec> eta_current; //funcion de x 
	std::vector<prec> eta_prev; //funcion de x
	Ensemble_param eta_0;  //funcion de t

	prec dt;
	prec dx;
	int N_L;
	int N_T;

	boost::mt19937 &rng;

	void first_tic()
	{
		eta_current[0]=eta_0[0];
		for (int i = 1; i < eta_current.size()-1; ++i)
		{
			//eta_current[i]=m_eta_prev[i]+0.5*pow(m_DT/m_DX,2)*m_Y[i]*(m_eta_prev[i+1]-2*m_eta_prev[i]+m_eta_prev[i-1]);
		}
		eta_current[eta_current.size()-1]=0;
	}

	void evolve()
	{
		int idx_max=eta_next.size()-1;
		int idt_max=eta_0.size()-1;
		for (int i = 1; i < idt_max; ++i)
		{
			eta_next[0]=eta_0[i];
			for (int j = 1; j < idx_max; ++j)
			{
				eta_next[j]=(pow(dt/dx,2)*Y[j]*((1-0.5*pow(0.25*0.0229*(idx_max - j+1)*100/(idx_max),2))*(eta_current[j+1]-eta_current[j])+(1-0.5*pow(0.25*0.0229*(idx_max - j)*100/(idx_max),2))*(eta_current[j-1]-eta_current[j]))+dt*0.5*G[j]*eta_prev[j]+2*eta_current[j]-eta_prev[j])/(1+0.5*dt*G[j]);
			}
			eta_next[idx_max]=0;
			eta_prev=eta_current;
			eta_current=eta_next;
		}
	}

	public:

	    wave_eq_evolve(prec in_dt, prec in_dx, boost::mt19937 &in_rng, prec T, prec L) : N_L(L/in_dx), N_T(T/in_dt), Y(L/in_dx,"Y_i",0,0,in_rng,true,true), G(L/in_dx,"G_i",0,0,in_rng,true,true), eta_current(L/in_dx), eta_prev(L/in_dx), eta_next(L/in_dx), eta_0(T/in_dt,"eta_i",0,0,in_rng,true,true), dt(in_dt), dx(in_dx), rng(in_rng)
	    {

	    }

	    void run()
	    {
	    	first_tic();
	    	evolve();
	    }





};

int main()
{
	std::cout << "Hello World" << std::endl;

	return 0;
}