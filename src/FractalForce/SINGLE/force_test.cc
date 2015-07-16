#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  template <class M>  void force_test(M& mem,Fractal& fractal)
  {
    static ofstream FileForce;
    if(!FileForce.is_open())
      FileForce.open("force.d");
    FileForce.precision(5);
    double sum_m=0.0;
    double sum_fx1=0.0;
    double sum_fy1=0.0;
    double sum_fz1=0.0;
    double sum_fx2=0.0;
    double sum_fy2=0.0;
    double sum_fz2=0.0;
    vector <double> force(3,0.0);
    double mass=0.0;
    for(int n=0;n < fractal.get_number_particles();++n)
      {
	Particle* p=fractal.particle_list[n];
	mass=p->get_mass();
	p->get_force(force);
	sum_m+=mass;
	sum_fx1+=mass*force[0];
	sum_fy1+=mass*force[1];
	sum_fz1+=mass*force[2];
	sum_fx2+=pow(mass*force[0],2);
	sum_fy2+=pow(mass*force[1],2);
	sum_fz2+=pow(mass*force[2],2);
      }
    FileForce << fractal.get_steps() << " " << sum_m << " ";
    FileForce << sum_fx1 << " " << sqrt(sum_fx2) << " ";
    FileForce << sum_fy1 << " " << sqrt(sum_fy2) << " ";
    FileForce << sum_fz1 << " " << sqrt(sum_fz2) << endl;
  }
}
namespace FractalSpace
{
template void force_test(Fractal_Memory& mem, Fractal& fractal);
}
