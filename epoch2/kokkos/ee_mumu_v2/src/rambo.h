#include <vector>
#include "Kokkos_Core.hpp"


#define PI 3.14159265358628
#define ACC 1e-14

template <typename ExecSpace>
Kokkos::View<double ***,ExecSpace> get_initial_momenta(
                                  Kokkos::View<double***,ExecSpace> d_p,
                                  const int nexternal,const double energy,
                                  Kokkos::View<double*,ExecSpace> masses,
                                  int events_per_iter){
  // printf("initialize momenta\n");
  // Set momenta for incoming particles
  Kokkos::parallel_for("get_initial_momenta",Kokkos::RangePolicy<ExecSpace>(0,events_per_iter),
    KOKKOS_LAMBDA(const int& i){

      auto e2 = pow(energy,2);
      double mom = sqrt((pow(e2, 2) - 2 * e2 * pow(masses(0), 2) + pow(masses(0), 4) -
          2 * e2 * pow(masses(0), 2) - 2 * pow(masses(0), 2) * pow(masses(1), 2) + pow(masses(1), 4)) /
         (4 * e2));
      auto energy1 = sqrt(pow(mom, 2) + pow(masses(0), 2));
      auto energy2 = sqrt(pow(mom, 2) + pow(masses(1), 2));

      for(int j=0;j<nexternal;++j)
        for(int k=0;k<4;++k)
          d_p(i,j,k) = 0.;
      // particle 1
      d_p(i,0,0) = energy1;
      d_p(i,0,3) = mom;
      // particle 2
      d_p(i,1,0) = energy2;
      d_p(i,1,3) = -mom;
    });
  return d_p;
}

template <typename ExecSpace>
void get_final_momenta(const int ninitial,const int nexternal,const double energy,
                                  Kokkos::View<double*,ExecSpace> masses,
                                  Kokkos::View<double ***,ExecSpace> d_p,
                                  Kokkos::View<double **,ExecSpace> random_numbers,
                                  Kokkos::View<double,ExecSpace> d_wgt, int dim) {
  //---- auxiliary function to change convention between MadGraph5_aMC@NLO and
  // rambo
  //---- four momenta.
  // int nfinal = nexternal - ninitial;
  
  // auto h_rns = Kokkos::create_mirror_view(random_numbers);
  // Kokkos::deep_copy(h_rns,random_numbers);

  Kokkos::View<int*,ExecSpace> iwarn("iwarn",5);
  
  // printf("initialize momenta\n");
  // Set momenta for incoming particles
  Kokkos::parallel_for("get_momenta_ninitial",Kokkos::RangePolicy<ExecSpace>(0,dim),
    KOKKOS_LAMBDA(const int& i){

      /**********************************************************************
       *                       rambo                                         *
       *    ra(ndom)  m(omenta)  b(eautifully)  o(rganized)                  *
       *                                                                     *
       *    a democratic multi-particle phase space generator                *
       *    authors:  s.d. ellis,  r. kleiss,  w.j. stirling                 *
       *    this is version 1.0 -  written by r. kleiss                      *
       *    -- adjusted by hans kuijf, weights are logarithmic (20-08-90)    *
       *                                                                     *
       *    n  = number of particles                                         *
       *    et = total centre-of-mass energy                                 *
       *    xm = particle masses ( dim=nexternal-nincoming )                 *
       *    p  = particle momenta ( dim=(4,nexternal-nincoming) )            *
       *    wt = weight of the event                                         *
       ***********************************************************************/
      constexpr int n_outgoing = 2; // TODO: hardcoded due to compiler not allowing dynamic array creation, really needs to be (nexternal - ninitial);
      double z[n_outgoing], r[4] = {0.,0.,0.,0.}, b[3] = {0.,0.,0.}, p2[n_outgoing];
      double xm2[n_outgoing], e[n_outgoing], v[n_outgoing];
      for(int j = 0;j<n_outgoing;++j){
      	p2[j] = 0;
      	xm2[j] = 0;
      	e[j] = 0;
      	v[j] = 0;
      }
      constexpr int itmax = 6;
      const double po2log = log(PI / 2.);

      double q[n_outgoing][4];
      double p[n_outgoing][4];

      // initialization step: factorials for the phase space weight
      z[0] = 0;
      z[1] = po2log;
      for (int j = 2; j < n_outgoing; j++)
        z[j] = z[j - 1] + po2log - 2. * log(double(j - 1));
      for (int j = 2; j < n_outgoing; j++)
        z[j] = (z[j] - log(double(j)));

      // check whether total energy is sufficient; count nonzero masses
      double xmt = 0.;
      int nm = 0;
      for (int j = 0; j < n_outgoing; ++j) {
        if (masses(j+ninitial) != 0.)
          nm = nm + 1;
        xmt = xmt + abs(masses(j+ninitial));
      }
      if (xmt > energy) {
        printf("Too low energy: %f needed %f\n",energy,xmt);
        return;
      }

      // the parameter values are now accepted

      // generate n massless momenta in infinite phase space
      for (int j = 0; j < n_outgoing; j++) {
        double r1 = random_numbers(i,j*4);
        double c = 2. * r1 - 1.;
        double s = sqrt(1. - c * c);
        double f = 2. * PI * random_numbers(i,j*4+1);
        r1 = random_numbers(i,j*4+2);
        double r2 = random_numbers(i,j*4+3);
        q[j][0] = -log(r1 * r2);
        q[j][3] = q[j][0] * c;
        q[j][2] = q[j][0] * s * cos(f);
        q[j][1] = q[j][0] * s * sin(f);
      }

      // calculate the parameters of the conformal transformatjon
      for (int j = 0; j < n_outgoing; j++) {
        for (int k = 0; k < 4; k++)
          r[k] = r[k] + q[j][k];
      }

      double rmas = sqrt(pow(r[0], 2) - pow(r[3], 2) - pow(r[2], 2) - pow(r[1], 2));
      for (int j = 1; j < 4; j++)
        b[j - 1] = -r[j] / rmas;
      double g = r[0] / rmas;
      double a = 1. / (1. + g);
      double x = energy / rmas;
      
      // transform the q's conformally into the p's
      for (int j = 0; j < n_outgoing; j++) {
        double bq = b[0] * q[j][1] + b[1] * q[j][2] + b[2] * q[j][3];
        for (int k = 1; k < 4; k++)
          p[j][k] = x * (q[j][k] + b[k - 1] * (q[j][0] + a * bq));
        p[j][0] = x * (g * q[j][0] + bq);
      }

      // calculate weight and possible warnings
      d_wgt() = po2log;
      if (n_outgoing != 2)
        d_wgt() = (2. * n_outgoing - 4.) * log(energy) + z[n_outgoing - 1];
      if (d_wgt() < -180.) {
        if (iwarn(0) <= 5)
        	printf("Too small wt, risk for underflow: %f\n",d_wgt());
        Kokkos::atomic_fetch_add(&iwarn(0),1);
      }
      if (d_wgt() > 174.) {
        if (iwarn(1) <= 5)
          printf("Too large wt, risk for overflow: %f\n",d_wgt());
        Kokkos::atomic_fetch_add(&iwarn(1),1);
      }

      // return for weighted massless momenta
      if (nm == 0) {
        for (int j = 0; j < n_outgoing; j++)
          for (int k = 0; k < 4; k++)
            d_p(i,j+ninitial,k) = p[j][k];
        return;
      }

      // massive particles: rescale the momenta by a factor x
      double xmax = sqrt(1. - pow(xmt / energy, 2));
      for (int j = 0; j < n_outgoing; j++) {
        xm2[j] = pow(masses(j+ninitial), 2);
        p2[j] = pow(p[j][0], 2);
      }

      int iter = 0;
      x = xmax;
      double accu = energy * ACC;
      while (true) {
        double f0 = -energy;
        double g0 = 0.;
        double x2 = x * x;
        for (int j = 0; j < n_outgoing; j++) {
          e[j] = sqrt(xm2[j] + x2 * p2[j]);
          f0 = f0 + e[j];
          g0 = g0 + p2[j] / e[j];
        }
        if (abs(f0) <= accu)
          break;
        iter = iter + 1;
        if (iter > itmax) {
          printf("Too many iterations without desired accuracy: %d\n", itmax);
          break;
        }
        x = x - f0 / (x * g0);
      }

      for (int j = 0; j < n_outgoing; j++) {
        v[j] = x * p[j][0];
        for (int k = 0; k < 4; k++){
          if(k>0)
            p[j][k] = x * p[j][k];
          printf("i = %03d  p[%03d][%03d] = %10.5f\n",i,j,k,p[j][k]);
          d_p(i,j+ninitial,k) = p[j][k];
        }
        p[j][0] = e[j];
      }

      // calculate the mass-effect weight factor
      double wt2 = 1.;
      double wt3 = 0.;
      for (int j = 0; j < n_outgoing; j++) {
        wt2 = wt2 * v[j] / e[j];
        wt3 = wt3 + pow(v[j], 2) / e[j];
      }
      double wtm = (2. * n_outgoing - 3.) * log(x) + log(wt2 / wt3 * energy);

      // return for  weighted massive momenta
      d_wgt() = d_wgt() + wtm;
      if (d_wgt() < -180.) {
        if (iwarn(2) <= 5)
    	   printf("Too small wt, risk for underflow: %f\n",d_wgt());
        Kokkos::atomic_fetch_add(&iwarn(2),1);
      }
      if (d_wgt() > 174.) {
        if (iwarn(3) <= 5)
    	  printf("Too large wt, risk for overflow: %f\n",d_wgt());
        Kokkos::atomic_fetch_add(&iwarn(3),1);
      }

      for (int j = 0; j < n_outgoing; j++)
        for (int k = 0; k < 4; k++)
          d_p(i,j+ninitial,k) = p[j][k];

    });

}