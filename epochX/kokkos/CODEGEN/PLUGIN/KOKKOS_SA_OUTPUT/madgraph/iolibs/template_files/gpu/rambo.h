#ifndef RAMBO_H
#define RAMBO_H 1
#include <vector>
#include "Kokkos_Core.hpp"
#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

#define ACC 1e-14


void get_initial_momenta(
      Kokkos::View<double***> d_p,
      const int nexternal,const double energy,
      const int& league_size,
      const int& team_size)
{
  
  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(const member_type& team_member){
    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    const fptype energy1 = energy/2;
    const fptype energy2 = energy/2;
    const fptype mom = energy/2;

    // particle 1
    d_p(ievt,0,0) = energy1;
    d_p(ievt,0,1) = 0.;
    d_p(ievt,0,2) = 0.;
    d_p(ievt,0,3) = mom;
    // particle 2
    d_p(ievt,1,0) = energy2;
    d_p(ievt,1,1) = 0.;
    d_p(ievt,1,2) = 0.;
    d_p(ievt,1,3) = -mom;
  });
}

void get_final_momenta(const int ninitial,const int nexternal,const double energy,
    // const Kokkos::View<double*>& masses,
    Kokkos::View<double ***>& d_p,
    const Kokkos::View<double **>& random_numbers,
    Kokkos::View<double*>& d_wgt,
    const int& league_size,
    const int& team_size)
  {

  Kokkos::View<int*> iwarn("iwarn",5);
  
  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
    KOKKOS_LAMBDA(const member_type& team_member){
      const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();

      auto lp = Kokkos::subview(d_p,ievt,Kokkos::ALL,Kokkos::ALL);

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
      
      // initialization step: factorials for the phase space weight
      const fptype twopi = 8. * atan(1.);
      const fptype po2log = log(twopi / 4.);
      fptype z[mgOnGpu::nparf];
      z[1] = po2log;
      for (int kpar = 2; kpar < mgOnGpu::nparf; kpar++)
        z[kpar] = z[kpar - 1] + po2log - 2. * log(fptype(kpar - 1));
      for (int kpar = 2; kpar < mgOnGpu::nparf; kpar++)
        z[kpar] = (z[kpar] - log(fptype(kpar)));

      fptype& wt = d_wgt[ievt];

      // generate n massless momenta in infinite phase space
      fptype q[mgOnGpu::nparf][mgOnGpu::np4];
      for (int iparf = 0; iparf < mgOnGpu::nparf; iparf++) {
        const fptype r1 = random_numbers(ievt,iparf*mgOnGpu::nparf + 0);
        const fptype r2 = random_numbers(ievt,iparf*mgOnGpu::nparf + 1);
        const fptype r3 = random_numbers(ievt,iparf*mgOnGpu::nparf + 2);
        const fptype r4 = random_numbers(ievt,iparf*mgOnGpu::nparf + 3);
        // printf("ievt=%02d  iparf = %02d  r = %6.3f,%6.3f,%6.3f,%6.3f\n",ievt,iparf,r1,r2,r3,r4);
        
        const fptype c = 2. * r1 - 1.;
        const fptype s = sqrt(1. - c * c);
        const fptype f = twopi * r2;
        q[iparf][0] = -log(r3 * r4);
        q[iparf][3] = q[iparf][0] * c;
        q[iparf][2] = q[iparf][0] * s * cos(f);
        q[iparf][1] = q[iparf][0] * s * sin(f);
      }

      // calculate the parameters of the conformal transformation
      fptype r[mgOnGpu::np4];
      fptype b[mgOnGpu::np4-1];
      for (int i4 = 0; i4 < mgOnGpu::np4; i4++)
        r[i4] = 0.;
      for (int iparf = 0; iparf < mgOnGpu::nparf; iparf++) {
        for (int i4 = 0; i4 < mgOnGpu::np4; i4++)
          r[i4] = r[i4] + q[iparf][i4];
      }
      const fptype rmas = sqrt(pow(r[0], 2) - pow(r[3], 2) - pow(r[2], 2) - pow(r[1], 2));
      for (int i4 = 1; i4 < mgOnGpu::np4; i4++)
        b[i4-1] = -r[i4] / rmas;
      const fptype g = r[0] / rmas;
      const fptype a = 1. / (1. + g);
      const fptype x0 = energy / rmas;

      // transform the q's conformally into the p's (i.e. the 'momenta')
      for (int iparf = 0; iparf < mgOnGpu::nparf; iparf++) {
        fptype bq = b[0] * q[iparf][1] + b[1] * q[iparf][2] + b[2] * q[iparf][3];
        for (int i4 = 1; i4 < mgOnGpu::np4; i4++)
          lp(iparf+mgOnGpu::npari,i4) = x0 * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
        lp(iparf+mgOnGpu::npari,0) = x0 * (g * q[iparf][0] + bq);

        // printf("ievt=%02d  iparf = %02d  p = %6.3f,%6.3f,%6.3f,%6.3f\n",ievt,iparf,
        //   lp(iparf+mgOnGpu::npari,0),
        //   lp(iparf+mgOnGpu::npari,1),
        //   lp(iparf+mgOnGpu::npari,2),
        //   lp(iparf+mgOnGpu::npari,3));
      }


      // calculate weight (NB return log of weight)
      wt = po2log;
      if (mgOnGpu::nparf != 2)
        wt = (2. * mgOnGpu::nparf - 4.) * log(energy) + z[mgOnGpu::nparf-1];

      // return for weighted massless momenta
      // nothing else to do in this event if all particles are massless (nm==0)

      return; // TODO Fix for massive partons
/*
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
          //printf("i = %03d  p[%03d][%03d] = %10.5f\n",i,j,k,p[j][k]);
          //lp(j+ninitial,k) = p[j][k];
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
      d_wgt(i) = d_wgt(i) + wtm;
      if (d_wgt(i) < -180.) {
        if (iwarn(2) <= 5)
         printf("Too small wt, risk for underflow: %f\n",d_wgt(i));
        Kokkos::atomic_fetch_add(&iwarn(2),1);
      }
      if (d_wgt(i) > 174.) {
        if (iwarn(3) <= 5)
        printf("Too large wt, risk for overflow: %f\n",d_wgt(i));
        Kokkos::atomic_fetch_add(&iwarn(3),1);
      }

      for (int j = 0; j < n_outgoing; j++)
        for (int k = 0; k < 4; k++)
          lp(j+ninitial,k) = p[j][k];
*/
    }
  );
}

#endif // RAMBO_H 1
