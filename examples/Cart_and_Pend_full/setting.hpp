#ifndef SETTING_HPP
#define SETTING_HPP

namespace sac {

  /*********************************************/
  /* Parameters */
  double T = 0.28;  // 1.2; 2.0 // prediction horizon
  double lam = -10; 
  double maxdt = .2;
  double ts = 0.001; // .0167;
  double usat[1][2] = { {25, -25} };
  double calc_tm = ts;
  /**/

  const size_t xlen = 4, ulen = 1;

}

#endif  // SETTING_HPP
