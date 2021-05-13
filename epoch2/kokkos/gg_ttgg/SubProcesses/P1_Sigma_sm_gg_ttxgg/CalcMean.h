#ifndef __CALCMEAN__H
#define __CALCMEAN__H
#include <cmath>
#include <limits>

template <typename float_precision,typename decimal_precision>
class CalcMean{
private:
   float_precision m_mean;
   float_precision m_sigma;
   float_precision m_sum;
   float_precision m_sum2;
   float_precision m_min;
   float_precision m_max;
   decimal_precision m_n;
public:
   CalcMean(void):m_mean(0.),m_sigma(0.),
      m_sum(0.),m_sum2(0.),m_n(0),m_min(std::numeric_limits<float_precision>::max()),
      m_max(std::numeric_limits<float_precision>::min()){};

   void add_value(float_precision value){
      m_n += 1;
      m_sum += value;
      m_sum2 += value*value;
      m_mean = 0;
      m_sigma = 0;
      if(value < m_min) m_min = value;
      if(value > m_max) m_max = value;
   }

   float_precision mean(){
      if(m_mean != 0)
         return m_mean;
      if(m_n == 0)
         return 0;

      m_mean = float_precision(m_sum)/float_precision(m_n);
      return m_mean;
   }

   float_precision sigma(){
      if(m_sigma != 0)
         return m_sigma;
      if(m_n == 0)
         return 0;
      m_mean = mean();
      m_sigma = sqrt( (1./m_n)*m_sum2 - m_mean*m_mean);
      return m_sigma;
   }

   decimal_precision n(){return m_n;};
   void n(decimal_precision value){m_n = value;};

   float_precision sum(){return m_sum;};
   void sum(float_precision value){m_sum = value;};

   float_precision sum2(){return m_sum2;}
   void sum2(float_precision value){m_sum2 = value;};

   float_precision min(){return m_min;};
   float_precision max(){return m_max;};

};

#endif
