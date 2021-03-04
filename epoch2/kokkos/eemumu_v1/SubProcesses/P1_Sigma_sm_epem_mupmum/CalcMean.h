#ifndef __CALCMEAN__H
#define __CALCMEAN__H
#include <cmath>
#include <limits>

class CalcMean{
private:
   float m_mean;
   float m_sigma;
   float m_sum;
   float m_sum2;
   float m_min;
   float m_max;
   unsigned int m_n;
public:
   CalcMean(void):m_mean(0.),m_sigma(0.),
      m_sum(0.),m_sum2(0.),m_n(0),m_min(std::numeric_limits<float>::max()),
      m_max(std::numeric_limits<float>::min()){};

   void add_value(float value){
      m_n += 1;
      m_sum += value;
      m_sum2 += value*value;
      m_mean = 0;
      m_sigma = 0;
      if(value < m_min) m_min = value;
      if(value > m_max) m_max = value;
   }

   float mean(){
      if(m_mean != 0)
         return m_mean;
      if(m_n == 0)
         return 0;

      m_mean = float(m_sum)/float(m_n);
      return m_mean;
   }

   float sigma(){
      if(m_sigma != 0)
         return m_sigma;
      if(m_n == 0)
         return 0;
      m_mean = mean();
      m_sigma = sqrt( (1./m_n)*m_sum2 - m_mean*m_mean);
      return m_sigma;
   }

   int n(){return m_n;};
   void n(int value){m_n = value;};

   float sum(){return m_sum;};
   void sum(float value){m_sum = value;};

   float sum2(){return m_sum2;}
   void sum2(float value){m_sum2 = value;};

   float min(){return m_min;};
   float max(){return m_max;};

};

#endif
