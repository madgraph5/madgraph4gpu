# This file was automatically created  by FeynRules 2.0 (static file)
# Mathematica version: 8.0 for Mac OS X x86 (64-bit) (November 6, 2010)
# Date: Mon 1 Oct 2012 14:58:26

from object_library import all_propagators, Propagator


# define only once the denominator since this is always the same
denominator = "(P('mu', id) * P('mu', id) - Mass(id) * Mass(id) + complex(0,1) * Mass(id) * Width(id))"
denominatorSq = denominator + "**2"

numV = "(- Metric(1, 2) + Metric(1,'mu')* P('mu', id) * P(2, id) / Mass(id)**2) "
numF = "(Gamma('mu', 1, 2) * P('mu', id) + Mass(id) * Identity(1, 2))"

# propagator for the scalar
S = Propagator(name = "S",
               numerator = "1",
               denominator = denominator
               )

# propagator for the incoming fermion # the one for the outcomming is computed on the flight
F = Propagator(name = "F",
                numerator = numF,
                denominator = denominator
               )

# massive vector in the unitary gauge, can't be use for massless particles
V1 = Propagator(name = "V1",
                numerator = numV,
                denominator = denominator
               )

# massless vector and massive vector in unitary gauge
V2 = Propagator(name = "V2",
                numerator = "- Metric(1, 2)",
                denominator =  "P('mu', id) * P('mu', id)"
               )

# propagators with width corrections
Z1 =  Propagator(name = "Z1",
                numerator = "-" +  numV + "* complex(0,1) * Mass(id) * dWZ",
                denominator = denominatorSq
               )

W1 =  Propagator(name = "W1",
                numerator = "-" +  numV + "* complex(0,1) * Mass(id) * dWW",
                denominator = denominatorSq
               )

H1 =  Propagator(name = "H1",
                numerator = "- complex(0,1) * Mass(id) * dWH",
                denominator = denominatorSq
               )

T1 =  Propagator(name = "T1",
                numerator = "-" +  numF + "* complex(0,1) * Mass(id) * dWT",
                denominator = denominatorSq
               )

