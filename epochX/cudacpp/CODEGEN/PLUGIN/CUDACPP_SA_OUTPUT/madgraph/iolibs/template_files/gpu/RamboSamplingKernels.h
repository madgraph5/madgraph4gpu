// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: J. Teig, A. Valassi (2021-2024) for the MG5aMC CUDACPP plugin.

#ifndef RAMBOSAMPLINGKERNELS_H
#define RAMBOSAMPLINGKERNELS_H 1

#include "mgOnGpuConfig.h"

#include "MemoryBuffers.h"

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  // A base class encapsulating phase space sampling on a CPU host or on a GPU device
  class SamplingKernelBase //: virtual public ISamplingKernel
  {
  protected:

    // Constructor from existing input and output buffers
    SamplingKernelBase( const fptype energy,               // input: energy
                        const BufferRndNumMomenta& rndmom, // input: random numbers in [0,1]
                        BufferMomenta& momenta,            // output: momenta
                        BufferWeights& weights )           // output: weights
      : m_energy( energy )
      , m_rndmom( rndmom )
      , m_momenta( momenta )
      , m_weights( weights )
    {
    }

  public:

    // Destructor
    virtual ~SamplingKernelBase() {}

    // Get momenta of initial state particles
    virtual void getMomentaInitial() = 0;

    // Get momenta of final state particles and weights
    virtual void getMomentaFinal() = 0;

    // Is this a host or device kernel?
    virtual bool isOnDevice() const = 0;

  protected:

    // The energy
    const fptype m_energy;

    // The buffer for the input random numbers
    const BufferRndNumMomenta& m_rndmom;

    // The buffer for the output momenta
    BufferMomenta& m_momenta;

    // The buffer for the output weights
    BufferWeights& m_weights;
  };

  //--------------------------------------------------------------------------

  // A class encapsulating RAMBO phase space sampling on a CPU host
  class RamboSamplingKernelHost final : public SamplingKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    RamboSamplingKernelHost( const fptype energy,               // input: energy
                             const BufferRndNumMomenta& rndmom, // input: random numbers in [0,1]
                             BufferMomenta& momenta,            // output: momenta
                             BufferWeights& weights,            // output: weights
                             const size_t nevt );

    // Destructor
    virtual ~RamboSamplingKernelHost() {}

    // Get momenta of initial state particles
    void getMomentaInitial() override final;

    // Get momenta of final state particles and weights
    void getMomentaFinal() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return false; }
  };

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  // A class encapsulating RAMBO phase space sampling on a GPU device
  class RamboSamplingKernelDevice final : public SamplingKernelBase, public NumberOfEvents
  {
  public:

    // Constructor from existing input and output buffers
    RamboSamplingKernelDevice( const fptype energy,               // input: energy
                               const BufferRndNumMomenta& rndmom, // input: random numbers in [0,1]
                               BufferMomenta& momenta,            // output: momenta
                               BufferWeights& weights,            // output: weights
                               const size_t gpublocks,
                               const size_t gputhreads );

    // Destructor
    virtual ~RamboSamplingKernelDevice() {}

    // Get momenta of initial state particles
    void getMomentaInitial() override final;

    // Get momenta of final state particles and weights
    void getMomentaFinal() override final;

    // Is this a host or device kernel?
    bool isOnDevice() const override final { return true; }

  private:

    // The number of blocks in the GPU grid
    size_t m_gpublocks;

    // The number of threads in the GPU grid
    size_t m_gputhreads;
  };
#endif

  //--------------------------------------------------------------------------
}
#endif // RAMBOSAMPLINGKERNELS_H
