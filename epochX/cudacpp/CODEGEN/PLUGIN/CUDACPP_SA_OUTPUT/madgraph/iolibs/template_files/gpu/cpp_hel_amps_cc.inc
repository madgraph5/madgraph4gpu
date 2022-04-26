%(function_definitions)s
  // Compute the output couplings gc10 and gc11 from the input gs
  template<class G_ACCESS, class C_ACCESS>
  __device__ void
  G2COUP( const fptype gs[],
          fptype couplings[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& gs_sv = G_ACCESS::kernelAccessConst( gs );
    fptype* gc10s = C_ACCESS::idcoupAccessBuffer( couplings, Parameters_sm_dependentCouplings::idcoup_GC_10 );
    fptype* gc11s = C_ACCESS::idcoupAccessBuffer( couplings, Parameters_sm_dependentCouplings::idcoup_GC_11 );
    cxtype_sv_ref gc10s_sv = C_ACCESS::kernelAccess( gc10s );
    cxtype_sv_ref gc11s_sv = C_ACCESS::kernelAccess( gc11s );
    gc10s_sv = Parameters_sm_dependentCouplings::GC_10_fromG( gs_sv );
    gc11s_sv = Parameters_sm_dependentCouplings::GC_11_fromG( gs_sv );
    /*
    printf( "G2COUP: pgs=%%p pgc10s=%%p pgc11s=%%p\n", gs, gc10s, gc11s );
#ifdef __CUDACC__
    printf( "G2COUP: pgssv=%%p pgc10ssv=%%p pgc11ssv=%%p\n", &gs_sv, &gc10s_sv, &gc11s_sv );
    //printf( "G2COUP: pgc10ssvr=%%p pgc11ssvr=%%p\n", gc10s_sv.m_preal, gc11s_sv.m_preal ); // needs a hack in cxtype_ref
    printf( "G2COUP: %%s=%%2d gs=%%f gc10s=(%%f, %%f) gc11s=(%%f, %%f)\n",
            "ievt", blockDim.x * blockIdx.x + threadIdx.x, gs_sv,
            cxreal( gc10s_sv ), cximag( gc10s_sv ), cxreal( gc11s_sv ), cximag( gc11s_sv ) );
#else
#ifndef MGONGPU_CPPSIMD
    printf( "G2COUP: %%s=%%2d gs=%%f gc10s=(%%f, %%f) gc11s=(%%f, %%f)\n",
            "ieppV", 0, gs_sv,
            cxreal( gc10s_sv ), cximag( gc10s_sv ), cxreal( gc11s_sv ), cximag( gc11s_sv ) );
#else
    for( int ieppV = 0; ieppV < neppV; ieppV++ )
      printf( "G2COUP: %%s=%%2d gs=%%f gc10s=(%%f, %%f) gc11s=(%%f, %%f)\n", "ieppV",
              ieppV, gs_sv[ieppV],
              cxreal( gc10s_sv )[ieppV], cximag( gc10s_sv )[ieppV],
              cxreal( gc11s_sv )[ieppV], cximag( gc11s_sv )[ieppV] );
#endif
#endif
    */
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

} // end namespace

#endif // HelAmps_%(model_name)s_H
