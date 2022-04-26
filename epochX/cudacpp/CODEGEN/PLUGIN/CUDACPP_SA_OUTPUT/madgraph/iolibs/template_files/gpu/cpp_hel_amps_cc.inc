%(function_definitions)s
  // Compute the output couplings (e.g. gc10 and gc11) from the input gs
  template<class G_ACCESS, class C_ACCESS>
  __device__ void
  G2COUP( const fptype gs[],
          fptype couplings[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const fptype_sv& gs_sv = G_ACCESS::kernelAccessConst( gs );
    //%%(dcoupaccessbuffer)s
    //%%(dcoupkernelaccess)s
    //%%(dcoupcompute)s
    fptype* GC_10s = C_ACCESS::idcoupAccessBuffer( couplings, Parameters_sm_dependentCouplings::idcoup_GC_10 );
    fptype* GC_11s = C_ACCESS::idcoupAccessBuffer( couplings, Parameters_sm_dependentCouplings::idcoup_GC_11 );
    cxtype_sv_ref GC_10s_sv = C_ACCESS::kernelAccess( GC_10s );
    cxtype_sv_ref GC_11s_sv = C_ACCESS::kernelAccess( GC_11s );
    GC_10s_sv = Parameters_sm_dependentCouplings::GC_10_fromG( gs_sv );
    GC_11s_sv = Parameters_sm_dependentCouplings::GC_11_fromG( gs_sv );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

} // end namespace

#endif // HelAmps_%(model_name)s_H
