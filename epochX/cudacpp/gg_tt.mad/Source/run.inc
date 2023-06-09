C     **************************************************************
C     run.inc 
C     **************************************************************
C     
C     Scales
C     
      REAL*8          SCALE,SCALEFACT,ALPSFACT,MUE_REF_FIXED
     $ ,MUE_OVER_REF
      LOGICAL         FIXED_REN_SCALE,FIXED_FAC_SCALE1,
     $  FIXED_FAC_SCALE2,FIXED_COUPLINGS,HMULT
      LOGICAL         FIXED_EXTRA_SCALE
      INTEGER         ICKKW,NHMULT,ASRWGTFLAVOR,
     $  DYNAMICAL_SCALE_CHOICE,IEVO_EVA

      COMMON/TO_SCALE/SCALE,SCALEFACT,ALPSFACT, MUE_REF_FIXED,
     $  MUE_OVER_REF,
     $                fixed_ren_scale,fixed_fac_scale1, fixed_fac_scale2,
     $                fixed_couplings, fixed_extra_scale,ickkw,nhmult,hmult,asrwgtflavor,
     $                dynamical_scale_choice
      COMMON/TO_EVA/IEVO_EVA
C     
C     Collider
C     
      INTEGER                                        LPP(2)
      DOUBLE PRECISION    EBEAM(2), XBK(2),Q2FACT(2)
      COMMON/TO_COLLIDER/ EBEAM   , XBK   ,Q2FACT,   LPP

      INTEGER FRAME_ID
      COMMON/TO_FRAME_ME/FRAME_ID
C     
C     Number of events
C     
      INTEGER NEVENTS, GEVENTS
C     
C     BW treatment
C     
      DOUBLE PRECISION BWCUTOFF
      COMMON/TO_BWCUTOFF/ BWCUTOFF
C     
C     Sudakov grid file name
C     
      CHARACTER*130 ISSGRIDFILE
      COMMON/TO_SGRIDFILE/ISSGRIDFILE
C     
C     kT/pT scheme for xqcut, clustering according to channel
C     
      INTEGER KTSCHEME
      LOGICAL CHCLUSTER,PDFWGT
      COMMON/TO_CLUSTER/KTSCHEME,CHCLUSTER,PDFWGT

C     
C     Parameters for systematics variations info
C     
      LOGICAL USE_SYST
C     Common block for systematics variations
      DOUBLE PRECISION S_SCALE
      INTEGER N_QCD,N_ALPSEM
      DOUBLE PRECISION S_QALPS(MAX_PARTICLES-2)
      INTEGER N_PDFRW(2),I_PDGPDF(MAX_PARTICLES-2,2)
      DOUBLE PRECISION S_XPDF(MAX_PARTICLES-2,2),S_QPDF(MAX_PARTICLES
     $ -2,2)
      DOUBLE PRECISION S_RWFACT
      COMMON/TO_SYST/USE_SYST,N_QCD,N_ALPSEM,N_PDFRW,I_PDGPDF,
     $               s_scale,s_qalps,s_xpdf,s_qpdf,s_rwfact
C     
C     Flag on how to write the LHE events 
C     Include <clustering> tag for Pythia 8 CKKW-L matching
C     
      LOGICAL CLUSINFO
      DOUBLE PRECISION LHE_VERSION
      COMMON/TO_LHEFORMAT/LHE_VERSION,CLUSINFO

C     
C     Controls wheter to perform Monte-Carlo sampling over grouped
C      subprocesses
C     
      LOGICAL MC_GROUPED_SUBPROC
      COMMON/TO_MC_GROUPED_SUBPROC/MC_GROUPED_SUBPROC

C     
C     block for heavy ion beam
C     
      INTEGER NB_PROTON(2), NB_NEUTRON(2)
      COMMON/TO_HEAVYION_PDG/ NB_PROTON, NB_NEUTRON
      DOUBLE PRECISION MASS_ION(2)
      COMMON/TO_HEAVYION_MASS/MASS_ION

C     
C     Controls what are the PDGs included in the CKKWl merging
C      procedure, i.e. what
C     are the PDGs subject to the ktdurham cut
C     
      INTEGER PDGS_FOR_MERGING_CUT(0:1000)
      COMMON/TO_MERGE/PDGS_FOR_MERGING_CUT
C     
C     
C     
      INTEGER PDG_CUT(0:25)
      DOUBLE PRECISION PTMIN4PDG(0:25)
      DOUBLE PRECISION PTMAX4PDG(0:25)
      DOUBLE PRECISION EMIN4PDG(0:25)
      DOUBLE PRECISION EMAX4PDG(0:25)
      DOUBLE PRECISION ETAMIN4PDG(0:25)
      DOUBLE PRECISION ETAMAX4PDG(0:25)
      DOUBLE PRECISION MXXMIN4PDG(0:25)
      LOGICAL MXXPART_ANTIPART(1:25)
      COMMON/TO_PDG_SPECIFIC_CUT/PDG_CUT, PTMIN4PDG,PTMAX4PDG,
     $  EMIN4PDG, EMAX4PDG, ETAMIN4PDG,
     &etamax4pdg, mxxmin4pdg,mxxpart_antipart

      DOUBLE PRECISION SMALL_WIDTH_TREATMENT
      COMMON/NARROW_WIDTH/SMALL_WIDTH_TREATMENT

      DOUBLE PRECISION TMIN_FOR_CHANNEL
      INTEGER SDE_STRAT  ! 1 means standard single diagram enhancement strategy,
C     2 means approximation by the denominator of the propagator
      COMMON/TO_CHANNEL_STRAT/TMIN_FOR_CHANNEL, SDE_STRAT
      INTEGER EXEC_MODE  ! added by autodef
      INTEGER BRIDGE_MODE  ! added by autodef
      INTEGER VECTOR_SIZE  ! added by autodef
C     START USER COMMON BLOCK
      COMMON/USER_CUSTOM_RUN/EXEC_MODE,BRIDGE_MODE,VECTOR_SIZE
C     STOP USER COMMON BLOCK
