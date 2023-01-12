      // PART 0 - INITIALISATION (before calculate_wavefunctions)
      // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      fptype_sv allMEs = FPZERO_SV;
      #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
          fptype_sv allNumerators = FPZERO_SV;
          fptype_sv allDenominators = FPZERO_SV;
      #endif

      // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
      // (using precomputed good helicities)
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      
      fptype_sv jamp2_sv[ncolor] = { 0 }; // Running sum of partial amplitudes squared for event by event color selection (#402)
      fptype_sv MEs_ighel[ncomb] = { 0 }; // sum of MEs for all good helicities up to ighel (for this event)
      for (size_t ighel = 0; ighel < cNGoodHel[0]; ighel++) {
          const size_t ihel = cGoodHel[ighel];
          #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
              allMEs += calculate_wavefunctions( allmomenta, &allNumerators, &allDenominators, channelId, cHel + ihel*npar, COUPs, cIPD, jamp2_sv );
          #else
              allMEs += calculate_wavefunctions( allmomenta, cHel + ihel*npar, COUPs, cIPD, jamp2_sv );
          #endif
      }

      // Event-by-event random choice of helicity #403
      bool_sv selhel_unset = bool_sv(true);
      for (size_t ighel = 0; ighel < cNGoodHel[0]; ighel++) {
          if (FPANY_SV(selhel_unset)) {
              bool_sv selhel_flip = selhel_unset && (rndhel[0] < (MEs_ighel[ighel]/MEs_ighel[cNGoodHel[0] - 1]));
              selhel[0] = FPCONDITIONAL_SV(selhel[0], int_sv(cGoodHel[ighel] + 1), selhel_flip);
              selhel_unset = selhel_unset && !(selhel_flip);
          }
          else {
              break;
          }
      }

      #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
          auto l_icolamp = mgOnGpu::icolamp<bool>;

          // Event-by-event random choice of color #402
          const size_t channelIdC = channelId - 1; // coloramps.h uses the C array indexing starting at 0
          fptype_sv targetamp[ncolor];
          for (size_t icolC = 0; icolC < ncolor; icolC++) {
              if (icolC == 0) {
                  targetamp[icolC] = FPZERO_SV;
              }
              else {
                  targetamp[icolC] = targetamp[icolC - 1];
              }
              if (l_icolamp[ncolor*channelIdC + icolC]) { targetamp[icolC] += jamp2_sv[icolC]; }
          }

          bool_sv selcol_unset = bool_sv(true);
          for (size_t icolC = 0; icolC < ncolor; icolC++) {
              if (FPANY_SV(selcol_unset)) {
                  bool_sv selcol_flip = selcol_unset && (rndcol[0] < (targetamp[icolC]/targetamp[ncolor - 1]));
                  selcol[0] = FPCONDITIONAL_SV(selcol[0], int_sv(icolC + 1), selcol_flip);
                  selcol_unset = selcol_unset && !(selcol_flip);
              }
              else {
                  break;
              }
          }
      #endif

      // PART 2 - FINALISATION (after calculate_wavefunctions)
      // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
      // [NB 'sum over final spins, average over initial spins', eg see
      // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
          if( channelId > 0 ) allMEs *= allNumerators/allDenominators; // FIXME (#343): assume nprocesses == 1
      #endif
      return allMEs/denominators;
  }

  //--------------------------------------------------------------------------
