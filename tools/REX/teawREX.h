/***
 *     _            ______          
 *    | |           | ___ \         
 *    | |_ ___  __ _| |_/ /_____  __
 *    | __/ _ \/ _` |    // _ \ \/ /
 *    | ||  __/ (_| | |\ \  __/>  < 
 *     \__\___|\__,_\_| \_\___/_/\_\
 *                                              
 ***/
//
// *t*ensorial *e*vent *a*daption with *Rex* Version 0.9.0
// teaRex is an extension to the Rex C++ library for parsing and manipulating Les Houches Event-format (LHE) files,
// designed for leading order event reweighting based on input LHE file(s) and scattering amplitude functions.
// teaRex is in development and may not contain all features necessary for all desired features,
// and does not have documentation beyond the code itself.
//
// Copyright © 2023-2024 CERN, CERN Author Zenny Wettersten. 
// Licensed under the GNU Lesser General Public License (version 3 or later).
// All rights not expressly granted are reserved.
//

#ifndef _TEAREX_H_
#define _TEAREX_H_

#include <unistd.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <variant>
#include <stdarg.h>
#include "REX.h"

namespace REX::teaw
{

    using amplitude = std::function<std::shared_ptr<std::vector<double>>(std::vector<double>&, std::vector<double>&)>;
    using vecMap = std::map<REX::event, std::shared_ptr<std::vector<double>>, REX::eventComp>;

    using iterator = std::function<bool()>;
    using weightor = std::function<std::shared_ptr<std::vector<double>>( REX::procSoA& )>;

    // ZW: wrapper for procSoA and the reweighting process
    // Not a derived class of procSoA, such that the full reweight
    // class can be derived from lheSoA without having to deal with
    // the procSoA objects stored in the lheSoA class
    struct procRwgt {
        std::shared_ptr<REX::procSoA> proc;
        std::vector<weightor> amplitude;
        weightor originalAmp;
        std::shared_ptr<std::vector<double>> invOriginalAmp;
        std::vector<std::shared_ptr<std::vector<double>>> weights;
        std::shared_ptr<std::vector<double>> backlog; // FIXME: this should be a vector of vectors if ampsPerIter > 1
        size_t iteration;
        void uniqueReset();
        void reset();
        procRwgt();
        procRwgt( const REX::relEvArgs& relData );
        procRwgt( const procSoA& process );
        procRwgt( procSoA* process );
        procRwgt( std::shared_ptr<procSoA> process );
        procRwgt( const procRwgt& process );
        procRwgt( procRwgt* process );
        procRwgt( std::shared_ptr<procRwgt> process );
        procRwgt( std::vector<std::shared_ptr<REX::event>> lheFile, REX::relEvArgs relArgs = REX::relEvArgs(), 
            std::vector<int> relevStats = {-1,1}, std::function<bool( REX::event& )> relFcn = nullptr );
        procRwgt( std::vector<weightor> ampFcns );
        procRwgt& setAmplitude( weightor amp );
        procRwgt& setAmplitude( std::vector<weightor> amp );
        procRwgt& setOriginalAmp( weightor amp );
        bool initialise();
        bool evaluate();
        bool normalise();
        void doBacklog( bool pass );
    };

    struct reweightor : REX::lheSoA {
        std::vector<std::shared_ptr<procRwgt>> amps;
        std::vector<iterator> iterators;
        iterator originator;
        iterator terminator;
        std::vector<std::shared_ptr<std::vector<double>>> wgts;
        std::vector<double> xSecs;
        std::vector<double> xErrs;
        std::vector<bool> success;
        double ampNorm;
        size_t ampsPerIter = 1;
        void uniqueReset();
        void reset() override;
        void setAmpsFromSubprocs();
        reweightor();
        reweightor( const lheSoA& lheFile );
        reweightor( const reweightor& lheFile );
        reweightor( std::vector<std::shared_ptr<event>> lheFile );
        reweightor( std::vector<std::shared_ptr<event>> lheFile, std::vector<std::function<bool( event& )>> evSort );
        reweightor( lheNode& lheFile, std::vector<std::function<bool( event& )>> evSort );
        reweightor( std::vector<std::shared_ptr<event>> lheFile, std::vector<std::function<bool( event& )>> evSort, 
            std::vector<relEvArgs> relData );
        reweightor( lheNode& lheFile, std::vector<std::function<bool( event& )>> evSort, 
            std::vector<relEvArgs> relData );
        reweightor( std::vector<std::shared_ptr<event>> lheFile, std::vector<std::function<bool( event& )>> evSort, 
            std::vector<std::vector<int>> relStats );
        reweightor( lheNode& lheFile, std::vector<std::function<bool( event& )>> evSort,
            std::vector<std::vector<int>> relStats );
        reweightor( std::vector<std::shared_ptr<event>> lheFile, std::vector<std::function<bool( event& )>> evSort, 
            std::vector<relEvArgs> relData, std::vector<std::vector<int>> relStats );
        reweightor( lheNode& lheFile, std::vector<std::function<bool( event& )>> evSort,
            std::vector<relEvArgs> relData, std::vector<std::vector<int>> relStats );   
        reweightor& setAmps( std::vector<std::shared_ptr<procRwgt>> newAmps ); 
        reweightor& setAmps( std::vector<weightor> ampFcns );
        reweightor& setAmps( std::vector<std::vector<weightor>> ampFcns );
        reweightor& setIterators( std::vector<iterator> iters );
        reweightor& setOriginator( iterator origin );
        reweightor& setTerminator( iterator term );
        // wrapper for running amps (skips iteration if not all amps succeed?)
        bool runAmps();
        // wrapper for running iterator and then amps
        void runBacklog( bool success );
        bool doIterate( size_t index );
        void doAllIterations();
        void doInit();
        void doFin();
        void flattenWeights();
        void calcAmpNorm();
        bool setAmpNorm( bool hard = false );
        void calcXSecs();
        void calcXErrs();
        void run();
        void appendWgtsSimple( lheNode& lheFile, std::vector<std::string_view> procs, std::shared_ptr<std::vector<std::string>> names );
        void appendWgts( lheNode& lheFile, std::vector<std::string_view> procs, std::shared_ptr<std::vector<std::string>> names );
    };

    struct rwgtVal : REX::paramVal{
    public:
        std::string_view blockName;
        bool allStat;
        bool isAll();
        rwgtVal();
        rwgtVal( std::string_view paramLine );
        std::string_view getLine();
        void outWrite( REX::paramBlock& srcBlock );
    };

    struct rwgtBlock {
    public:
        std::string name;
        std::vector<rwgtVal> rwgtVals;
        rwgtBlock( std::vector<std::string_view> values = {}, std::string_view title = "" );
        rwgtBlock( const std::vector<rwgtVal>& vals, std::string_view title = "" );
        std::string_view getBlock();
        //void outWrite( REX::paramBlock& srcBlock, const std::map<std::string_view, int>& blocks );
        void outWrite( REX::paramBlock& srcBlock );
    protected:
        std::string runBlock;
        bool written = false;
    };

    struct rwgtProc {
    public:
        std::vector<rwgtBlock> rwgtParams;
        std::string_view procString;
        std::string_view rwgtName;
        std::vector<std::string_view> rwgtOpts;
        void parse();
        rwgtProc( REX::lesHouchesCard slhaSet, std::string_view rwgtSet = "", bool parseOnline = false );
        std::shared_ptr<REX::lesHouchesCard> outWrite( const REX::lesHouchesCard& paramOrig );
        std::string_view comRunProc();
    };

    struct rwgtCard{
    public:
        REX::lesHouchesCard slhaCard;
        std::vector<iterator> cardWriters;
        std::vector<rwgtProc> rwgtRuns;
        std::vector<std::string_view> rwgtProcs;
        std::vector<std::string_view> opts;
        std::vector<std::shared_ptr<REX::lesHouchesCard>> writtenCards;
        std::shared_ptr<std::vector<std::string>> rwgtNames;
        std::string_view srcCard;
        void parse( bool parseOnline = false );
        rwgtCard( std::string_view reweight_card );
        rwgtCard( std::string_view reweight_card, REX::lesHouchesCard slhaParams, bool parseOnline = false );
        std::vector<std::shared_ptr<REX::lesHouchesCard>> writeCards( REX::lesHouchesCard& slhaOrig );
        std::shared_ptr<std::vector<std::string>> getNames();
        std::vector<std::string_view> getProcs();
        std::vector<iterator> getIterators( std::string path );
    };

    
    struct rwgtCollection {
    public:
        void setRwgt( std::shared_ptr<rwgtCard> rwgts );
        void setRwgt( rwgtCard rwgts );
        void setSlha( std::shared_ptr<REX::lesHouchesCard> slha );
        void setSlha( REX::lesHouchesCard slha );
        void setLhe( std::shared_ptr<REX::lheNode> lhe );
        void setLhe( REX::lheNode& lhe );
        void setLhe( std::string_view lhe_file );
        std::shared_ptr<rwgtCard> getRwgt();
        std::shared_ptr<REX::lesHouchesCard> getSlha();
        std::shared_ptr<REX::lheNode> getLhe();
        REX::transSkel& getSkeleton();
        REX::transSkel& getSkeleton( std::vector<REX::eventSet>& evSets );
        rwgtCollection();
        rwgtCollection( std::shared_ptr<REX::lheNode> lhe, std::shared_ptr<REX::lesHouchesCard> slha, std::shared_ptr<rwgtCard> rwgts );
        rwgtCollection( const rwgtCollection& rwgts );
        std::shared_ptr<std::vector<std::string>> getNames();
    protected:
        template<class... Args>
        void setDoubles(Args&&... args);
        void setSkeleton( std::vector<REX::eventSet>& evSets);
        void setDoublesFromSkeleton();
        std::shared_ptr<rwgtCard> rwgtSets;
        std::shared_ptr<REX::lesHouchesCard> slhaParameters;
        std::shared_ptr<REX::lheNode> lheFile;
        std::vector<std::shared_ptr<std::vector<double>>> wgts;
        std::vector<std::shared_ptr<std::vector<double>>> gS;
        std::vector<std::shared_ptr<std::vector<double>>> momenta;
        std::shared_ptr<std::vector<double>> flatWgts;
        bool lheFileSet = false;
        bool slhaSet = false;
        bool rwgtSet = false;
        bool skeleton = false;
        bool doublesSet = false;
        REX::transLHE eventFile;
        REX::transSkel lheSkeleton;
    };

    struct rwgtFiles : rwgtCollection {
        void setRwgtPath( std::string_view path );
        void setSlhaPath( std::string_view path );
        void setLhePath( std::string_view path );
        rwgtFiles();
        rwgtFiles( std::string_view lhe_card, std::string_view slha_card, std::string_view reweight_card );
        rwgtFiles( const rwgtFiles& rwgts );
        REX::transSkel& initCards( std::vector<REX::eventSet>& evSets);
        void initDoubles();
        template<class... Args>
        void initCards(Args&&... args);
        template<class... Args>
        void initCards( std::string_view lhe_card, std::string_view slha_card, std::string_view reweight_card, Args&&... args );
    protected:
        bool rwgtPulled();
        bool slhaPulled();
        bool lhePulled();
        void pullRwgt();
        void pullSlha();
        void pullLhe();
        bool initialised = false;
        std::string rwgtPath;
        std::string lhePath;
        std::string slhaPath;
        std::shared_ptr<std::string> lheCard;
        std::shared_ptr<std::string> slhaCard;
        std::shared_ptr<std::string> rewgtCard;
    };

    struct rwgtRunner : rwgtFiles{
    public:
        void setMeEval( amplitude eval );
        void addMeEval( const REX::event& ev, const amplitude& eval );
        rwgtRunner();
        rwgtRunner( rwgtFiles& rwgts );
        rwgtRunner( rwgtFiles& rwgts, amplitude meCalc );
        rwgtRunner( rwgtFiles& rwgts, std::vector<amplitude>& meCalcs );
        rwgtRunner( std::string_view lhe_card, std::string_view slha_card, std::string_view reweight_card,
        amplitude meCalc );
        rwgtRunner(const rwgtRunner& rwgts);
        bool oneME();
        bool singAmp();
    protected:
        bool meInit = false;
        bool meCompInit = false;
        bool meSet = false;
        bool normWgtSet = false;
        amplitude meEval;
        //ampCall meEvals;
        std::vector<amplitude> meVec;
        std::vector<std::shared_ptr<std::vector<double>>> initMEs;
        std::vector<std::shared_ptr<std::vector<double>>> meNormWgts;
        std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> reWgts;
        std::shared_ptr<std::vector<double>> normWgt;
        double ampNorm = 0.0;
        std::shared_ptr<std::vector<double>> normXSecs;
        std::shared_ptr<std::vector<double>> errXSecs;
        std::shared_ptr<REX::weightGroup> rwgtGroup;
        template<class... Args>
        void setMEs(Args&&... args);
        void setAmpNorm( double precision );
        bool setParamCard( std::shared_ptr<REX::lesHouchesCard> slhaParams );
        void setNormWgtsSingleME();
        void setNormWgtsMultiME();
        bool calcXSecs();
        bool calcXErrs();
        template<class... Args>
        void setNormWgts(Args&&... args);
        bool singleRwgtIter( std::shared_ptr<REX::lesHouchesCard> slhaParams, std::shared_ptr<REX::lheNode> lheFile, size_t currId );
        bool singleRwgtIter( std::shared_ptr<REX::lesHouchesCard> slhaParams, std::shared_ptr<REX::lheNode> lheFile, size_t currId, std::string& id );
        bool singleRwgtIter( std::shared_ptr<REX::lesHouchesCard> slhaParams, std::shared_ptr<REX::lheNode> lheFile, size_t currId, REX::event& ev );
        bool singleRwgtIter( std::shared_ptr<REX::lesHouchesCard> slhaParams, std::shared_ptr<REX::lheNode> lheFile, size_t currId, 
        std::string& id, REX::event& ev );
        bool lheFileWriter( std::shared_ptr<REX::lheNode> lheFile, std::string outputDir = "rwgt_evts.lhe" );
    public:
        void runRwgt( const std::string& output, double precision = 1e-6 );
        std::shared_ptr<std::vector<double>> getReXSecs();
        std::shared_ptr<std::vector<double>> getReXErrs();
    };


    void rwgtRun( rwgtRunner& rwgt, const std::string& path );

}

#endif