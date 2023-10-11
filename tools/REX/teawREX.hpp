/***
 *     _                     ______ _______   __
 *    | |                    | ___ \  ___\ \ / /
 *    | |_ ___  __ ___      _| |_/ / |__  \ V / 
 *    | __/ _ \/ _` \ \ /\ / /    /|  __| /   \ 
 *    | ||  __/ (_| |\ V  V /| |\ \| |___/ /^\ \
 *     \__\___|\__,_| \_/\_/ \_| \_\____/\/   \/
 *                                              
 ***/

// THIS IS NOT A LICENSED RELEASE
// IF YOU SEE THIS FILE, IT HAS BEEN SPREAD
// FROM AN IMPROPER RELEASE.

// Copyright © 2023 CERN, CERN Author Zenny Wettersten. 
// All rights reserved.

#include <unistd.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <variant>
#include "REX.hpp"

namespace REX::teaw
{
    template<typename T1, typename T2>
    std::shared_ptr<std::vector<T1>> scatAmpEval(std::vector<T2>& momenta, std::function<std::shared_ptr<std::vector<T1>>(std::vector<T2>&)> evalFunc)
    { return evalFunc(momenta); }

    template<typename T1, typename T2>
    std::shared_ptr<std::vector<T1>> scatAmpEval(std::vector<T2>& momenta, std::function<std::vector<T1>(std::vector<T2>&)> evalFunc)
    { return evalFunc(momenta); }

    template<typename T1, typename T2>
    std::shared_ptr<std::vector<T1>> scatAmpEval(std::vector<T2>& momenta, std::function<std::shared_ptr<std::vector<T1>>(std::vector<T2>&, std::vector<T2>&)> evalFunc)
    { return evalFunc(momenta); }

    template<typename T1, typename T2>
    std::shared_ptr<std::vector<T1>> scatAmpEval(std::vector<T2>& momenta, std::function<std::vector<T1>(std::vector<T2>&, std::vector<T2>&)> evalFunc)
    { return evalFunc(momenta); }

    struct rwgtVal : REX::paramVal{
    public:
        std::string_view blockName;
        bool allStat;
        bool isAll(){ return (idStr == "all"); }
        rwgtVal() : paramVal(){ return; }
        rwgtVal( std::string_view paramLine )
        : paramVal( paramLine, false ){if( paramLine.size() == 0 ){ return; }
            realLine = paramLine;
            auto vals = *REX::nuBlankSplitter( realLine );
            blockName = vals[1];
            idStr = vals[2];
            valStr = vals[3];
        }
        std::string_view getLine(){ return realLine; }
        void outWrite( REX::paramBlock& srcBlock ){
            if ( isAll() )
            {
                for( auto param : srcBlock.params )
                {
                    param.valStr = valStr;
                    param.modded = true;
                } 
                return;
            } 
            auto currPar = std::find_if( srcBlock.params.begin(), srcBlock.params.end(), 
            [&]( const REX::paramVal& parPar ){ return (parPar.idStr == idStr ); } );
            if( currPar == srcBlock.params.end() ){ 
                srcBlock.params.push_back( REX::paramVal( realLine.substr(realLine.find("set") + 4) ) );
                srcBlock.params[ srcBlock.params.size() - 1 ].modded = true; 
                srcBlock.modded = true;
                return;
            }
            currPar->valStr = valStr;
            currPar->modded = true;
            srcBlock.modded = true;
            return;
        }
    };

    struct rwgtBlock {
    public:
        std::string_view name;
        std::vector<rwgtVal> rwgtVals;
        rwgtBlock( std::vector<std::string_view> values = {}, std::string_view title = "" )
        {
            name = title;
            rwgtVals.resize( values.size() );
            for( int k = 0 ; k < values.size() ; ++k )
            {
                rwgtVals[k] = rwgtVal( values[k] );
            }
        }
        rwgtBlock( const std::vector<rwgtVal>& vals, std::string_view title = "" )
        {
            name = title;
            rwgtVals = vals;
        }
        std::string_view getBlock(){
            if( written ){ return runBlock; }
            runBlock = "";
            for( auto val : rwgtVals ){
                runBlock += std::string(val.getLine()) + "\n";
            }
            written = true;
            return runBlock;
        }
        void outWrite( REX::paramBlock& srcBlock, const std::map<std::string_view, int>& blocks )
        {
            for( auto parm : rwgtVals )
            {
                parm.outWrite( srcBlock );
            }
            srcBlock.modded = true;
            return;
        }
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
        void parse(){
            std::vector<std::string_view> blocks;
            std::vector<std::shared_ptr<std::vector<rwgtVal>>> params;
            auto procLines = *REX::nuLineSplitter( procString );
            for( auto line : procLines )
            {
                auto strtPt = line.find("set");
                auto words = *REX::nuWordSplitter( line );
                auto currBlock = words[1]; 
                auto loc = std::find_if( blocks.begin(), blocks.end(), 
                [&]( std::string_view block ){ return (block == currBlock); } );
                if( loc == blocks.end() ){ 
                    blocks.push_back( currBlock ); 
                    params.push_back( std::make_shared<std::vector<rwgtVal>>( std::vector<rwgtVal>({rwgtVal( line )} ) )); }
                else { 
                    params[ std::distance( blocks.begin(), loc ) - 1 ]->push_back( rwgtVal( line ) );
                }
            }
            rwgtParams.reserve(blocks.size());
            for( int k = 0 ; k < blocks.size() ; ++k )
            {
                rwgtParams.push_back( rwgtBlock( *params[k], blocks[k] ) );
            }
        }
        rwgtProc( REX::lesHouchesCard slhaSet, std::string_view rwgtSet = "", bool parseOnline = false )
        {
            if( rwgtSet == "" ){ return; }
            auto strtLi = rwgtSet.find( "\n", rwgtSet.find("launch") ) + 1;
            auto endLi = rwgtSet.find("\n", strtLi);
            while( rwgtSet[rwgtSet.find_first_not_of("\n ", endLi)] == 's' )
            { endLi = rwgtSet.find( "\n", endLi + 1 ); }
            procString = rwgtSet.substr( strtLi, endLi - strtLi );
            if( parseOnline ){ parse(); }
        }
        std::shared_ptr<REX::lesHouchesCard> outWrite( const REX::lesHouchesCard& paramOrig ){
            auto slhaOrig = std::make_shared<REX::lesHouchesCard>( paramOrig );
            std::map<std::string_view, int> blockIds;
            for( int k = 0 ; k < slhaOrig->blocks.size() ; ++k )
            {   slhaOrig->blocks[k].parse( true );
                auto nyama = std::pair<std::string_view, int>( slhaOrig->blocks[k].name, k);
                blockIds.insert( nyama ); }
            for( auto rwgts : rwgtParams )
            { rwgts.outWrite( slhaOrig->blocks[ blockIds.at( rwgts.name ) ], blockIds ); }
            slhaOrig->modded = true;
            return slhaOrig;
        }
        std::string_view comRunProc(){ return procString; }
    };

    struct rwgtCard{
    public:
        REX::lesHouchesCard slhaCard;
        std::vector<rwgtProc> rwgtRuns;
        std::vector<std::string_view> rwgtProcs;
        std::vector<std::string_view> opts;
        std::vector<std::string> rwgtNames;
        std::string_view srcCard;
        void parse( bool parseOnline = false ) {
            auto strt = srcCard.find("launch");
            while( auto commPos = srcCard.find_last_of("#", strt) > srcCard.find_last_of("\n", strt) ){ 
                if( commPos == REX::npos ){
                    break;
                }
                strt = srcCard.find("launch", strt + 6 );
            }
            while( auto chPos = srcCard.find( "set" ) < strt ){
                if( srcCard.find_last_of("#", chPos) > srcCard.find_last_of("\n", chPos) ){ chPos = srcCard.find("change", strt + 6 ); continue; }
                opts.push_back( srcCard.substr( chPos, srcCard.find("\n", chPos) - chPos ) );
            }
            std::vector<size_t> lnchPos({strt}); 
            auto nuLnch = srcCard.find( "launch", strt + 6 );
            while ( nuLnch != std::string_view::npos )
            {
                if( srcCard.find_last_of("#", nuLnch) < srcCard.find_last_of("\n", nuLnch) ){ lnchPos.push_back(nuLnch); }
                nuLnch = srcCard.find( "launch", nuLnch + 6 );
            }
            for( int k = 0 ; k < lnchPos.size() - 1 ; ++k )
            {
                auto strtLi = srcCard.find( "set", lnchPos[k] );
                rwgtRuns.push_back( rwgtProc( slhaCard, srcCard.substr( strtLi, lnchPos[k+1] - strtLi ), parseOnline ) );
                if( srcCard.find( "--", lnchPos[k] ) < strtLi ){
                    auto strtPos = srcCard.find( "--", lnchPos[k] );
                    while( (strtPos < strtLi ) && (strtPos!= std::string_view::npos) ){
                        auto nuStrtPos = std::min( srcCard.find( "\n", strtPos ), srcCard.find( "--", strtPos + 1 ));
                        rwgtRuns[ rwgtRuns.size() - 1 ].rwgtOpts.push_back( srcCard.substr( strtPos, nuStrtPos - strtPos ) );
                        if( rwgtRuns[ rwgtRuns.size() - 1 ].rwgtOpts[ rwgtRuns[ rwgtRuns.size() - 1 ].rwgtOpts.size() - 1 ].substr(2,11) == "rwgt_name"){
                            rwgtRuns[ rwgtRuns.size() - 1 ].rwgtName = rwgtRuns[ rwgtRuns.size() - 1 ].
                                rwgtOpts[ rwgtRuns[ rwgtRuns.size() - 1 ].rwgtOpts.size() - 1 ].substr( 11, nuStrtPos - strtPos - 11 );
                        }
                        if( nuStrtPos == srcCard.find( "\n", strtPos ) ){ break; }
                        strtPos = nuStrtPos;
                    }
                }
            }
            size_t endLi = srcCard.find( "\n", lnchPos[ lnchPos.size() - 1 ] );
            if( srcCard.substr( endLi + 1, 3 ) == "set" ){
                while( srcCard.substr( endLi + 1, 3 ) == "set" )
                {
                    endLi = srcCard.find( "\n", endLi + 1 );
                }
                rwgtRuns.push_back( rwgtProc( slhaCard, srcCard.substr( lnchPos[lnchPos.size()-1], endLi - lnchPos[lnchPos.size()-1] ), parseOnline ) );
            }
            rwgtProcs = std::vector<std::string_view>(); rwgtProcs.reserve( rwgtRuns.size() );
            rwgtNames.reserve( rwgtRuns.size() );
            int p = 1;
            for( auto run : rwgtRuns ){
                rwgtProcs.push_back( run.comRunProc() );
                if( run.rwgtName == "" ){
                    rwgtNames.push_back( "rwgt_" + std::to_string( p++ ) );
                } else {
                    rwgtNames.push_back( std::string(run.rwgtName) );
                }
            }
        }
        rwgtCard( std::string_view reweight_card ){
            srcCard = reweight_card;
        }
        rwgtCard( std::string_view reweight_card, REX::lesHouchesCard slhaParams, bool parseOnline = false ){
            srcCard = reweight_card;
            slhaCard = slhaParams;
            if( parseOnline ){ parse( parseOnline ); }
        }
        std::vector<std::shared_ptr<REX::lesHouchesCard>> writeCards( REX::lesHouchesCard& slhaOrig ){
            std::vector<std::shared_ptr<REX::lesHouchesCard>> cardVec;
            slhaOrig.parse();
            cardVec.reserve( rwgtRuns.size() );
            for( auto rwgt : rwgtRuns )
            {
                cardVec.push_back( rwgt.outWrite( slhaOrig ) );
            }
            return cardVec;
        }
    };

    struct rwgtCollection {
    public:
        void setRwgt( std::shared_ptr<rwgtCard> rwgts ){ 
            if( rwgtSet ){ return; }
            rwgtSets = rwgts; 
            rwgtSet = true;
        }
        void setRwgt( rwgtCard rwgts ){ 
            if( rwgtSet ){ return; }
            setRwgt( std::make_shared<rwgtCard>( rwgts ) ); rwgtSet = true; 
        }
        void setSlha( std::shared_ptr<REX::lesHouchesCard> slha ){ 
            if( slhaSet ){ return; }
            slhaParameters = slha; 
            slhaParameters->parse(); 
            slhaSet = true; 
        }
        void setSlha( REX::lesHouchesCard slha ){ 
            if( slhaSet ){ return; }
            setSlha( std::make_shared<REX::lesHouchesCard>( slha ) ); 
            slhaSet = true;
        }
        void setLhe( std::shared_ptr<REX::lheNode> lhe ){ 
            if( lheFileSet ){ return; }
            lheFile = lhe;
            lheFileSet = true;
        }
        void setLhe( REX::lheNode lhe ){ 
            if( lheFileSet ){ return; } 
            setLhe( std::make_shared<REX::lheNode>( lhe ) ); 
            lheFileSet = true;
        }
        void setLhe( std::string_view lhe_file ){
            if( lheFileSet ){ return; }
            size_t strt = 0;
            size_t post = *REX::nodeEndFind( lhe_file, strt );
            lheFile = REX::lheParser( lhe_file, strt, post );
            lheFileSet = true;
        }
        std::shared_ptr<rwgtCard> getRwgt(){ return rwgtSets; }
        std::shared_ptr<REX::lesHouchesCard> getSlha(){ return slhaParameters; }
        std::shared_ptr<REX::lheNode> getLhe(){ return lheFile; }
        rwgtCollection(){ return; }
        rwgtCollection( std::shared_ptr<REX::lheNode> lhe, std::shared_ptr<REX::lesHouchesCard> slha, std::shared_ptr<rwgtCard> rwgts ){
            setLhe( lhe );
            setSlha( slha );
            setRwgt( rwgts );
        }
    protected:
        void setDoubles(){
            if( lheFile == nullptr || rwgtSets == nullptr || slhaParameters == nullptr )
                throw std::runtime_error( "One or more of the necessary files (SLHA parameter card, LHE event storage file, and MadGraph-format reweight card) have not been initialised." );
            REX::lheRetDs returnBools; returnBools.xwgtup = true; returnBools.aqcdup = true; returnBools.pup = true;
            auto vecOfVecs = REX::lheValDoubles( *lheFile, returnBools );
            if( vecOfVecs->size() != 3 )
                throw std::runtime_error( "LHE file appears to contain multiple types of processes. This has not yet been implemented." );
            wgts = vecOfVecs->at( 0 ); gS = vecOfVecs->at( 1 ); momenta = vecOfVecs->at( 2 );
        }
        std::shared_ptr<rwgtCard> rwgtSets;
        std::shared_ptr<REX::lesHouchesCard> slhaParameters;
        std::shared_ptr<REX::lheNode> lheFile;
        std::shared_ptr<std::vector<double>> wgts;
        std::shared_ptr<std::vector<double>> gS;
        std::shared_ptr<std::vector<double>> momenta;
        bool lheFileSet = false;
        bool slhaSet = false;
        bool rwgtSet = false;
    };

    struct rwgtFiles : rwgtCollection {
        void setRwgtPath( std::string_view path ){ rwgtPath = path; }
        void setSlhaPath( std::string_view path ){ slhaPath = path; }
        void setLhePath( std::string_view path ){ lhePath = path; }
        rwgtFiles() : rwgtCollection(){ return; }
        rwgtFiles( std::string_view lhe_card, std::string_view slha_card, std::string_view reweight_card ) : rwgtCollection(){
            setRwgtPath( reweight_card );
            setSlhaPath( slha_card );
            setLhePath( lhe_card );
        }
        void initCards(){
            if( rwgtPath == "" || slhaPath == "" || lhePath == "" )
                throw std::runtime_error( "Paths to reweight card, parameter card, or LHE file have not been set" );
            pullRwgt(); pullSlha(); pullLhe();
            setLhe( *lheCard );
            setSlha( std::make_shared<REX::lesHouchesCard>( *slhaCard ) );
            setRwgt( std::make_shared<rwgtCard>( *rewgtCard, *slhaParameters, true ) );
            setDoubles();
        }
        void initCards( std::string_view lhe_card, std::string_view slha_card, std::string_view reweight_card ){
            setLhePath( lhe_card );
            setSlhaPath( slha_card );
            setRwgtPath( reweight_card );
            initCards();
        }
    protected:
        void pullRwgt(){
            rewgtCard = REX::filePuller( rwgtPath );
        }
        void pullSlha(){
            slhaCard = REX::filePuller( slhaPath );
        }
        void pullLhe(){
            lheCard = REX::filePuller( lhePath );
        }
        std::string rwgtPath;
        std::string lhePath;
        std::string slhaPath;
        std::shared_ptr<std::string> lheCard;
        std::shared_ptr<std::string> slhaCard;
        std::shared_ptr<std::string> rewgtCard;
    };

    struct rwgtRunner : rwgtFiles{
    public:
        void setMeEval( std::function<std::shared_ptr<std::vector<double>>(std::vector<double>&, std::vector<double>&)> eval ){ meEval = eval; meInit = true; }
        rwgtRunner() : rwgtFiles(){ return; }
        rwgtRunner( rwgtFiles& rwgts ) : rwgtFiles( rwgts ){ return; }
        rwgtRunner( rwgtFiles& rwgts, std::function<std::shared_ptr<std::vector<double>>(std::vector<double>&, std::vector<double>&)> meCalc ) : rwgtFiles( rwgts ){
            meEval = meCalc;
            meInit = true;
        }
        rwgtRunner( std::string_view lhe_card, std::string_view slha_card, std::string_view reweight_card,
        std::function<std::shared_ptr<std::vector<double>>(std::vector<double>&, std::vector<double>&)> meCalc ) : rwgtFiles( lhe_card, slha_card, reweight_card ){
            meEval = meCalc;
            meInit = true;
        }
    protected:
        bool meInit = false;
        bool meSet = false;
        bool normWgtSet = false;
        std::function<std::shared_ptr<std::vector<double>>(std::vector<double>&, std::vector<double>&)> meEval;
        std::shared_ptr<std::vector<double>> initMEs;
        std::shared_ptr<std::vector<double>> meNormWgts;
        std::shared_ptr<REX::weightGroup> rwgtGroup;
        void setMEs(){
            initCards(); 
            if( !meInit )
                throw std::runtime_error( "No function for evaluating scattering amplitudes has been provided." );
            auto ins = meEval( *momenta, *gS );
            initMEs = std::make_shared<std::vector<double>>( ins->begin(), ins->begin() + wgts->size() );
            meSet = true;
        }
        bool setParamCard( std::shared_ptr<REX::lesHouchesCard> slhaParams ){
            if( slhaPath == "" )
                throw std::runtime_error( "No parameter card path has been provided." );
            if( slhaParameters == nullptr )
                throw std::runtime_error( "No SLHA parameter card has been provided." );
            if( !REX::filePusher( slhaPath, *slhaParams->selfWrite() ) )
                throw std::runtime_error( "Failed to overwrite parameter card." );
            return true;
        }
        void setNormWgts(){
            if( !meSet ){ setMEs(); } 
            if( initMEs->size() != wgts->size() )
                throw std::runtime_error( "Inconsistent number of events and event weights." );
            meNormWgts = std::make_shared<std::vector<double>>( wgts->size() );
            for( size_t k = 0; k < initMEs->size(); k++ ){
                meNormWgts->at( k ) = wgts->at( k ) / initMEs->at( k );
            }
            normWgtSet = true;
        }
        bool singleRwgtIter( std::shared_ptr<REX::lesHouchesCard> slhaParams, std::shared_ptr<REX::lheNode> lheFile, size_t currId ){
            if( !normWgtSet )
                throw std::runtime_error( "Normalised original weights (wgt/|ME|) not evaluated -- new weights cannot be calculated." );
            if( !setParamCard( slhaParams ) )
                throw std::runtime_error( "Failed to rewrite parameter card." );
            auto newMEs = meEval( *momenta, *gS );
            auto newWGTs = REX::vecElemMult( *newMEs, *meNormWgts );
            REX::newWgt nuWgt( rwgtSets->rwgtRuns[currId].comRunProc(), newWGTs );
            lheFile->addWgt( 0, nuWgt );
            return true;
        }
        bool singleRwgtIter( std::shared_ptr<REX::lesHouchesCard> slhaParams, std::shared_ptr<REX::lheNode> lheFile, size_t currId, std::string& id ){
            if( !normWgtSet )
                throw std::runtime_error( "Normalised original weights (wgt/|ME|) not evaluated -- new weights cannot be calculated." );
            if( !setParamCard( slhaParams ) )
                throw std::runtime_error( "Failed to rewrite parameter card." );
            auto newMEs = meEval( *momenta, *gS );
            auto newWGTs = REX::vecElemMult( *newMEs, *meNormWgts );
            REX::newWgt nuWgt( rwgtSets->rwgtRuns[currId].comRunProc(), newWGTs, id );
            lheFile->addWgt( 0, nuWgt );
            return true;
        }
        bool lheFileWriter( std::shared_ptr<REX::lheNode> lheFile, std::string outputDir = "rwgt_evts.lhe" ){
            bool writeSuccess = REX::filePusher( outputDir, *lheFile->nodeWriter() );
            if( !writeSuccess )
                throw std::runtime_error( "Failed to write LHE file." );
            return true;
        }
    public:
        void runRwgt( const std::string& output ){
            setMEs();
            setNormWgts();
            rwgtGroup = std::make_shared<REX::weightGroup>();
            auto currInd = lheFile->header->addWgtGroup( rwgtGroup );
            auto paramSets = rwgtSets->writeCards( *slhaParameters );
            for( int k = 0 ; k < paramSets.size(); k++ ){
                singleRwgtIter( paramSets[k], lheFile, k, rwgtSets->rwgtNames[k] );
                std::cout << ".";
            }
            lheFileWriter( lheFile, output );
            REX::filePusher( slhaPath, *slhaCard );
            std::cout << "\nReweighting done.\n";
        }
    };
}