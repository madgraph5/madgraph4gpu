/***
 *      _____  ______ _____  _____  ______ _____  
 *     |  __ \|  ____|  __ \|  __ \|  ____|  __ \ 
 *     | |__) | |__  | |__) | |__) | |__  | |__) |
 *     |  ___/|  __| |  ___/|  ___/|  __| |  _  / 
 *     | |    | |____| |    | |    | |____| | \ \ 
 *     |_|    |______|_|    |_|    |______|_|  \_\                              
 *                                                
 ***/
#include <unistd.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <variant>
#include "PEP.hpp"

namespace PEP::PER
{
    template<typename T1, typename T2>
    std::shared_ptr<std::vector<T1>> scatAmpEval(std::vector<T2>& momenta, std::function<std::shared_ptr<std::vector<T1>>(std::vector<T2>&)> evalFunc)
    { return evalFunc(momenta); }

    template<typename T1, typename T2>
    std::shared_ptr<std::vector<T1>> scatAmpEval(std::vector<T2>& momenta, std::function<std::vector<T1>(std::vector<T2>&)> evalFunc)
    { return evalFunc(momenta); }

    struct rwgtVal : PEP::paramVal{
    public:
        std::string_view blockName;
        bool allStat;
        bool isAll(){ return (idStr == "all"); }
        rwgtVal( std::string_view paramLine = "" )
        : paramVal( paramLine, false ){if( paramLine.size() == 0 ){ return; }
            realLine = paramLine;
            auto vals = *PEP::nuBlankSplitter( realLine );
            blockName = vals[1];
            idStr = vals[2];
            valStr = vals[3];
        }
        void selfWrite( PEP::paramBlock& srcBlock )
        {
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
            [&]( const PEP::paramVal& parPar ){ return (parPar.idStr == idStr ); } );
            if( currPar == srcBlock.params.end() ){ 
                srcBlock.params.push_back( PEP::paramVal( realLine.substr(realLine.find("set") + 4) ) );
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
        void selfWrite( PEP::paramBlock& srcBlock, const std::map<std::string_view, int>& blocks )
        {
            for( auto parm : rwgtVals )
            {
                parm.selfWrite( srcBlock );
            }
            srcBlock.modded = true;
            return;
        }
    };

    struct rwgtProc {
    public:
        std::vector<rwgtBlock> rwgtParams;
        std::string_view procString;
        void parse(){
            std::vector<std::string_view> blocks;
            std::vector<std::shared_ptr<std::vector<rwgtVal>>> params;
            auto procLines = *PEP::nuLineSplitter( procString );
            for( auto line : procLines )
            {
                auto strtPt = line.find("set");
                auto words = *PEP::nuWordSplitter( line );
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
        rwgtProc( PEP::lesHouchesCard slhaSet, std::string_view rwgtSet = "", bool parseOnline = false )
        {
            if( rwgtSet == "" ){ return; }
            auto strtLi = rwgtSet.find( "\n", rwgtSet.find("launch") ) + 1;
            auto endLi = rwgtSet.find("\n", strtLi);
            while( rwgtSet[rwgtSet.find_first_not_of("\n ", endLi)] == 's' )
            { endLi = rwgtSet.find( "\n", endLi + 1 ); }
            procString = rwgtSet.substr( strtLi, endLi - strtLi );
            if( parseOnline ){ parse(); }
        }
        std::shared_ptr<PEP::lesHouchesCard> selfWrite( const PEP::lesHouchesCard& paramOrig ){
            auto slhaOrig = std::make_shared<PEP::lesHouchesCard>( paramOrig );
            std::map<std::string_view, int> blockIds;
            for( int k = 0 ; k < slhaOrig->blocks.size() ; ++k )
            { blockIds.insert({slhaOrig->blocks[k].name, k }); }
            for( auto rwgts : rwgtParams )
            { rwgts.selfWrite( slhaOrig->blocks[ blockIds.at( rwgts.name ) ], blockIds ); }
            slhaOrig->modded = true;
            return slhaOrig;
        }
    };

    struct rwgtCard{
    public:
        PEP::lesHouchesCard slhaCard;
        std::vector<rwgtProc> rwgtRuns;
        std::vector<std::string_view> opts;
        std::string_view srcCard;
        void parse( bool parseOnline = false ) {
            auto strt = srcCard.find("launch");
            while( srcCard.find_last_of("#", strt) > srcCard.find_last_of("\n", strt) ){ strt = srcCard.find("launch", strt + 6 ); }
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
            }
            size_t endLi = srcCard.find( "\n", lnchPos[ lnchPos.size() - 1 ] );
            if( srcCard.substr( endLi + 1, 3 ) == "set" ){
                while( srcCard.substr( endLi + 1, 3 ) == "set" )
                {
                    endLi = srcCard.find( "\n", endLi + 1 );
                }
                rwgtRuns.push_back( rwgtProc( slhaCard, srcCard.substr( lnchPos[lnchPos.size()-1], endLi - lnchPos[lnchPos.size()-1] ), parseOnline ) );
            }
        }
        rwgtCard( std::string_view reweight_card ){
            srcCard = reweight_card;
        }
        rwgtCard( std::string_view reweight_card, PEP::lesHouchesCard slhaParams, bool parseOnline = false ){
            srcCard = reweight_card;
            slhaCard = slhaParams;
            if( parseOnline ){ parse( parseOnline ); }
        }
        std::vector<std::shared_ptr<PEP::lesHouchesCard>> writeCards( PEP::lesHouchesCard& slhaOrig ){
            std::vector<std::shared_ptr<PEP::lesHouchesCard>> cardVec;
            slhaOrig.parse();
            cardVec.reserve( rwgtRuns.size() );
            for( auto rwgt : rwgtRuns )
            {
                cardVec.push_back( rwgt.selfWrite( slhaOrig ) );
            }
            return cardVec;
        }
    };

    struct rwgtCollection {
    public:
        void setRwgt( std::shared_ptr<rwgtCard> rwgts ){ rwgtSets = rwgts; }
        void setRwgt( rwgtCard rwgts ){ setRwgt( std::make_shared<rwgtCard>( rwgts ) ); }
        void setSlha( std::shared_ptr<PEP::lesHouchesCard> slha ){ slhaParameters = slha; }
        void setSlha( PEP::lesHouchesCard slha ){ setSlha( std::make_shared<PEP::lesHouchesCard>( slha ) ); }
        void setLhe( std::shared_ptr<PEP::lheNode> lhe ){ lheFile = lhe; }
        void setLhe( PEP::lheNode lhe ){ setLhe( std::make_shared<PEP::lheNode>( lhe ) ); }
        rwgtCollection(){ return; }
        rwgtCollection( std::shared_ptr<PEP::lheNode> lhe, std::shared_ptr<PEP::lesHouchesCard> slha, std::shared_ptr<rwgtCard> rwgts ){
            setLhe( lhe );
            setSlha( slha );
            setRwgt( rwgts );
        }
    protected:
        std::shared_ptr<rwgtCard> rwgtSets;
        std::shared_ptr<PEP::lesHouchesCard> slhaParameters;
        std::shared_ptr<PEP::lheNode> lheFile;
        std::shared_ptr<std::vector<double>> wgts;
        std::shared_ptr<std::vector<double>> momenta;
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
                throw std::runtime_error( "rwgtPath, slhaPath and lhePath must have not been set" );
            setLhe( std::make_shared<PEP::lheNode> ( *lheCard ) );
            setSlha( std::make_shared<PEP::lesHouchesCard>( *slhaCard ) );
            setRwgt( std::make_shared<rwgtCard>( *rewgtCard, *slhaParameters, true ) );
        }
        void initCards( std::string_view lhe_card, std::string_view slha_card, std::string_view reweight_card ){
            setLhePath( lhe_card );
            setSlhaPath( slha_card );
            setRwgtPath( reweight_card );
            setLhe( std::make_shared<PEP::lheNode> ( *lheCard ) );
            setSlha( std::make_shared<PEP::lesHouchesCard>( *slhaCard ) );
            setRwgt( std::make_shared<rwgtCard>( *rewgtCard, *slhaParameters, true ) );
        }
    protected:
        void pullRwgt(){
            rewgtCard = PEP::filePuller( rwgtPath );
        }
        void pullSlha(){
            rewgtCard = PEP::filePuller( slhaPath );
        }
        void pullLhe(){
            rewgtCard = PEP::filePuller( lhePath );
        }
        std::string rwgtPath;
        std::string lhePath;
        std::string slhaPath;
        std::shared_ptr<std::string> lheCard;
        std::shared_ptr<std::string> slhaCard;
        std::shared_ptr<std::string> rewgtCard;
    };

}