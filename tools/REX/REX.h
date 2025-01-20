/***
 *    ______          
 *    | ___ \         
 *    | |_/ /_____  __
 *    |    // _ \ \/ /
 *    | |\ \  __/>  < 
 *    \_| \_\___/_/\_\
 *                                             
 ***/
//
// *R*apid *e*vent e*x*traction Version 0.9.0
// Rex is a C++ library for parsing and manipulating Les Houches Event-format (LHE) files.
// It is designed to fast and lightweight, in comparison to internal parsers in programs like MadGraph.
// Currently, Rex is in development and may not contain all features necessary for full LHE parsing;
// particularly, it can only parse existing LHE files, rather than writing completely new ones.
//
// Copyright © 2023-2024 CERN, CERN Author Zenny Wettersten. 
// Licensed under the GNU Lesser General Public License (version 3 or later).
// All rights not expressly granted are reserved.
//

#ifndef _REX_H_
#define _REX_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>  
#include <string>
#include <set>
#include <cmath>
#include <utility>
#include <memory>
#include <map>
#include <algorithm>
#include <cctype>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <queue>
#include <charconv>

// ZW: all fcns within the REX standard sit in the
// namespace REX
// Note that as a convention, std::string_view objects will be
// referred to as strings unless the difference is relevant
namespace REX
{
    //#pragma warning( push )
    //#pragma warning( disable : 4101)
    static const size_t npos = (size_t)-1;
    #define UNUSED(x) (void)(x)
    //#pragma warning( pop ) 

    using sortFcn = std::function<std::shared_ptr<std::vector<size_t>>(std::vector<int>)>;
    using statSort = std::function<std::shared_ptr<std::vector<size_t>>(int, std::vector<int>)>;

    // ZW: generic warning function for printing warnings without throwing anything
    void warning( std::string message );

    // ZW: generic fcns for converting string-like objects to integers and doubles
    // Assumes input has no leading blankspace to check for a leading +,
    // but should not fail as long as there is no + even if there is blankspace
    // Note that Str needs to have a .compare(), .data() and .size() method
    template <typename Str>
    int ctoi(  Str str );
    extern template int ctoi<std::string>( std::string str );
    extern template int ctoi<std::string_view>( std::string_view str );
    template <typename Str>
    double ctod( Str str );
    extern template double ctod<std::string>( std::string str );
    extern template double ctod<std::string_view>( std::string_view str );

    // ZW: index sorting function, which returns vector
    // of the indices of the original vector sorted 
    // by default in ascending order
    // ie, for [5.0, 0.25, 2.0, 9.2] returns [1, 2, 0, 3]
    template <typename T>
    std::shared_ptr<std::vector<size_t>> indSort(const std::vector<T> &vector, std::function<bool(const T&, const T&)> comp = std::less<T>());
    extern template std::shared_ptr<std::vector<size_t>> indSort<int>(const std::vector<int> &vector, std::function<bool(const int&, const int&)> comp = std::less<int>());
    extern template std::shared_ptr<std::vector<size_t>> indSort<double>(const std::vector<double> &vector, std::function<bool(const double&, const double&)> comp = std::less<double>());

    template <typename T>
    std::shared_ptr<std::vector<size_t>> stoiSort(const std::vector<T> &vector);
    extern template std::shared_ptr<std::vector<size_t>> stoiSort<std::string_view>(const std::vector<std::string_view> &vector);

    template <typename T>
    std::shared_ptr<std::vector<size_t>> getRefOrder(const std::vector<T>& reference, const std::vector<T>& to_sort);
    extern template std::shared_ptr<std::vector<size_t>> getRefOrder<std::string_view>(const std::vector<std::string_view>& reference, const std::vector<std::string_view>& to_sort);

    std::shared_ptr<std::vector<size_t>> findEach( std::string_view textFile, std::string_view searchTerm );
    std::shared_ptr<std::vector<std::string_view>> lineSplitter( std::string_view currEvt );
    std::shared_ptr<std::vector<std::string_view>> blankSplitter( std::string_view currEvt );
    std::shared_ptr<std::string> filePuller( const std::string& fileLoc );
    bool filePusher( std::string fileLoc, std::string fileCont );

    // ZW: templated fcn for multiplying two vectors elementwise,
    // assuming T has a multiplication operator*
    template<typename T>
    std::shared_ptr<std::vector<T>> vecElemMult( const std::vector<T>& vec1, const std::vector<T>& vec2){
        if( vec1.size() < vec2.size() ){ return vecElemMult( vec2, vec1 ); }
        auto valVec = std::make_shared<std::vector<T>>( vec1.size() );
        std::transform( vec1.begin(), vec1.end(), vec2.begin(), valVec->begin(), []( const T& v1, const T& v2 ){
            return v1 * v2;
        } );
        return valVec;
    }
    struct xmlTree;

// ZW: struct for handling tags in XML node opening tags
    struct xmlTag {
    public:
        void setVal( std::string_view valSet );
        void setId( std::string_view idSet );
        std::string_view getVal();
        std::string_view getId();
        bool isModded();
        xmlTag();
        xmlTag( xmlTag& oldTag );
        xmlTag( std::string_view initId, std::string_view initVal);
    protected:
        bool modded;
        std::string_view val;
        std::string_view id;
    };

    struct xmlTree{
    public:
        xmlTree();
        xmlTree( std::string_view file );
        xmlTree( std::string_view file, size_t& strt, size_t& nd );
        auto getChildren(){ return children; }
        std::string_view& getOrigin(){ return origin; }
        size_t getStart(){ return start; }
        size_t getEnd(){ return end; }
        size_t getContStart(){ return contSt; }
        size_t getContEnd(){ return contEnd; }
        std::pair<size_t, size_t> getCont(){ if (contSt == npos || contEnd == npos){ return {0,0}; } return {contSt, contEnd - contSt}; }
        std::pair<size_t, size_t> getHeader(){ // ZW: pointers to beginning and end of node header, with some safety checks
            if( start == npos || end == npos || faux ){ return {0,0}; }
            if( contSt == npos ){ return {start, end - start }; }
            return {start, contSt - start};
        }
        bool isFaux(){ return faux; }
        bool isInit(){ return initialised; }
        bool hasChildren(){ return children->size() > 0; }
    protected:
        std::shared_ptr<std::vector<std::shared_ptr<xmlTree>>> children; // vector of pointers to children nodes
        std::string_view origin;
        size_t start; // position of opening bracket of node opening
        size_t end; // position of final character of ending node, including trailing blankspace
        size_t contSt;
        size_t contEnd;
        bool faux = false; // boolean showing whether this item is a true node or content squeezed between nodes
        bool initialised;
    };

    struct xmlNode {
    public:
        xmlNode();
        xmlNode( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} );
        xmlNode( xmlTree &tree );
        xmlNode( const xmlNode& original );
        std::vector<std::shared_ptr<xmlNode>> getChildren();
        std::vector<std::shared_ptr<xmlTag>> getTags();
        std::string_view getFile();
        std::string_view getName();
        std::string_view getContent();
        size_t getStart();
        size_t getEnd();
        xmlTree getTree();
        virtual bool isModded();
        virtual bool isModded( bool deep );
        bool isWritten();
        bool isParsed();
        bool isFaux();
        bool hasChildren();
        void setModded( bool mod );
        bool deepModded();
        bool deepParse();
        void parser( bool recursive );
        void addChild( std::shared_ptr<xmlNode> child );
        void addTag( std::shared_ptr<xmlTag> tag );
        void setFile( std::string_view file );
        void setName( std::string_view newName );
        void setCont( std::string_view cont );
    protected:
        virtual bool parse();
        virtual bool parse( bool recurs );
        bool parseTop();
        virtual bool parseContent();
        bool parseChildren( bool recursive );
        std::string nodeHeader;
        std::string nodeContent;
        std::string nodeEnd;
        xmlTree structure;
        std::vector<std::shared_ptr<xmlNode>> children;
        std::vector<std::shared_ptr<xmlTag>> tags;
        std::shared_ptr<std::string> writtenSelf; 
        bool deepMod = false;
        std::string_view xmlFile;
        std::string_view name;
        std::string_view content;
        size_t start;
        size_t end = npos;
        bool modded = false;
        bool written = false;
        bool parsed = false;
        bool deepParsed = false;
        bool faux = false;
        virtual void headWriter();
        virtual void endWriter();
        virtual void contWriter();
        virtual void childWriter();
        virtual void endFinder();
        virtual void fullWriter();
    public:
        virtual int childCounter();
        virtual void childCounter( int& noChilds );
        virtual std::shared_ptr<std::string> nodeWriter();
    };

    struct lhePrt{
    public:
        std::string_view getLine();
        std::string_view getComment();
        std::vector<double> getMom();
        double getE();
        double getMass();
        double getVTim();
        double getSpin();
        int getPDG();
        int getStatus();
        std::vector<int> getMothers();
        std::vector<int> getColor();
        void setComment( std::string_view nuCom );
        void setMom( std::vector<double> nuMom );
        void setEnergy( double nuE );
        void setMass( double nuM );
        void setVTim( double nuVTim );
        void setSpin( double nuSpin );
        void setPDG( int nuPDG );
        void setStatus( int nuSt );
        void setMothers( std::vector<int> nuMum );
        void setColors( std::vector<int> nuCol );
        bool isModded();
        bool isWritten();
        std::shared_ptr<std::string> getContent();
        lhePrt();
        lhePrt( std::pair<int,int>& prtInfo );
        lhePrt( int idup, int istup, int moth1, int moth2, int icol1, int icol2, double px, double py, double pz, double e, double m, double vt, double sp );
        lhePrt( std::pair<std::string,std::string>& prtInfo );
        lhePrt( const std::string_view originFile, const size_t& beginLine = 0, const size_t& endLine = npos );
    protected:
        std::shared_ptr<std::string> content;
        std::string_view sourceFile;
        std::string_view comment;
        double mom[3];
        double energy;
        double mass;
        double vtim;
        double spin;
        int pdg;
        int status;
        int mothers[2];
        int icol[2];
        bool modded = false;
        bool written = false;
        void writer();
    };

    struct evHead {
    public:
        std::string_view getComment();
        double getWeight();
        double getScale();
        double getAQED();
        double getAQCD();
        int getNprt();
        int getProcID();
        bool isModded();
        bool isWritten();
        void setComment( std::string_view nuCom );
        void setWeight( double nuWgt );
        void setScale( double nuScale );
        void setAQED( double nuAQED );
        void setAQCD( double nuAQCD );
        void setNprt( int nuNprt );
        void setProcID( int nuProcID );
        std::shared_ptr<std::string> getContent();
        evHead();
        evHead( const std::string_view originFile, size_t beginLine = 0, size_t endLine = npos );
    protected:
        std::shared_ptr<std::string> content;
        std::string_view sourceFile;
        std::string_view comment;
        double weight;
        double scale;
        double aqed;
        double aqcd;
        int nprt;
        //int nprtint;
        //std::string nprtstr;
        int procid;
        bool modded = false;
        bool written = false;
        void writer();
    };

    struct bodyWgt : public xmlNode {
    public:
        void setComment( std::string_view nuComment );
        void setVal( std::string nuVal );
        void setVal( std::string_view nuVal );
        void setVal( double nuVal );
        void setId( std::string nuId );
        void setModded( bool nuModded );
        std::string_view getComment();
        std::string_view getValS();
        double getValD();
        bodyWgt();
        bodyWgt( std::string_view value );
        bodyWgt( double value );
        bodyWgt( std::string_view value, xmlTag rwgtId );
        bodyWgt( double value, xmlTag rwgtId );
        bodyWgt( std::string_view value, std::shared_ptr<xmlTag> rwgtId );
        bodyWgt( double value, std::shared_ptr<xmlTag> rwgtId );
        bodyWgt( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} );
        bodyWgt( xmlNode& wgtNode );
        bodyWgt( xmlNode* wgtNode );
        bodyWgt( std::shared_ptr<xmlNode> wgtNode );
        bodyWgt( xmlTree& wgtTree );
        bodyWgt( xmlTree* wgtTree );
        bodyWgt( std::shared_ptr<xmlTree> wgtTree );
        bodyWgt( double value, std::string& idTag );
        void appendWgt( std::shared_ptr<std::string> document );
        void appendWgt( std::string* document );
        std::shared_ptr<std::string> appendWgt( std::string_view document );
    protected:
        std::string_view comment;
        std::string valS;
        std::string id;
        double valD;
        void fullWriter() override;
    };

    struct event : public xmlNode {
    public:
        evHead getHead();
        std::vector<std::shared_ptr<lhePrt>> getPrts();
        std::map<int,std::vector<std::shared_ptr<lhePrt>> > getSortedPrts();
        std::vector<std::shared_ptr<bodyWgt>> getWgts();
        void setHead( evHead head );
        void addPrt( std::shared_ptr<lhePrt> prtcl );
        void addPrt( lhePrt prtcl );
        void setPrts( std::vector<std::shared_ptr<lhePrt>> prtcls );
        void addWgt( bodyWgt nuWgt );
        void addWgt( std::shared_ptr<bodyWgt> nuWgt );
        void addWgt( bodyWgt nuWgt, std::string& id );
        void addWgt( std::shared_ptr<bodyWgt> nuWgt, std::string& id );
        bool newWeight();
        int getNprt();
        bool isModded() override;
        bool isModded( bool deep ) override ;
        event();
        event( std::vector<std::pair<int,int>>& prtInfo );
        event( std::vector<std::pair<std::string,std::string>>& prtInfo );
        event( std::vector<std::shared_ptr<lhePrt>> prtInfo );
        event( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} ) ;
        event( const xmlNode& originFile );
        event( const xmlNode* originFile );
        event( const std::shared_ptr<xmlNode>& originFile );
        event( xmlTree& originFile );
        event( xmlTree* originFile );
        event( std::shared_ptr<xmlTree> originFile );
        event( const event& original );
        event( event* original );
        event( std::shared_ptr<event> original );
        bool prtsAreMod();
        bool headIsMod();
        bool isSpecSort() const;
        sortFcn getSortFcn() const;
        statSort getStatSort() const;
    protected:
        std::vector<std::shared_ptr<bodyWgt>> rwgt;
        std::shared_ptr<xmlNode> childRwgt;
        bool hasRwgt();
        bool rwgtChild();
        bool bothRwgt();
        bool eitherRwgt();
        evHead header;
        bool hasBeenProc = false;
        std::vector<std::shared_ptr<lhePrt>> prts;
        std::map<int, std::vector<int>> procMap;
        std::map<int, std::vector<size_t>> procOrder;
        std::map<int, std::vector<std::shared_ptr<lhePrt>> > sortPrts;
        sortFcn eventSort = []( std::vector<int> vec ){ return indSort( vec ); };
        statSort specSort = []( int stat, std::vector<int> vec ){ 
            UNUSED(stat); 
            return indSort( vec ); 
        };
        bool specSorted = false;
        bool initProcMap(bool hard = false);
        bool initProcMap( sortFcn sorter, bool hard = false );
        bool initProcMap( statSort sorter, bool hard = false );
        bool inRwgtChild( std::string_view name );
        bool checkRwgtOverlap();
        void childRwgtWriter();
        void vecRwgtWriter( bool midNode = false );
        void rwgtWriter();
        void contWriter() override;
        void childWriter() override;
        bool addedWgt = false;
        void fullWriter() override;
        void fullWriter( bool deep );
        void appendWgts();
    public:
        std::shared_ptr<std::string> nodeWriter() override;
        std::shared_ptr<std::string> nodeWriter( bool recursive );
        std::map<int, std::vector<int>> &getProc( bool hard = false );
        std::map<int, std::vector<size_t>> &getProcOrder( bool hard = false );
        std::map<int, std::vector<int>> getProc() const;
        std::map<int, std::vector<size_t>> getProcOrder() const;
        std::map<int, std::vector<int>> &getProc(sortFcn sorter, bool hard = true);
        std::map<int, std::vector<size_t>> &getProcOrder(sortFcn sorter, bool hard = true);
        std::map<int, std::vector<int>> &getProc(statSort sorter, bool hard = true);
        std::map<int, std::vector<size_t>> &getProcOrder(statSort sorter, bool hard = true);
    };

    using eventComparison = std::function<bool(event&, event&, std::vector<int>&)>;

    using eventSetComp = std::function<bool(event&, std::vector<int>&)>;

    struct eventSet{
        eventSet();
        eventSet( const eventSet& nuEvents );
        eventSet( std::vector<event>& nuEvents );
        eventSet( std::vector<std::shared_ptr<event>>& nuEvents );
        void setRelStats( std::vector<int>& nuStats );
        void addEvent( event& nuEvent );
        void addEvent( std::shared_ptr<event> nuEvent );
        void addEvent( std::vector<event>& nuEvents );
        void addEvent( std::vector<std::shared_ptr<event>> nuEvents );
        void setComp( eventSetComp nuComp );
        bool belongs( event& nuEvent );
        bool belongs( std::shared_ptr<event> nuEvent );
    protected:
        std::vector<event> events;
        std::vector<int> relStats = {-1, 1};
        eventSetComp comp;
    };

    struct paramVal{
    public:
        double value = 0;
        int id = 0;
        std::string_view realLine;
        std::string_view comment;
        std::string_view idStr;
        std::string_view valStr;
        virtual void parse();
        paramVal();
        paramVal( std::string_view paramLine, bool parseOnline = false );
        bool isMod();
        bool modded = false;
        virtual std::shared_ptr<std::string> selfWrite();
    };

    struct paramBlock {
    public:
        std::string_view realBlock;
        size_t startPt;
        std::string_view comment;
        std::string_view initComm;
        std::string_view name;
        std::vector<paramVal> params;
        virtual void parse( bool parseOnline = false );
        paramBlock();
        paramBlock( std::string_view paramSet, bool parseOnline = false );
        bool isMod();
        bool modded = false;
        virtual std::shared_ptr<std::string> selfWrite();
    };

    struct decVal : public paramVal{
    public:
        void parse() override;
        decVal( std::string_view paramLine = "", bool parseOnline = false );
        std::shared_ptr<std::string> selfWrite() override;
    };

    struct decBlock : public paramBlock {
    public:
        std::vector<decVal> decays;
        void parse( bool parseOnline = false ) override;
        void parse( std::shared_ptr<std::vector<size_t>> decLines, bool parseOnline = false );
        decBlock( std::string_view paramSet = "", bool parseOnline = false );
        std::shared_ptr<std::string> selfWrite() override;
    };

    bool clStringComp( std::string_view str1, std::string str2 );
    bool clStringComp( std::string_view str1, std::string_view str2 );
    bool clStringComp( std::string str1, std::string str2 );
    bool clStringComp( std::string str1, std::string_view str2 );

    struct lesHouchesCard {
    public:
        decBlock decays;
        std::string_view xmlFile;
        size_t start;
        size_t end;
        bool modded;
        bool parsed;
        std::string_view header;
        std::vector<paramBlock> blocks;
        size_t blockStart;
        std::function<bool(size_t&, const std::string_view&)> lambda = [&]( size_t& conPt, const std::string_view& file )
            { return !( file[conPt+1] == ' ' || file[conPt+1] == '#' || file[conPt+1] == '\n' ); };
        std::function<bool(size_t&, const std::string_view&)> lambdaNu = [&]( size_t& conPt, const std::string_view& file )
            { return !( file[conPt+1] == ' ' || file[conPt+1] == '\n' || file[conPt+1] == '<'); };
        std::function<bool(size_t&, const std::string_view&)> lambdaD = [&]( size_t& conPt, const std::string_view& file )
            { return !(  clStringComp(file.substr(conPt+1, 1), std::string("d") ) ); };
        void parse( bool parseOnline = false );
        lesHouchesCard( const std::string_view originFile = "", const size_t& begin = 0, bool parseOnline = false );
        bool isMod();
        std::shared_ptr<std::string> selfWrite();
    };


    struct headWeight : public xmlNode {
    public:
        int getId();
        std::string_view getTag();
        bool hasTag();
        headWeight();
        headWeight( std::string_view paramSet, const size_t& begin = 0 );
        headWeight( std::string_view paramSet, std::string_view idText, int idNo, const size_t& begin = 0 );
        headWeight( xmlNode& node );
        headWeight( xmlNode* node );
        headWeight( std::shared_ptr<xmlNode> node );
        headWeight( xmlTree& tree );
        headWeight( xmlTree* tree );
        headWeight( std::shared_ptr<xmlTree> tree );
        headWeight( std::string_view paramSet, std::string& idText, unsigned int idNo, const size_t& begin = 0 );
        headWeight( std::string_view paramSet, std::string& idText);
        void setId( std::string identity );
    protected:
        std::string idTag;
        long unsigned int id = npos;
        void headWriter() override;
        void headWriter( bool incId );
        void endWriter() override;
        void contWriter() override;
        void childWriter() override;
        void childWriter( bool hasChildren );
        void fullWriter() override;
        void fullWriter( bool incId, bool hasChildren=true );
    };


    // ZW: struct for handling rwgt groups
    // in the LHE header initrwgt node
    struct weightGroup : public xmlNode {
    public:
        bool getIncId();
        void setIncId( bool nuIncId );
        std::vector<std::shared_ptr<headWeight>> getWgts();
        void addWgt( headWeight nuWgt );
        void addWgt( std::shared_ptr<headWeight> nuWgt );
        weightGroup();
        weightGroup( std::vector<std::shared_ptr<headWeight>> nuWgts );
        weightGroup( std::vector<std::string> nuWgts );
        weightGroup( xmlNode& wgtNode );
        weightGroup( xmlNode* wgtNode );
        weightGroup( xmlTree& wgtTree );
        weightGroup( xmlTree* wgtTree );
        weightGroup( std::shared_ptr<xmlTree> wgtTree );
        weightGroup( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} );
    protected:
        std::string_view rwgtName;
        std::string_view wgtNamStrat;
        bool includeId = false;
        std::vector<std::shared_ptr<headWeight>> paramSets;
        bool nu;
        std::string_view idTag;
        int id;
        void headWriter() override;
        void contWriter() override;
        void childWriter() override;
        void childWriter( bool hasChildren );
        void endWriter() override;
    };


    struct initRwgt : public xmlNode {
    public:
        std::vector<std::shared_ptr<weightGroup>> getGroups();
        size_t noGrps();
        void addGroup( weightGroup nuGroup );
        void addGroup( std::shared_ptr<weightGroup> nuGroup );
        void addWgt( unsigned int index, std::shared_ptr<headWeight> nuWgt );
        void addWgt( unsigned int index, headWeight nuWgt );
        initRwgt();
        initRwgt( std::vector<std::shared_ptr<xmlNode>> nuGroups );
        initRwgt( xmlNode& wgtNode );
        initRwgt( xmlNode* wgtNode );
        initRwgt( std::shared_ptr<xmlNode> wgtNode );
        initRwgt( xmlTree& wgtTree );
    protected:
        bool grpIsInit = false;
        bool grpInit( std::shared_ptr<weightGroup>& wgt );
        std::vector<std::shared_ptr<weightGroup>> groups = {};
        void contWriter() override;
        void childWriter() override;
        void childWriter( bool hasChildren );
    };

    struct lheInitHead{
    public:
        std::string_view idbmup[2];
        std::string_view ebmup[2];
        std::string_view pdfgup[2];
        std::string_view pdfsup[2];
        std::string_view idwtup;
        std::string_view nprup;
        bool isWritten();
        bool isModded();
        std::shared_ptr<std::string> getContent();
        lheInitHead( std::string_view initHead );
        lheInitHead( xmlNode& initNode );
    protected:
        std::shared_ptr<std::string> content;
        bool written = false;
        bool modded = false;
        void writer();
    };

    struct lheInitLine {
    public:
        std::string_view xsecup;
        std::string_view xerrup;
        std::string_view xmaxup;
        std::string_view lprup;
        bool isWritten();
        bool isModded();
        std::shared_ptr<std::string> getContent();
        lheInitLine();
        lheInitLine( std::string_view procLine );
    protected:
        std::shared_ptr<std::string> content;
        bool written = false;
        bool modded = false;
        void writer();
    };


    struct slhaNode : public xmlNode {
    public:
        std::shared_ptr<lesHouchesCard> getParameters();
        slhaNode();
        slhaNode( lesHouchesCard parameters );
        slhaNode( std::shared_ptr<lesHouchesCard> parameters );
        slhaNode( xmlNode& node, bool parseOnline = false );
        slhaNode( xmlNode* node, bool parseOnline = false );
        slhaNode( std::shared_ptr<xmlNode> node, bool parseOnline = false );
        slhaNode( xmlTree tree, bool parseOnline = false );
        slhaNode( std::shared_ptr<xmlTree> tree, bool parseOnline = false );
        slhaNode( xmlTree* tree, bool parseOnline = false );
        slhaNode( const std::string_view originFile, const size_t& begin = 0, bool parseOnline = false );
    protected:
        std::shared_ptr<lesHouchesCard> parameterCard;
        bool pCardInit = false;
        void headWriter() override;
        void endWriter() override;
        void contWriter() override;
    };

    struct initNode : public xmlNode {
    public:
        std::shared_ptr<lheInitHead> getHead();
        std::vector<std::shared_ptr<lheInitLine>> getLines();
        void setHead( std::shared_ptr<lheInitHead> head );
        void setLines( std::vector<std::shared_ptr<lheInitLine>> lines );
        void addLine( std::shared_ptr<lheInitLine> line );
        initNode();
        initNode( const std::string_view originFile, const size_t& begin = 0, bool parseOnline = false );
        initNode( xmlNode& node, bool parseOnline = false );
        initNode( xmlNode* node, bool parseOnline = false );
        initNode( std::shared_ptr<xmlNode> node, bool parseOnline = false );
        initNode( xmlTree tree, bool parseOnline = false );
        initNode( std::shared_ptr<xmlTree> tree, bool parseOnline = false );
        initNode( xmlTree* tree, bool parseOnline = false );
    protected:
        std::shared_ptr<lheInitHead> initHead;
        std::vector<std::shared_ptr<lheInitLine>> initLines;
        bool parseContent() override;
        void contWriter() override;
    };

    struct lheHead : public xmlNode {
    public:
        size_t addWgtGroup( std::shared_ptr<weightGroup>& wgtGroup );
        size_t addWgtGroup( weightGroup wgtGroup );
        void addWgt( size_t index, std::shared_ptr<headWeight> nuWgt );
        void addWgt( size_t index, headWeight nuWgt );
        void addWgt( size_t index, std::shared_ptr<headWeight> nuWgt, std::string idTagg );
        void addWgt( size_t index, headWeight nuWgt, std::string idTagg );
        void setInitRwgt( initRwgt initWgt );
        void setInitRwgt( std::shared_ptr<initRwgt> initWgt );
        std::vector<std::shared_ptr<weightGroup>> getWgtGroups();
        std::shared_ptr<initRwgt> getInitRwgt();
        std::shared_ptr<slhaNode> getParameters();
        void setParameters( std::shared_ptr<slhaNode> params );
        bool rwgtInc();
        lheHead();
        lheHead( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} );
        lheHead( xmlNode& node );
        lheHead( xmlNode* node );
        lheHead( std::shared_ptr<xmlNode> node );
        lheHead( xmlTree tree );
        lheHead( std::shared_ptr<xmlTree> tree );
        lheHead( xmlTree* tree );
    protected:
        bool wgtGrpIsInit = false;
        bool wgtGrpInit( std::shared_ptr<weightGroup>& wgtGrp );
        std::shared_ptr<slhaNode> parameters;
        bool hasRwgt = false;
        std::shared_ptr<initRwgt> rwgtNodes;
        std::vector<std::shared_ptr<weightGroup>> initrwgt;
        bool relChildSet = false;
        std::vector<int> relChild;
        void setRelChild();
        bool parseChildren( bool recursive );
        void headWriter() override;
        void childWriter() override;
        void fullWriter() override;
    };

    struct newWgt{
    protected:
        std::shared_ptr<headWeight> headWgt;
        std::vector<std::shared_ptr<bodyWgt>> bodyWgts;
    public:
        newWgt( std::shared_ptr<headWeight> heaWgt, std::vector<std::shared_ptr<bodyWgt>> bodWgts );
        newWgt( std::shared_ptr<headWeight> heaWgt, std::shared_ptr<std::vector<double>> wgts );
        newWgt( std::string_view parameters, std::shared_ptr<std::vector<double>> wgts, std::string idTag = "rex_rwgt" );
        newWgt( std::string_view parameters, int idNum, std::shared_ptr<std::vector<double>> wgts, std::string idTag = "rex_rwgt" );
        newWgt( std::string& parameters );
        newWgt( std::string& parameters, std::string& idTag );
        std::shared_ptr<headWeight> getHeadWgt();
        std::vector<std::shared_ptr<bodyWgt>> getBodyWgts();
        void addBdyWgts( std::shared_ptr<std::vector<double>> wgts );
    };


    struct lheNode : public xmlNode {
    public:
        lheNode();
        lheNode( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} );
        std::shared_ptr<lheHead> getHeader();
        std::shared_ptr<initNode> getInit();
        std::vector<std::shared_ptr<event>> getEvents();
        bool isModded() override;
        bool isModded( bool deep ) override;
        void setInit( std::shared_ptr<initNode> initNod );
        void setHeader( std::shared_ptr<lheHead> headNod );
        void addWgt( size_t index, newWgt& addedWgt );
        void addWgt( size_t index, newWgt& addedWgt, std::string& idTag );
        void setRelStats( std::vector<std::string_view>& particles );
        std::vector<std::string_view>& getRelStats();
        void setSameSort( sortFcn& sortF );
        sortFcn& getSameSort();
        void setStatSort( statSort& statS );
        statSort& getStatSort();
    protected:
        std::vector<std::shared_ptr<event>> events = {};
        std::shared_ptr<lheHead> header =  std::make_shared<lheHead>(xmlFile, start);
        std::shared_ptr<initNode> init = std::make_shared<initNode>(xmlFile, start);
        std::vector<std::string_view> relStat = {"-1", "1"};
        sortFcn particleSort = []( std::vector<int> prts ){ return indSort(prts); };
        statSort statParticleSort = []( int dummy, std::vector<int> prts ){ UNUSED(dummy); return indSort(prts); };
        virtual void headerWriter();
        virtual void initWriter();
        virtual void eventWriter();
        void contWriter() override;
        void fullWriter() override;
    public:    
        virtual std::shared_ptr<std::string> nodeWriter();
    };

    struct evtInfo {
    public:
        std::vector<double> wgts;
        std::vector<double> scales;
        std::vector<double> aQEDs;
        std::vector<double> aQCDs;
        std::vector<int> nprts;
        std::vector<size_t> relNPrts;
        std::vector<int> procIDs;
        evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile = {} );
        evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const std::vector<int>& statVec );
        evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const std::vector<int>& statVec, 
        sortFcn sorter );
        evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const std::vector<int>& statVec, 
        statSort sorter );
    };

    struct  prtInfo {
    public:
        std::vector<double> moms;
        std::vector<double> masses;
        std::vector<double> vtims;
        std::vector<double> spins;
        std::vector<int> statuses;
        std::vector<int> mothers;
        std::vector<int> icols;
        std::vector<int> pdgs;
        prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile = {}, const int nPrt = 8 );
        prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, const std::vector<int>& statVec );
        prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, const std::vector<int>& statVec, 
        sortFcn sorter );
        prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, const std::vector<int>& statVec, 
        statSort sorter );
    };

    struct transSkel {
    public:
        std::vector<std::vector<std::shared_ptr<REX::event>>> procSets;
        std::vector<std::shared_ptr<std::vector<bool>>> relProcs;
        std::vector<bool> relEvSet;
        transSkel();
        transSkel( transSkel& skeleton );
        transSkel( lheNode& lheFile, std::vector<eventSet>& evSet );
        transSkel( std::shared_ptr<lheNode> lheFile, std::vector<eventSet>& evSet );
    };

    struct transMonoLHE {
    public:
        evtInfo evtsHead;
        prtInfo evtsData;
        std::shared_ptr<event> process;
        transMonoLHE( const std::vector<std::shared_ptr<REX::event>> lheFile = {}, const int nPrt = 8 );
        transMonoLHE( const std::vector<std::shared_ptr<REX::event>> lheFile, const int nPrt, const std::vector<int>& statVec );
        transMonoLHE( const std::vector<std::shared_ptr<REX::event>> lheFile, const int nPrt, 
        sortFcn sorter,
        std::vector<int> statVec = { -1, 1 } );
        transMonoLHE( const std::vector<std::shared_ptr<REX::event>> lheFile, const int nPrt, 
        statSort sorter,
        std::vector<int> statVec = { -1, 1 } );
        transMonoLHE( const transMonoLHE& lheFile );
    };

    struct transLHE {
    public:
        std::string_view xmlFile;
        std::vector<std::shared_ptr<transMonoLHE>> subProcs;
        std::vector<std::shared_ptr<event>> procSets;
        std::vector<std::shared_ptr<std::vector<bool>>> relProcs;
        std::vector<bool> relEvSets;
        void setRelEvSets();
        transLHE();
        transLHE( lheNode& lheFile );
        transLHE( lheNode& lheFile, 
        sortFcn sorter, 
        const std::vector<int>& statVec = { -1, 1 } );
        transLHE( lheNode& lheFile, 
        statSort sorter, 
        const std::vector<int>& statVec = { -1, 1 } );
        transLHE( lheNode& lheFile, const std::vector<int>& statVec );
        transLHE( transSkel& skeleton );
        transLHE( const transLHE& lheFile );
        std::shared_ptr<std::vector<double>> vectorFlat( std::vector<std::shared_ptr<std::vector<double>>> vecVec );
    };

    struct lheRetDs{
    public:
        bool ebmup = false;
        bool xsecup = false;
        bool xerrup = false;
        bool xmaxup = false;
        bool xwgtup = false;
        bool scalup = false;
        bool aqedup = false;
        bool aqcdup = false;
        bool pup = true;
        bool mass = false;
        bool vtimup = false;
        bool spinup = false;
        std::vector<bool> getBools();
    };

    // ZW: bool struct to define which int values
    // to extract transposed from LHE file
    struct lheRetInts{
    public:
        //bool maxpup = false;
        bool idbmup = false;
        bool pdfgup = false;
        bool pdfsup = false;
        bool idwtup = false;
        bool nprup = false;
        bool lprup = false;
        //bool maxnup = false;
        bool nup = true;
        bool idprup = false;
        bool idup = true;
        bool istup = true;
        bool mothup = false;
        bool icolup = false;
        std::vector<bool> getBools();
    };

    struct eventComp{
        bool operator()( event& firstEv, event& secEv);
        bool operator()( const event& firstEv, const event& secEv) const;
        bool operator()(event& firstEv, event& secEv, std::vector<int> statVec);
    };

    using evComp = std::function<bool(event&, event&)>;
    using evCheck = std::function<bool(event&)>;

    // ZW: struct to define which values to extract from events in
    // an AoS LHE structure to an SoA structure
    // Note: has two separate args for momenta,
    // momUp denotes the physical 4-momentum (E, px, py, pz),
    // whereas pUp denotes the SLHA/LHE style lab frame momenta (px, py, pz, E, m)
    struct relEvArgs{
        bool nUp = true; // number of particles in this event
        bool idPrUp = true; // process ID
        bool xWgtUp = true; // event weight
        bool scalUp = false; // event scale
        bool aQEDUp = false; // alpha QED
        bool aQCDUp = true; // alpha QCD
        bool idUp = true; // particle PDG code
        bool iStUp = true; // particle status
        bool mothUp[2] = {false, false}; // particle mothers
        bool iColUp[2] = {false, false}; // particle colour flow
        bool massUp = false; // particle mass in GeV
        bool momUp = true; // generic flag for particle momenta; momUp denotes full 4-momentum (excl mass), following flags denote individual components
        bool pUp[5] = {false, false, false, false, false}; // lab frame momenta (px, py, pz, E, m) --- separates the different components
        bool vTimUp = false; // invariant lifetime in mm
        bool spinUp = false; // cosine of angle between spin and 3-momentum of particle in lab frame
        relEvArgs();
        relEvArgs( const relEvArgs& relData );
        // ZW: self-returning setters for each flag,
        // to allow for easy chaining of multiple flags
        // in a single line (ie relEvArgs().setNUp(true).setIdPrUp(true) etc)
        relEvArgs& setNUp( bool nuNUp );
        relEvArgs& setIdPrUp( bool nuIdPrUp );
        relEvArgs& setXWgtUp( bool nuXWgtUp );
        relEvArgs& setScalUp( bool nuScalUp );
        relEvArgs& setAQEDUp( bool nuAQEDUp );
        relEvArgs& setAQCDUp( bool nuAQCDUp );
        relEvArgs& setIdUp( bool nuIdUp );
        relEvArgs& setIStUp( bool nuIStUp );
        relEvArgs& setMothUp( bool nuMothUp[2] );
        relEvArgs& setMothUp( bool nuMothUp );
        relEvArgs& setMothUp( std::vector<bool> nuMothUp );
        relEvArgs& setIColUp( bool nuIColUp[2] );
        relEvArgs& setIColUp( bool nuIColUp );
        relEvArgs& setIColUp( std::vector<bool> nuIColUp );
        relEvArgs& setMassUp( bool nuMassUp );
        relEvArgs& setMomUp( bool nuMomUp );
        relEvArgs& setPUp( bool nuPUp[5] );
        relEvArgs& setPUp( bool nuPUp );
        relEvArgs& setPUp( std::vector<bool> nuPUp );
        relEvArgs& setVTimUp( bool nuVTimUp );
        relEvArgs& setSpinUp( bool nuSpinUp );
        relEvArgs& setAll( bool nuAll );
    };

    // ZW: individual process lines
    // keeping object oriented format due to size, 
    // but could easily be replaced with an SoA format if needed
    // (turn arguments here into vectors, and replace the lheSoA->procLines vector with a single struct)
    struct procLine{
        double xSecUp; // cross section in pb
        double xErrUp; // cross section error in pb
        double xMaxUp; // maximum event weight
        int lPrUp; // process ID label (corresponds to idPrUp in the events)
        procLine();
        procLine( const procLine& process );
        procLine( double xSec, double xErr, double xMax, int lPr );
        procLine( lheInitLine& process );
        procLine( std::shared_ptr<lheInitLine> process );
    };

    std::vector<procLine> lheProcLines( lheNode& lheFile );

    struct procSoA{
        relEvArgs relData;
        std::vector<int> relStats;
        evCheck relevantEvent;
        std::vector<std::shared_ptr<event>> events;
        std::vector<bool> relEvMap;
        std::vector<std::shared_ptr<event>> relEvents;
        std::vector<int> nUp;
        std::vector<int> idPrUp;
        std::vector<double> xWgtUp;
        std::vector<double> scalUp;
        std::vector<double> aQEDUp;
        std::vector<double> aQCDUp;
        std::vector<int> idUp;
        std::vector<int> iStUp;
        std::vector<std::vector<int>> mothUp;
        std::vector<std::vector<int>> iColUp;
        std::vector<double> massUp;
        std::vector<double> momUp;
        std::vector<std::vector<double>> pUp;
        std::vector<double> vTimUp;
        std::vector<double> spinUp;
        virtual void reset();
        void setRelEvs( std::vector<std::shared_ptr<event>> lheFile );
        bool relevant( event& ev );
        bool relevant( std::shared_ptr<event> ev );
        bool extract( std::vector<std::shared_ptr<event>> lheFile, std::vector<int> relevStats = {-1,1} );
        bool empty();
        procSoA();
        procSoA( const relEvArgs& relData );
        procSoA( const procSoA& process );
        procSoA( procSoA* process );
        procSoA( std::shared_ptr<procSoA> process );
        procSoA( std::vector<std::shared_ptr<event>> lheFile, relEvArgs relArgs = relEvArgs(), 
            std::vector<int> relevStats = {-1,1}, std::function<bool( event& )> relFcn = nullptr );
        procSoA& setRelData( relEvArgs& relArgs );
        procSoA& setRelStats( std::vector<int>& relevStats );
        procSoA& setRelevant( evCheck& relFcn );
        procSoA& setEvents( std::vector<std::shared_ptr<event>> lheFile );
    };

    struct lheSoA{
        std::vector<std::shared_ptr<procSoA>> subProcesses; // transposed events grouped by arbitrary sorting functions
        // ZW: note, at extraction events which fail all relEvFcns not extracted but are still stored in the events vectors
        std::vector<std::shared_ptr<event>> events; // all events in the LHE file
        std::vector<std::vector<std::shared_ptr<event>>> sortedEvents; // events grouped by process before extraction. events that fail all relEvFcns are stored in the last entry
        std::vector<procLine> procLines; // individual process lines from the LHE file header
        std::vector<evCheck> relEvFcns; // vector of event classifiers
        std::vector<size_t> eventGrouping; // indices to which subProcess each event belongs to
        std::vector<relEvArgs> relEvData; // vector of relEvArgs for each subProcess
        std::vector<std::vector<int>> relEvStats; // vector of relevant particle statuses for each subProcess
        int idBmUp[2]; // beam IDs
        double eBmUp[2]; // beam energies in GeV
        int pdfGUp[2]; // Cernlib PDFlib specification authors
        int pdfSUp[2]; // Cernlib PDFlib specification PDF sets
        int idWtUp; // event weight model
        int nPrUp; // number of different user subprocesses; note, not necessarily the same as subProcess vector
        virtual void reset(); // reset all data but maintain instantiation
        size_t eventIndex( event& ev ); // map an event to the correct subProcess, based on the order in relEvFcns
        bool sortEvents( bool hard = false ); // sort events based on relEvFcns into sortedEvents
        bool extractEvents( bool hard = false ); // extract events from sortedEvents into subProcesses
        std::function<bool( event& )> evFcnIndex( size_t index ); // return the event classifier function at index if it exists, otherwise return 0-index
        relEvArgs evDataIndex( size_t index ); // return the relEvArgs at index if it exists, otherwise return 0-index
        std::vector<int> evStatsIndex( size_t index ); // return the relevant particle statuses at index if it exists, otherwise return 0-index
        lheSoA();
        lheSoA( const lheSoA& lheFile );
        lheSoA( std::vector<std::shared_ptr<event>> lheFile );
        lheSoA( std::vector<std::shared_ptr<event>> lheFile, std::vector<std::function<bool( event& )>> evSort );
        lheSoA( lheNode& lheFile, std::vector<std::function<bool( event& )>> evSort );
        lheSoA( std::vector<std::shared_ptr<event>> lheFile, std::vector<std::function<bool( event& )>> evSort, 
            std::vector<relEvArgs> relData );
        lheSoA( lheNode& lheFile, std::vector<std::function<bool( event& )>> evSort, 
            std::vector<relEvArgs> relData );
        lheSoA( std::vector<std::shared_ptr<event>> lheFile, std::vector<std::function<bool( event& )>> evSort, 
            std::vector<std::vector<int>> relStats );
        lheSoA( lheNode& lheFile, std::vector<std::function<bool( event& )>> evSort,
            std::vector<std::vector<int>> relStats );
        lheSoA( std::vector<std::shared_ptr<event>> lheFile, std::vector<std::function<bool( event& )>> evSort, 
            std::vector<relEvArgs> relData, std::vector<std::vector<int>> relStats );
        lheSoA( lheNode& lheFile, std::vector<std::function<bool( event& )>> evSort,
            std::vector<relEvArgs> relData, std::vector<std::vector<int>> relStats );
        lheSoA& setInit( lheInitHead& init );
        lheSoA& setInit(  lheNode& lheFile );
        lheSoA& setEvents( std::vector<std::shared_ptr<event>> lheFile );
        lheSoA& setEvents( lheNode& lheFile );
        lheSoA& setProcLines( std::vector<procLine> procLines );
        lheSoA& setProcLines(  std::vector<std::shared_ptr<lheInitLine>> procLines );
        lheSoA& setProcLines( initNode& init );
        lheSoA& setProcLines( lheNode& lheFile );
        lheSoA& setRelEvFcns( std::vector<std::function<bool( event& )>> evSort );
        lheSoA& setRelEvData( std::vector<relEvArgs> relData );
    };

std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> lheValDoubles( lheNode& lheFile, lheRetDs vals = lheRetDs() );

std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> lheValDoubles(transLHE& lheAOS, lheRetDs vals = lheRetDs() );

}

#endif
