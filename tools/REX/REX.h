/***
 *    ______ _______   __
 *    | ___ \  ___\ \ / /
 *    | |_/ / |__  \ V / 
 *    |    /|  __| /   \ 
 *    | |\ \| |___/ /^\ \
 *    \_| \_\____/\/   \/
 *                                             
 ***/

// THIS IS NOT A LICENSED RELEASE
// IF YOU SEE THIS FILE, IT HAS BEEN SPREAD
// FROM AN IMPROPER RELEASE.

// Copyright © 2023-2024 CERN, CERN Author Zenny Wettersten. 
// All rights reserved.

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

    using sortFcn = std::function<std::shared_ptr<std::vector<size_t>>(std::vector<std::string_view>)>;
    using statSort = std::function<std::shared_ptr<std::vector<size_t>>(std::string_view, std::vector<std::string_view>)>;

    template <typename T>
    std::shared_ptr<std::vector<size_t>> stoiSort(const std::vector<T> &vector);
    extern template std::shared_ptr<std::vector<size_t>> stoiSort<std::string_view>(const std::vector<std::string_view> &vector);

    template <typename T>
    std::shared_ptr<std::vector<size_t>> getRefOrder(const std::vector<T>& reference, const std::vector<T>& to_sort);
    extern template std::shared_ptr<std::vector<size_t>> getRefOrder<std::string_view>(const std::vector<std::string_view>& reference, const std::vector<std::string_view>& to_sort);

    std::shared_ptr<std::vector<std::string_view>> nuWordSplitter( std::string_view line );

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
        std::vector<std::string_view> getMom();
        std::string_view getE();
        std::string_view getMass();
        std::string_view getVTim();
        std::string_view getSpin();
        std::string_view getPDG();
        std::string_view getStatus();
        std::vector<std::string_view> getMothers();
        std::vector<std::string_view> getColor();
        void setComment( std::string_view nuCom );
        void setMom( std::vector<std::string_view> nuMom );
        void setEnergy( std::string_view nuE );
        void setMass( std::string_view nuM );
        void setVTim( std::string_view nuVTim );
        void setSpin( std::string_view nuSpin );
        void setPDG( std::string_view nuPDG );
        void setStatus( std::string_view nuSt );
        void setMothers( std::vector<std::string_view> nuMum );
        void setColors( std::vector<std::string_view> nuCol );
        bool isModded();
        bool isWritten();
        std::shared_ptr<std::string> getContent();
        lhePrt();
        lhePrt( std::pair<int,int>& prtInfo );
        lhePrt( std::pair<std::string,std::string>& prtInfo );
        lhePrt( const std::string_view originFile, const size_t& beginLine = 0, const size_t& endLine = npos );
    protected:
        std::shared_ptr<std::string> content;
        std::string_view sourceFile;
        std::string_view comment;
        std::string_view mom[3];
        std::string_view energy;
        std::string_view mass;
        std::string_view vtim;
        std::string_view spin;
        std::string_view pdg;
        std::string_view status;
        std::string_view mothers[2];
        std::string_view icol[2];
        bool modded = false;
        bool written = false;
        void writer();
    };

    struct evHead {
    public:
        std::string_view getComment();
        std::string_view getWeight();
        std::string_view getScale();
        std::string_view getAQED();
        std::string_view getAQCD();
        std::string_view getNprt();
        std::string_view getProcID();
        bool isModded();
        bool isWritten();
        void setComment( std::string_view nuCom );
        void setWeight( std::string_view nuWgt );
        void setScale( std::string_view nuScale );
        void setAQED( std::string_view nuAQED );
        void setAQCD( std::string_view nuAQCD );
        void setNprt( std::string_view nuNprt );
        void setNprt( int nuNprt );
        void setProcID( std::string_view nuProcID );
        std::shared_ptr<std::string> getContent();
        evHead();
        evHead( const std::string_view originFile, size_t beginLine = 0, size_t endLine = npos );
    protected:
        std::shared_ptr<std::string> content;
        std::string_view sourceFile;
        std::string_view comment;
        std::string_view weight;
        std::string_view scale;
        std::string_view aqed;
        std::string_view aqcd;
        std::string_view nprt;
        int nprtint;
        std::string nprtstr;
        std::string_view procid;
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
        std::map<std::string_view, std::vector<std::string_view>> procMap;
        std::map<std::string_view, std::vector<size_t>> procOrder;
        sortFcn eventSort = []( std::vector<std::string_view> vec ){ return stoiSort( vec ); };
        statSort specSort = []( std::string_view stat, std::vector<std::string_view> vec ){ UNUSED(stat); return stoiSort( vec ); };
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
        std::map<std::string_view, std::vector<std::string_view>> &getProc();
        std::map<std::string_view, std::vector<size_t>> &getProcOrder();
        std::map<std::string_view, std::vector<std::string_view>> getProc() const;
        std::map<std::string_view, std::vector<size_t>> getProcOrder() const;
        std::map<std::string_view, std::vector<std::string_view>> &getProc(sortFcn sorter);
        std::map<std::string_view, std::vector<size_t>> &getProcOrder(sortFcn sorter);
        std::map<std::string_view, std::vector<std::string_view>> &getProc(statSort sorter);
        std::map<std::string_view, std::vector<size_t>> &getProcOrder(statSort sorter);
    };

    using eventComparison = std::function<bool(event&, event&, std::vector<std::string>&)>;

    using eventSetComp = std::function<bool(event&, std::vector<std::string>&)>;

    struct eventSet{
        eventSet();
        eventSet( const eventSet& nuEvents );
        eventSet( std::vector<event>& nuEvents );
        eventSet( std::vector<std::shared_ptr<event>>& nuEvents );
        void setRelStats( std::vector<std::string>& nuStats );
        void addEvent( event& nuEvent );
        void addEvent( std::shared_ptr<event> nuEvent );
        void addEvent( std::vector<event>& nuEvents );
        void addEvent( std::vector<std::shared_ptr<event>> nuEvents );
        void setComp( eventSetComp nuComp );
        bool belongs( event& nuEvent );
        bool belongs( std::shared_ptr<event> nuEvent );
    protected:
        std::vector<event> events;
        std::vector<std::string> relStats = {"-1", "1"};
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
        std::vector<std::shared_ptr<weightGroup>> groups;
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
        sortFcn particleSort = []( std::vector<std::string_view> prts ){ return stoiSort(prts); };
        statSort statParticleSort = []( std::string_view dummy, std::vector<std::string_view> prts ){ UNUSED(dummy); return stoiSort(prts); };
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
        std::vector<std::string_view> wgts;
        std::vector<std::string_view> scales;
        std::vector<std::string_view> aQEDs;
        std::vector<std::string_view> aQCDs;
        std::vector<std::string_view> nprts;
        std::vector<size_t> relNPrts;
        std::vector<std::string_view> procIDs;
        evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile = {} );
        evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const std::vector<std::string_view>& statVec );
        evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const std::vector<std::string_view>& statVec, 
        sortFcn sorter );
        evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const std::vector<std::string_view>& statVec, 
        statSort sorter );
    };

    struct  prtInfo {
    public:
        std::vector<std::string_view> moms;
        std::vector<std::string_view> masses;
        std::vector<std::string_view> vtims;
        std::vector<std::string_view> spins;
        std::vector<std::string_view> statuses;
        std::vector<std::string_view> mothers;
        std::vector<std::string_view> icols;
        std::vector<std::string_view> pdgs;
        prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile = {}, const int nPrt = 8 );
        prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, const std::vector<std::string_view>& statVec );
        prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, const std::vector<std::string_view>& statVec, 
        sortFcn sorter );
        prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, const std::vector<std::string_view>& statVec, 
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
        transMonoLHE( const std::vector<std::shared_ptr<REX::event>> lheFile, const int nPrt, const std::vector<std::string_view>& statVec );
        transMonoLHE( const std::vector<std::shared_ptr<REX::event>> lheFile, const int nPrt, 
        sortFcn sorter,
        std::vector<std::string_view> statVec = { "-1", "1" } );
        transMonoLHE( const std::vector<std::shared_ptr<REX::event>> lheFile, const int nPrt, 
        statSort sorter,
        std::vector<std::string_view> statVec = { "-1", "1" } );
    };

    struct transLHE {
    public:
        std::string_view xmlFile;
        std::vector<std::shared_ptr<transMonoLHE>> subProcs;
        std::vector<std::shared_ptr<event>> procSets;
        std::vector<std::shared_ptr<std::vector<bool>>> relProcs;
        transLHE();
        transLHE( lheNode& lheFile );
        transLHE( lheNode& lheFile, 
        sortFcn sorter, 
        const std::vector<std::string_view>& statVec = { "-1", "1" } );
        transLHE( lheNode& lheFile, 
        statSort sorter, 
        const std::vector<std::string_view>& statVec = { "-1", "1" } );
        transLHE( lheNode& lheFile, const std::vector<std::string_view>& statVec );
        transLHE( transSkel& skeleton );
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
        bool operator()(event& firstEv, event& secEv, std::vector<std::string_view> statVec);
    };


std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> lheValDoubles( lheNode& lheFile, lheRetDs vals = lheRetDs() );

std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> lheValDoubles(transLHE& lheAOS, lheRetDs vals = lheRetDs() );

//    struct lhePrt;
//    struct xmlNode;
//    struct event : public xmlNode;
//    event& makeEv( std::vector<std::pair<int,int>>& particles );
//    std::vector<std::shared_ptr<lhePrt>> getParticles( event& ev );
//    struct eventComp;
}

#endif
