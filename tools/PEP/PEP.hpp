/***
 *      _____  ______ _____  
 *     |  __ \|  ____|  __ \ 
 *     | |__) | |__  | |__) |
 *     |  ___/|  __| |  ___/ 
 *     | |    | |____| |     
 *     |_|    |______|_|     
 *                                                     
 ***/
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
#include <stdexcept>

// ZW: all fcns within the PEP standard sit in the
// namespace PEP
// note that as a convention, std::string_view objects
// will be referred to as strings unless the difference
// is relevant
namespace PEP
{

    // ZW: minimal fcn for counting the amount
    // of times a given search term appears in
    // a string
    int nuStrCount( std::string_view searchString, std::string_view searchTerm )
    {
        int count = 0;
        size_t pos = 0;
        while((pos = searchString.find(searchTerm, pos)) != std::string::npos ){
            ++count;
            ++pos;
        }
        return count;
    }

    // ZW: fcn for finding the location of each
    // entry of seachTerm in the given string
    // textFile. Pre-allocates vector memory using
    // nuStrCount
    std::shared_ptr<std::vector<size_t>> nuFindEach( std::string_view textFile, std::string_view searchTerm )
    {
        auto eachPos = std::make_shared<std::vector<size_t>>();
        eachPos->reserve( nuStrCount(textFile, searchTerm) );
        eachPos->push_back( textFile.find( searchTerm ) );
        size_t currPos = textFile.find( searchTerm, eachPos->at(0) + 1 );
        while( currPos != std::string_view::npos )
        {
            eachPos->push_back( currPos );
            currPos = textFile.find( searchTerm, currPos + 1 );
        }
        return eachPos;
    }

    // ZW: fcn for splitting a string into a vector of strings,
    // each element differentiated by linebreaks
    // in the original string
    // removes sequential linebreaks, ie "\n\n\n" would
    // only result in a single element separation
    std::shared_ptr<std::vector<std::string_view>> nuLineSplitter( std::string_view currEvt )
    {
        auto lineBreaks = nuFindEach( currEvt, "\n" );
        std::vector<size_t> trueBreaks;
        trueBreaks.reserve( lineBreaks->size() );
        for( int k = 0 ; k < lineBreaks->size() - 1 ; ++k )
        {
            if( int( (*lineBreaks)[k+1] - (*lineBreaks)[k]) == 1){continue;}
            trueBreaks.push_back( (*lineBreaks)[k] );
        }
        auto splitLines = std::make_shared<std::vector<std::string_view>>();
        splitLines->reserve( trueBreaks.size() );
        size_t startPos = 0;
        for( auto k : trueBreaks )
        {
            splitLines->push_back( currEvt.substr( startPos + 1, k - startPos - 1) );
            startPos = k;
        }
        if( auto strung = currEvt.substr( startPos ).size() > 1 ){ splitLines->push_back( currEvt.substr( startPos ) ); }
        return splitLines;
    }

    // ZW: fcn for finding each linebreak in a string,
    // returning a vector of the positions of "\n" characters
    // Ignores sequential linebreaks, ie would only return { }
    // for the string "\n\n\n\n"
    std::shared_ptr<std::vector<size_t>> lineFinder( std::string_view currEvt, size_t startPos = 0, size_t endPos = std::string_view::npos )
    {
        auto lineBreaks = nuFindEach( currEvt.substr( startPos, endPos - startPos), "\n" );
        auto truBreaks = std::make_shared<std::vector<size_t>>();
        truBreaks->reserve( lineBreaks->size() );
        for( int k = 0 ; k < lineBreaks->size() ; ++k )
        {
            if( int( (*lineBreaks)[k+1] - (*lineBreaks)[k]) == 1){continue;}
            truBreaks->push_back( (*lineBreaks)[k] );
        }
        return truBreaks;
    }
    
    // ZW: fcn for splitting a string into a vector
    // of strings, each element separated by
    // blankspace (" ") in the original string
    // ignores sequential blankspaces, as well
    // as linebreaks, ie "hello     \n\n\n     world"
    // would return {"hello", "world"}
    // does not ignore linebreaks that are not
    // separated from words by anything other than
    // blankspace, however
    // ie "hello     \n\n\nworld   \n\n"
    // would return {"hello", "\n\nworld"}
    std::shared_ptr<std::vector<std::string_view>> nuWordSplitter( std::string_view currEvt )
    {
        std::vector<size_t> noSpace;
        size_t nuStart = currEvt.find_first_not_of( " " );
        size_t nuEnd = currEvt.find(" ", nuStart+1 );
        auto splitWords = std::make_shared<std::vector<std::string_view>>();
        splitWords->reserve(13);
        while( nuStart != std::string_view::npos )
        {
            std::string_view word = currEvt.substr( nuStart, nuEnd - nuStart );
            if( word == "" || word == "\n" || word == " " ){
                nuStart = currEvt.find_first_not_of(" ", nuEnd);
                nuEnd = currEvt.find( " ", nuStart + 1);
                continue; }
            splitWords->push_back( currEvt.substr( nuStart, nuEnd - nuStart ) );
            nuStart = currEvt.find_first_not_of(" ", nuEnd);
            nuEnd = currEvt.find( " ", nuStart + 1);
        }
        return splitWords;
    }

    // ZW: fcn for splitting a string into a vector of strings,
    // elements separated by any form of blankspace in the
    // original string
    // ignores sequential blankspaces of all forms
    std::shared_ptr<std::vector<std::string_view>> nuBlankSplitter( std::string_view currEvt )
    {
        auto lines = nuLineSplitter( currEvt );
        auto splitString = std::make_shared<std::vector<std::string_view>>();
        splitString->reserve( lines->size() * lines->at(0).size() );
        for( auto line : *lines )
        {
            auto words = nuWordSplitter(line);
            for( auto word : *words )
            {
                if( word == "" || word == "\n" || word == " " ){continue;}
                splitString->push_back( word );
            }
        }
        return splitString;
    }

    // ZW: templated fcn for comparing two
    // string-like objects, ignoring cases
    template<typename Str1, typename Str2>
    bool clStringComp( const Str1& org, const Str2& comp ){
        return std::equal( org.begin(), org.end(), comp.begin(), comp.end(), 
        []( const char& x, char y ){ return (std::toupper(x) == std::toupper(y)); } );
    }
    template<typename Str1Pt, typename Str2>
    bool clStringComp( const Str1Pt& orgStrt, const Str1Pt& orgEnd, const Str2& comp ){
        return std::equal( orgStrt, orgEnd, comp.begin(), comp.end(), 
        []( const char& x, char y ){ return (std::toupper(x) == std::toupper(y)); } );
    }

    // ZW: templated function for finding a 
    // caseless substring searchTerm in srcFile
    // with size_t(-1) returned on failure to find
    // the considered substring
    template<typename Str1, typename Str2>
    size_t clStringFind( const Str1& srcFile, const Str2& searchTerm, size_t strtPt = 0 ){
        size_t strLen = searchTerm.size();
        if( srcFile.size() == 0  || srcFile.size() < strLen ){ return size_t(-1); }
        for( size_t k = strtPt ; k < srcFile.size()  - strLen; ++k )
        {
            if( clStringComp( srcFile.substr(k, strLen), searchTerm ) ){ return k; }
        }
        return size_t(-1);
    }

    // ZW: templated fcn for finding a caseless
    // substring searchTerm of srcFile fulfilling
    // a particular predicate cond( size_t, string )
    template<typename Str1, typename Str2>
    size_t clStringFindIf( const Str1& srcFile, const Str2& searchTerm, std::function<bool(size_t&, const Str1&)>& cond, size_t strtPt = 0 )
    {
        auto currPt = clStringFind( srcFile, searchTerm, strtPt ); 
        bool condStat = cond( currPt, srcFile );
        while( !( condStat ) && currPt != size_t(-1))
        {
            currPt = clStringFind( srcFile, searchTerm, currPt + 1 );
            condStat = cond( currPt, srcFile );
        } 
        return currPt;
    }

    // ZW: templated fcn for counting the number
    // of occurances of the caseless substring 
    // searchTerm in string-like object srcFile
    template<typename Str1, typename Str2>
    int clStrCount( Str1 srcFile, Str2 searchTerm )
    {
        int count = 0;
        size_t pos = 0;
        while((pos = clStringFind( srcFile, searchTerm, pos ) ) != size_t(-1) ){
            ++count;
            ++pos;
        }
        return count;
    }

    // ZW: templated fcn for finding each instance
    // of substring searchTerm of string-like
    // object srcFile
    template<typename Str1, typename Str2>
    std::shared_ptr<std::vector<size_t>> clFindEach( Str1 srcFile, Str2 searchTerm )
    {
        auto eachPos = std::make_shared<std::vector<size_t>>();
        auto nos = clStrCount(srcFile, searchTerm);
        if( nos == 0 ){ return eachPos; }
        eachPos->reserve( nos );
        eachPos->push_back( clStringFind( srcFile, searchTerm ) );
        size_t currPos = clStringFind( srcFile, searchTerm, eachPos->at(0) + 1);
        while( currPos != size_t(-1) )
        {
            eachPos->push_back( currPos );
            currPos = clStringFind( srcFile, searchTerm, currPos + 1 );
        }
        return eachPos;
    }

    // ZW: fcn for finding left angle bracket
    // indicating the start of a new node
    // in an XML file
    std::shared_ptr<size_t> nodeStartFind( std::string_view parseFile, size_t strtPos )
    {
        auto retPtr = std::make_shared<size_t>(parseFile.find("<", strtPos));
        while( parseFile[*retPtr + 1] == '!' || parseFile[*retPtr +1] == '/' || parseFile[*retPtr +1] == '?' || parseFile[*retPtr +1] == '-' ){
            *retPtr = parseFile.find("<", *retPtr +1);
        }
        return retPtr;
    }

    // ZW: fcn for finding left angle bracket
    // indicating an end of a node
    // in an XML file
    std::shared_ptr<size_t> nodeEndFind( std::string_view parseFile, size_t strtPos )
    {
        auto retPtr = std::make_shared<size_t>(parseFile.find("<", strtPos));
        while( parseFile[*retPtr + 1] != '/' ){
            *retPtr = parseFile.find("<", *retPtr +1);
        }
        return retPtr;
    }

    // ZW: struct for handling tags in XML
    // node opening tags
    struct xmlTag {
    public:
        void setVal( std::string_view valSet ){ modded = true; val = valSet; }
        void setId( std::string_view idSet ){ modded = true; id = idSet; }
        std::string_view getVal(){ return val; }
        std::string_view getId(){ return id; }
        bool isModded(){ return modded; }
        xmlTag(){ modded = false; return; }
        xmlTag( xmlTag& oldTag ){
            modded = false; val = oldTag.getVal(); id = oldTag.getId();
        }
        xmlTag( std::string_view initId, std::string_view initVal){
            modded = false; val = initVal; id = initId;
        }
    protected:
        bool modded;
        std::string_view val;
        std::string_view id;
    };

    // ZW: function for parsing XML opening
    // tags and returning the next header tag
    std::shared_ptr<xmlTag> xmlTagParser( std::string_view tagLine, size_t& equPt )
    {
        //auto equPt = tagLine.find("=");
        auto tagBreaker = tagLine.find_first_not_of(" ", equPt+1);
        auto tagEnder = tagLine.find( tagLine[tagBreaker], tagBreaker+1);
        auto attrEnd = tagLine.find_last_not_of(" ", equPt - 1) ;
        auto attrStart = tagLine.find_last_of(" ", attrEnd) + 1;
        auto tagPtr = std::make_shared<xmlTag>(tagLine.substr(attrStart, attrEnd - attrStart + 1), tagLine.substr(tagBreaker + 1, tagEnder - tagBreaker - 1));
        equPt = tagLine.find("=", equPt + 1);
        return tagPtr;
    }

    // ZW: generic struct for handling XML
    // nodes in generic XML files
    struct xmlNode {
    public:
        xmlNode(){ modded = false; return; }
        xmlNode( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} ){
            modded = false; xmlFile = originFile; start = begin; children = childs;
            if( xmlFile.substr(start, 1) != "<" ){ start = *nodeStartFind( xmlFile, size_t(start) ); }
            size_t trueStart = xmlFile.find_first_not_of(" ", start+1);
            name = xmlFile.substr( trueStart, xmlFile.find_first_of(">/ ", trueStart) - trueStart );
            if( xmlFile.find( ">", trueStart ) < xmlFile.find( "/", trueStart ) ){
                content = xmlFile.substr( xmlFile.find( ">", trueStart ) + 1, xmlFile.find( "</", trueStart ) - xmlFile.find( ">", trueStart ) - 1 );
            }
        }
        std::vector<std::shared_ptr<xmlNode>> getChildren(){ return children; }
        std::vector<std::shared_ptr<xmlTag>> getTags(){ return tags; }
        std::string_view getFile(){ return xmlFile; }
        std::string_view getName(){ return name; }
        std::string_view getContent(){ return content; }
        size_t getStart(){ return start; }
        size_t getEnd(){ return end; }
        virtual bool isModded(){ return modded; }
        virtual bool isModded( bool deep ){
            bool modStat = isModded();
            if( !deep ){ return modStat; }
            for( auto child : children ){ modStat = (modStat || child->isModded( deep )); }
            return modStat;
        }
        bool isWritten(){ return written; }
        bool isParsed(){ return parsed; }
        void setModded( bool mod ){ modded = mod; }
        bool deepModded(){ return deepMod; }
        bool deepParse(){ return deepParsed; }
        void parser( bool recursive ){
            parsed = parse( recursive );
        }
        void addChild( std::shared_ptr<xmlNode> child ){ modded = true; children.push_back(child); }
        void addTag( std::shared_ptr<xmlTag> tag ){ modded = true; tags.push_back(tag); }
        void setFile( std::string_view file ){ modded = true; xmlFile = file; }
        void setName( std::string_view newName ){ modded = true; name = newName; }
        void setCont( std::string_view cont ){ modded = true; content = cont; }
    protected:
        virtual bool parse(){ 
            auto topStat = parseTop();
            auto contStat = parseContent();
            return ( topStat && contStat );
        }
        virtual bool parse( bool recurs )
        {
            bool parseSt = parse();
            if( !recurs ){ return parseSt; }
            bool childSt = parseChildren( recurs );
            deepMod = true;
            return (parseSt && childSt );
        }
        bool parseTop(){
            if( xmlFile == "" ){ return false; }
            size_t eqSgn = xmlFile.find( "=", start ); size_t nodeInitEnd = xmlFile.find( ">", start );
            while( eqSgn < nodeInitEnd ){ tags.push_back( xmlTagParser( xmlFile, eqSgn ) ); }
            return true;
        }
        virtual bool parseContent(){
            if( xmlFile == "" ){ return false; }
            auto firstR = xmlFile.find_first_of( ">/", start );
            auto nodeStrEnd = xmlFile.find(">", firstR);
            if( firstR < nodeStrEnd ){ content = ""; end = nodeStrEnd + 2; parsed = true; return true; }
            auto endNode = *nodeEndFind( xmlFile, start );
            auto startNode = *nodeStartFind( xmlFile, start + 1 );
            if( startNode > endNode ){end = xmlFile.find( ">", endNode ) + 1; content = xmlFile.substr( xmlFile.find( ">", start ) + 1, endNode - xmlFile.find( ">", start ) - 1  ); return true; }
            auto endPt = xmlFile.find( std::string("</") + std::string(name), start );
            content = xmlFile.substr( xmlFile.find(">", start) + 1, startNode - xmlFile.find(">") - 1 ); 
            end = xmlFile.find( ">", endPt ) + 2; 
            while( startNode < endNode ){
                auto nextNode = std::make_shared<xmlNode>( xmlFile, startNode );
                children.push_back( nextNode );
                int starts = 0;
                while( startNode < endNode )
                {
                    startNode = *nodeStartFind( xmlFile, startNode + 1 );
                    ++starts;
                }
                for( int k = 0 ; k < starts ; ++k ){ endNode = *nodeEndFind( xmlFile, endNode + 1 ); }
                if( endNode > end ){ break; }
            }
            return true;
        }
        bool parseChildren( bool recursive ){
            bool status = true;
            if( recursive ){
                for( auto child : children )
                {
                    status = (status && child->parse( true ));
                    deepParsed = true;
                }
            } else {
                for( auto child : children )
                {
                    status = (status && child->parse());
                    deepParsed = true;
                }
            }
            return status;
        }
        std::shared_ptr<std::string> writtenSelf; 
        bool deepMod = false;
        std::vector<std::shared_ptr<xmlNode>> children;
        std::vector<std::shared_ptr<xmlTag>> tags;
        std::string_view xmlFile;
        std::string_view name;
        std::string_view content;
        size_t start;
        size_t end = std::string_view::npos;
        bool modded = false;
        bool written = false;
        bool parsed = false;
        bool deepParsed = false;
        std::string nodeHeader;
        std::string nodeContent;
        std::string nodeEnd;
        virtual void headWriter() {
            nodeHeader =  "<" + std::string(name) ;
            for( auto tag : tags ){
                nodeHeader += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            nodeHeader += ">";
        }
        virtual void endWriter() {
            nodeEnd = "</" + std::string(name) + ">\n";
        }
        virtual void contWriter() {
            if( children.size() > 0 ){
            nodeContent = std::string(content.substr(0, children[0]->start - 1 ));
            } else {
            nodeContent = std::string(content);
            }
        }
        virtual void childWriter() {
            for(auto child : children){
                nodeContent += (*child->nodeWriter());
            }
        }
        virtual void endFinder(){
            auto headEnd = xmlFile.find(">", start);
            auto slashPos = xmlFile.find("/", start);
            if( headEnd > slashPos ){ end = headEnd; }
            else{ end = xmlFile.find( ">", xmlFile.find( "</" + std::string(name), start )); }
            if( end == std::string_view::npos ){ end = xmlFile.size(); return; }
            end += 2;
        }
        virtual void fullWriter(){
            if( isModded() ){
            headWriter();
            contWriter();
            childWriter();
            endWriter();
            writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
            written = true;
            modded = false;
            } else if( !isWritten() ){
            endFinder();
            if( start > xmlFile.size() ){ start = 0; }
            writtenSelf = std::make_shared<std::string>( xmlFile.substr( start, end - start ) );
            written = true;
            }
        }
    public:
        virtual void childCounter( int& noChilds )
        {
            for( auto child : children )
            {
                child->childCounter( noChilds );
                if( child->end == 0 ){ --noChilds; }
            }
            noChilds += children.size();
        }    
        virtual std::shared_ptr<std::string> nodeWriter() {
            if( isModded( true ) || !isWritten() ){ fullWriter(); }
            return writtenSelf;
        }
    };

    // ZW: function for large scale parsing of XML files
    // sequentially goes through the document and
    // recursively calls itself while the next node
    // beginning is closer than the next node ending
    std::shared_ptr<xmlNode> xmlPtrParser( std::string_view parseFile, size_t& initPos, size_t& endPos )
    {
        auto currNode = std::make_shared<xmlNode>(parseFile, initPos);
        size_t equalSign = parseFile.find("=", initPos);
        size_t nodeInitEnd = parseFile.find(">", initPos);
        initPos = *nodeStartFind( parseFile, initPos + 1 );
        while( equalSign < nodeInitEnd ){
            currNode->addTag( xmlTagParser(parseFile, equalSign) );
        }
        while( initPos < endPos )
        {
            currNode->addChild(xmlPtrParser( parseFile, initPos, endPos ));
        }
        
        initPos = *nodeStartFind( parseFile, endPos );
        endPos = *nodeEndFind( parseFile, endPos + 1 );
        return currNode;
    }

    // ZW: struct for handling rwgt parameter sets
    // in the LHE header initrwgt node
    struct headWeight : xmlNode {
    public:
        int getId(){ return id; }
        std::string_view getTag(){ return idTag; }
        headWeight(){ name = "weight"; return; }
        headWeight( std::string_view paramSet, const size_t& begin = 0 ) : xmlNode(){ name = "weight"; xmlFile = paramSet; content = paramSet; return; }
        headWeight( std::string_view paramSet, std::string_view idText, int idNo, const size_t& begin = 0 ) : xmlNode(){
            name = "weight"; xmlFile = paramSet; content = paramSet; idTag = idText; id = idNo;
        }
        headWeight( xmlNode& node ) : xmlNode( node ){
            parser( false );
            name = "weight";
            for (auto tag : tags ){
                if( tag->getId() == "id" ){
                    idTag = tag->getVal().substr(0, tag->getVal().find_last_of("_") - 1 );
                    id = std::stoi( std::string( tag->getVal().substr( idTag.size() + 1 ) ) );
                }
            }
        }
        headWeight( xmlNode* node ) : xmlNode( *node ){
            parser( false );
            name = "weight";
            for (auto tag : tags ){
                if( tag->getId() == "id" ){
                    idTag = tag->getVal().substr(0, tag->getVal().find_last_of("_") - 1 );
                    id = std::stoi( std::string( tag->getVal().substr( idTag.size() + 1 ) ) );
                }
            }
        }
        headWeight( std::shared_ptr<xmlNode> node ) : xmlNode( *node ){
            parser( false );
            name = "weight";
            for (auto tag : tags ){
                if( tag->getId() == "id" ){
                    idTag = tag->getVal().substr(0, tag->getVal().find_last_of("_") - 1 );
                    id = std::stoi( std::string( tag->getVal().substr( idTag.size() + 1 ) ) );
                }
            }
        }
        headWeight( std::string_view paramSet, std::string& idText, unsigned int idNo, const size_t& begin = 0 ) : xmlNode(){
            name = "weight"; xmlFile = paramSet; content = paramSet; idTag = idText; id = idNo;
        }
        headWeight( std::string_view paramSet, std::string& idText){
            name = "weight"; xmlFile = paramSet; content = paramSet; idTag = idText;
        }
    protected:
        std::string idTag;
        unsigned int id = -1;
        void headWriter() override{
            if( tags.size() == 0 ){
                if( idTag == "" ){ nodeHeader = "<weight>"; return; }
                if( id == -1 ){ nodeHeader = "<weight id=\"" + std::string(idTag) + "\">"; return; }
                nodeHeader = "<weight id=\"" + std::string(idTag) + std::to_string(id) + "\">";
            }
            nodeHeader = "<weight";
            for( auto tag : tags ){
                nodeHeader += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            nodeHeader += ">";
        }
        void headWriter( bool incId ){
            if( !incId ){ headWriter(); return; }
            if( idTag == "" ){ headWriter(); return; }
            if( id == -1 ){ nodeHeader = "<weight id=\"" + std::string( idTag ) + "\""; }
            else{ nodeHeader = "<weight id=\"" + std::string( idTag ) + "_" + std::to_string(id) + "\""; }
            for( auto tag : tags ){
                if( tag->getId() == "id" ){ continue; }
                nodeHeader += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            nodeHeader += ">";
        }
        void endWriter() override{
            nodeEnd = "</weight>\n";
        }
        void contWriter() override{ 
            nodeContent = std::string( content );
        }
        void childWriter() override{
            for( auto child : children){
                if( child->getName() == "weight" ){ continue; }
                nodeContent += *(child->nodeWriter());
            }
        }
        void childWriter( bool hasChildren ){
            if( hasChildren ){ childWriter(); }
        }
        void fullWriter() override{
            if( isModded() || !isWritten() ){
                headWriter();
                contWriter();
                childWriter();
                endWriter();
                writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
                writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
                written = true;
                modded = false;
            }
        }
        void fullWriter( bool incId, bool hasChildren=true ){
            if( isModded() || !isWritten() ){
            headWriter( incId );
            contWriter();
            childWriter( );
            endWriter();
            writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
            modded = false;
            written = true;
            }
        }
    };

    // ZW: struct for handling rwgt groups
    // in the LHE header initrwgt node
    struct weightGroup : xmlNode {
    public:
        bool getIncId(){ return includeId; }
        void setIncId( bool nuIncId ){ includeId = nuIncId; }
        std::vector<std::shared_ptr<headWeight>> getWgts(){ return paramSets; }
        void addWgt( headWeight nuWgt ){ modded = true; paramSets.push_back( std::make_shared<headWeight>( nuWgt ) ); }
        void addWgt( std::shared_ptr<headWeight> nuWgt ){ modded = true; paramSets.push_back( nuWgt); }
        weightGroup() : xmlNode(){ name = "weightgroup"; return; }
        weightGroup( std::vector<std::shared_ptr<headWeight>> nuWgts ) : xmlNode(){ name = "weightgroup"; paramSets = nuWgts; }
        weightGroup( std::vector<std::string> nuWgts ) : xmlNode(){
            name = "weightgroup";
            for( auto wgt : nuWgts ){
                paramSets.push_back( std::make_shared<headWeight>( wgt ) );
            }
        }
        weightGroup( xmlNode& wgtNode ) : xmlNode( wgtNode ){
            parser( true );
            name = "weightgroup";
            paramSets.reserve( children.size() );
            for(  auto child : children ){
                if( child->getName() == "weight" ){ paramSets.push_back( std::make_shared<headWeight>( *child ) ); }
            }
        }
        weightGroup( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} )
        : xmlNode( originFile, begin, childs ){
            name = "weightgroup";
            if( parseTop() ){
                int checker = 0;
                for( auto tag : tags ){
                    if( tag->getId() == "name" ){ ++checker; rwgtName = tag->getVal(); }
                    if( tag->getId() == "weight_name_strategy" ){ ++checker; wgtNamStrat = tag->getVal();
                        if(wgtNamStrat == "includeIdInWeightName"){ includeId = true; } }
                    if( checker == 2 ){ break; }
                }
            }
        }
    protected:
        std::string_view rwgtName;
        std::string_view wgtNamStrat;
        bool includeId = false;
        std::vector<std::shared_ptr<headWeight>> paramSets;
        bool nu;
        std::string_view idTag;
        int id;
        void headWriter() override{
            nodeHeader = "<weightgroup";
            if( rwgtName !="" ){ nodeHeader += " name=\"" + std::string( rwgtName ) +"\""; }else if( isModded() ){ nodeHeader += " name=\"pep_reweighting\""; }
            if( wgtNamStrat!="" ){ nodeHeader += " weight_name_strategy=\"" + std::string( wgtNamStrat ) +"\""; } 
            nodeHeader += ">";
        }
        void contWriter() override{
            nodeContent = "\n";
            for( auto wgt : paramSets ){
                nodeContent += (*wgt->nodeWriter());
            }
        }
        void childWriter( bool hasChildren = false ){
            if( hasChildren ){ childWriter(); }
            return;
        }
        void endWriter() override{ nodeEnd = "</weightgroup>\n"; }
    };

    struct initRwgt : xmlNode {
    public:
        std::vector<std::shared_ptr<weightGroup>> getGroups(){ return groups; }
        size_t noGrps(){ return groups.size(); }
        void addGroup( weightGroup nuGroup ){ 
            modded = true;
            auto nuGrpPtr = std::make_shared<weightGroup>( nuGroup );
            if( grpInit( nuGrpPtr ) ){ groups.push_back( std::make_shared<weightGroup>( nuGroup ) );  }
        }
        void addGroup( std::shared_ptr<weightGroup> nuGroup ){ 
            modded = true;
            if( grpInit( nuGroup ) ){ groups.push_back( nuGroup ); }
        }
        void addWgt( unsigned int index, std::shared_ptr<headWeight> nuWgt ){
            if( index < groups.size() ){ modded = true; groups[index]->addWgt( nuWgt ); }
            else throw std::out_of_range( "Appending weight to uninitialised weightgroup." );
        }
        void addWgt( unsigned int index, headWeight nuWgt ){
            if( index < groups.size() ){ modded = true; groups[index]->addWgt( nuWgt ); }
            else throw std::out_of_range( "Appending weight to uninitialised weightgroup." );
        }
        initRwgt() : xmlNode(){ name = "initrwgt"; return; }
        initRwgt( std::vector<std::shared_ptr<xmlNode>> nuGroups ) : xmlNode(){
            name = "initrwgt";
            for( auto group : nuGroups ){
                groups.push_back( std::make_shared<weightGroup>( *group ) );
            }
        }
        initRwgt( xmlNode& wgtNode ) : xmlNode( wgtNode ){
            parser( true );
            name = "initrwgt";
            groups.reserve( children.size() );
            for(  auto child : children ){
                groups.push_back( std::make_shared<weightGroup>( *child ) );
            }
        }
        initRwgt( std::shared_ptr<xmlNode> wgtNode ) : xmlNode( *wgtNode ){
            parser( true );
            name = "initrwgt";
            groups.reserve( children.size() );
            for(  auto child : children ){
                groups.push_back( std::make_shared<weightGroup>( *child ) );
            }
        }
    protected:
        bool grpIsInit = false;
        bool grpInit( std::shared_ptr<weightGroup>& wgt ){
            if( grpIsInit ){ return true; }
            else{
                groups = std::vector<std::shared_ptr<weightGroup>>( 1, wgt );
                grpIsInit = true;
                return false;
            }
        }
        std::vector<std::shared_ptr<weightGroup>> groups;
        void contWriter() override{
            nodeContent = "\n";
            for( auto group : groups ){
                nodeContent += (*group->nodeWriter());
            }
        }
        void childWriter() override{
            for( auto child : children ){
                if( child->getName() == "weightgroup" ){ continue; }
                nodeContent += (*child->nodeWriter());
            }
        }
        void childWriter( bool hasChildren ){
            if( hasChildren ){ childWriter(); }
            return;
        }
    };

    // ZW: struct for handling event
    // in event blocks of LHE files
    struct bodyWgt : xmlNode {
    public:
        void setComment( std::string_view nuComment ){ modded = true; comment = nuComment; }
        void setVal( std::string nuVal ){ modded = true; valS = nuVal; valD = std::stod(valS);}
        void setVal( std::string_view nuVal ){ modded = true; valS = std::string(nuVal); valD = std::stod(valS);}
        void setVal( double nuVal ){ modded = true; valD = nuVal; valS = std::to_string(valD);}
        void setModded( bool nuModded ){ modded = nuModded; }
        std::string_view getComment(){ return comment; }
        std::string_view getValS(){ return valS; }
        double getValD(){ return valD; }
        bodyWgt() : xmlNode(){ return; }
        bodyWgt( std::string_view value ) : xmlNode() { setVal( value ); modded = false; }
        bodyWgt( double value ) : xmlNode() { setVal( value ); modded = false; }
        bodyWgt( std::string_view value, xmlTag rwgtId ) : xmlNode() { setVal( value ); addTag( std::make_shared<xmlTag>(rwgtId) ); modded = false; }
        bodyWgt( double value, xmlTag rwgtId ) : xmlNode() { setVal( value ); addTag( std::make_shared<xmlTag>(rwgtId) ); modded = false; }
        bodyWgt( std::string_view value, std::shared_ptr<xmlTag> rwgtId ) : xmlNode() { setVal( value ); addTag( rwgtId ); modded = false; }
        bodyWgt( double value, std::shared_ptr<xmlTag> rwgtId ) : xmlNode() { setVal( value ); addTag( rwgtId ); modded = false; }
        bodyWgt( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} )
        : xmlNode( originFile, begin, childs ){
            auto strtPt = originFile.find_first_not_of(" >+", originFile.find(">", begin)+1);
            valS = originFile.substr( strtPt, originFile.find(" ", strtPt) - strtPt );
            valD = std::stod( valS );
        }
        bodyWgt( double value, std::string& idTag ){
            setVal( value );
            id = idTag;
            addTag( std::make_shared<xmlTag>("id",id) );
        }
        void appendWgt( std::shared_ptr<std::string> document ){
            if( !isWritten() ){ fullWriter(); }
            *document += *writtenSelf;
        }
        void appendWgt( std::string* document ){
            if( !isWritten() ){ fullWriter(); }
            *document += *writtenSelf;
        }
        std::shared_ptr<std::string> appendWgt( std::string_view document ){
            if(!isWritten() ){ fullWriter(); }
            auto retDoc = std::make_shared<std::string>( document );
            *retDoc += *writtenSelf;
            return retDoc;
        }
    protected:
        std::string_view comment;
        std::string valS;
        std::string id;
        double valD;
        void fullWriter() override {
            writtenSelf = std::make_shared<std::string>( "<wgt" );
            for( auto tag : tags ){
                *writtenSelf += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            *writtenSelf += ">" + std::string(valS) + "</wgt>\n";
            modded = false;
            written = true;
        }
    };

    // ZW: fcn for finding the next block in SLHA format
    // parameter cards
    size_t blockFinder( std::string_view parseFile, size_t startPt = 0 ){
        if( parseFile.size() > 5 ){ if( clStringComp(parseFile.substr(0,5), std::string("block") )){ return size_t(0); } }
        return clStringFind( parseFile, std::string("\nblock"), startPt );
    }

    // ZW: fcn for finding each decay line in SLHA format
    // parameter card
    std::vector<std::string_view> decBlockStractor( std::string_view parseFile ){
        auto allDs = nuFindEach( parseFile, "\nd" );
        std::vector<std::string_view> decLines;
        decLines.reserve( allDs->size() );
        for( auto pos : *allDs )
        {
            if( !(clStringComp(parseFile.substr( pos+1, 5 ), std::string("decay"))) ){ continue; }
            decLines.push_back( parseFile.substr( pos + 1, parseFile.find( "\n", pos + 1 ) - pos - 1 ) );
        }
        return decLines;
    }

    // ZW: fcn for extracting the relevant lines of
    // a block in SLHA format parameter card
    // removes any comments between start of this block and next
    // and also ignores lines with other information,
    // eg DECAY lines
    std::vector<std::string_view> blockLineStractor( std::string_view parseFile, size_t startPt = 0){
        auto blockStrt = blockFinder( parseFile, startPt );
        auto newBlock = blockFinder( parseFile, blockStrt + 1 );
        std::vector<std::string_view> paramLines;
        paramLines.reserve( nuStrCount( parseFile, "\n" ) );
        std::shared_ptr<std::vector<std::string_view>> parLines;
        if( newBlock == size_t(-1) ){ parLines = nuLineSplitter( parseFile.substr( blockStrt ) ); }
        else{ parLines = nuLineSplitter( parseFile.substr( blockStrt, newBlock - blockStrt ) ); }
        for( auto line : *parLines )
        {
            if( line.size() == 0 ){ continue; }
            if( line[0] != ' ' ){ continue; }
            paramLines.push_back( line );
        }
        return paramLines;
    }

    // ZW: struct for handling the first line of
    // LHE format event block
    struct evHead {
    public:
        std::string_view getComment(){ return comment; }
        std::string_view getWeight(){ return weight; }
        std::string_view getScale(){ return scale; }
        std::string_view getAQED(){ return aqed; }
        std::string_view getAQCD(){ return aqcd; }
        std::string_view getNprt(){ return nprt; }
        std::string_view getProcID(){ return procid; }
        bool isModded(){ return modded; }
        bool isWritten(){ return written; }
        void setComment( std::string_view nuCom ){ modded = true; comment = nuCom; }
        void setWeight( std::string_view nuWgt ){ modded = true; weight = nuWgt; }
        void setScale( std::string_view nuScale ){ modded = true; scale = nuScale; }
        void setAQED( std::string_view nuAQED ){ modded = true; aqed = nuAQED; }
        void setAQCD( std::string_view nuAQCD ){ modded = true; aqcd = nuAQCD; }
        void setNprt( std::string_view nuNprt ){ modded = true; nprt = nuNprt; }
        void setProcID( std::string_view nuProcID ){ modded = true; procid = nuProcID; }
        std::shared_ptr<std::string> getContent(){
            if( !isWritten() || isModded() ){ writer(); }
            return content;
        }
        evHead(){ return; }
        evHead( const std::string_view originFile, size_t beginLine = 0, size_t endLine = std::string_view::npos )
        {
            if( originFile.size() == 0){ return; }
            beginLine = originFile.find_first_not_of("\n ", beginLine);
            if( endLine == std::string_view::npos ){ endLine = originFile.find("\n", beginLine ) + 1; }
            sourceFile = originFile.substr( beginLine, endLine - beginLine );
            auto evLine = nuWordSplitter( sourceFile );
            nprt = evLine->at(0) ;
            procid = evLine->at(1);
            weight = evLine->at(2);
            scale = evLine->at(3);
            aqed = evLine->at(4);
            aqcd = evLine->at(5);
         }
    protected:
        std::shared_ptr<std::string> content;
        std::string_view sourceFile;
        std::string_view comment;
        std::string_view weight;
        std::string_view scale;
        std::string_view aqed;
        std::string_view aqcd;
        std::string_view nprt;
        std::string_view procid;
        bool modded = false;
        bool written = false;
        void writer(){
            if( isWritten() && !isModded() ){ return; }
            if( !isModded() ){ content = std::make_shared<std::string>( sourceFile ); return; }
            auto retText = std::make_shared<std::string>( " " );
            *content = " " + std::string( nprt );
            for( int k ; k < 8 - procid.length() ; ++k ){ *content += " "; }
            *content +=  std::string( procid ) + " " + std::string( weight ) + " " + std::string( scale ) + " " + std::string( aqed ) + " " + std::string( aqcd );
            if( comment != "" ){ *content += " # " + std::string( comment ); }
            *content += "\n";
            modded = false;
            written = true;
        }
    };

    // ZW: struct for handling particle lines
    // in LHE format event block
    struct lhePrt{
    public:
        std::string_view getLine(){ return sourceFile; }
        std::string_view getComment(){ return comment; }
        std::vector<std::string_view> getMom(){ return std::vector<std::string_view>( std::begin( mom ), std::end( mom ) ); }
        std::string_view getE(){ return energy; }
        std::string_view getMass(){ return mass; }
        std::string_view getVTim(){ return vtim; }
        std::string_view getSpin(){ return spin; }
        std::string_view getPDG(){ return pdg; }
        std::string_view getStatus(){ return status; }
        std::vector<std::string_view> getMothers(){ return std::vector<std::string_view>( std::begin( mothers ), std::end( mothers ) ); }
        std::vector<std::string_view> getColor(){ return std::vector<std::string_view>( std::begin( icol ), std::end( icol ) ); }
        void setComment( std::string_view nuCom ){ modded = true; comment = nuCom; }
        void setMom( std::vector<std::string_view> nuMom ){ modded = true; mom[0] = nuMom[0]; mom[1] = nuMom[1]; mom[2] = nuMom[2]; }
        void setEnergy( std::string_view nuE ){ modded = true; energy = nuE; }
        void setMass( std::string_view nuM ){ modded = true; mass = nuM; }
        void setVTim( std::string_view nuVTim ){ modded = true; vtim = nuVTim; }
        void setSpin( std::string_view nuSpin ){ modded = true; spin = nuSpin; }
        void setPDG( std::string_view nuPDG ){ modded = true; pdg = nuPDG; }
        void setStatus( std::string_view nuSt ){ modded = true; status = nuSt; }
        void setMothers( std::vector<std::string_view> nuMum ){ modded = true; mothers[0] = nuMum[0]; mothers[1] = nuMum[1]; }
        void setColors( std::vector<std::string_view> nuCol ){ modded = true; icol[0] = nuCol[0]; icol[1] = nuCol[1]; }
        bool isModded(){ return modded; }
        bool isWritten(){ return written; }
        std::shared_ptr<std::string> getContent(){
            if( !isWritten() || isModded() ){ writer(); }
            return content;
        }
        lhePrt(){ return; }
        lhePrt( const std::string_view originFile, const size_t& beginLine = 0, const size_t& endLine = std::string_view::npos )
        {
            sourceFile = originFile.substr( beginLine, endLine - beginLine );
            auto evLine = nuWordSplitter( sourceFile );
            pdg = evLine->at(0);
            status = evLine->at(1);
            mothers[0] = evLine->at(2); mothers[1] = evLine->at(3);
            icol[0] = evLine->at(4); icol[1] = evLine->at(5);
            for( int k = 6 ; k < 9 ; ++k){
                mom[k-6] = evLine->at(k);
            }
            energy = evLine->at(9);
            mass = evLine->at(10);
            vtim = evLine->at(11);
            spin = evLine->at(12);
            if( evLine->size() > 13 ){ comment = sourceFile.substr( sourceFile.find( "#" ) ); }
        }
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
        void writer(){
            if( isWritten() && !isModded() ){ return; }
            if( !isModded() ){ content = std::make_shared<std::string>( sourceFile ); return; }
            *content = "";
            for( int k = 0; k < 10 - pdg.length() ; ++k ){ *content += " "; }
            *content += std::string(pdg) + " " + std::string(status);
            for( auto mum : mothers ){ *content += "    " + std::string( mum ); }
            for( auto col : icol ){ *content += "  " + std::string( col ); }
            for( auto pval : mom ){ *content += " " + std::string(pval); }
            *content += " " + std::string( energy ) + " " + std::string( mass ) + " " + std::string( vtim ) + " " + std::string( spin );
            if( comment != "" ){ *content += " # " + std::string( comment ); }
            *content += "\n";
            modded = false;
            written = true;
        }
    };

    // ZW: struct for handling LHE format event block
    struct event : xmlNode {
    public:
        evHead getHead(){ return header; }
        std::vector<std::shared_ptr<lhePrt>> getPrts(){ return prts; }
        std::vector<std::shared_ptr<bodyWgt>> getWgts(){ return rwgt; }
        void setHead( evHead head ){ modded = true; header = head; }
        void addPrt( std::shared_ptr<lhePrt> prtcl ){ modded = true; prts.push_back( prtcl ); }
        void addPrt( lhePrt prtcl ){ modded = true; prts.push_back( std::make_shared<lhePrt>(prtcl) ); }
        void setPrts( std::vector<std::shared_ptr<lhePrt>> prtcls ){ modded = true; prts = prtcls; }
        void addWgt( bodyWgt nuWgt ){ addedWgt = true; rwgt.push_back( std::make_shared<bodyWgt>(nuWgt) ); }
        void addWgt( std::shared_ptr<bodyWgt> nuWgt ){ modded = true; rwgt.push_back( nuWgt ); }
        bool newWeight(){ return addedWgt; }
        int getNprt()
        {
            return prts.size();
        }
        bool isModded( bool deep ) override {
            bool modStat = modded;
            if( !deep ){ return modStat; }
            for( auto child : children ){ modStat = (modStat || child->isModded( deep )); }
            modStat = (modStat || header.isModded());
            for( auto prt : prts ){ modStat = (modStat || prt->isModded()); }
            for( auto wgt : rwgt ){ modStat = (modStat || wgt->isModded()); }
            return modStat;
        }
        event(){ return; }
        event( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} ) 
        : xmlNode(originFile, begin, childs) {
            xmlFile = originFile; start = begin; children = childs; size_t trueStart = originFile.find_first_not_of(" ", begin+1);
            if( trueStart != std::string_view::npos ){name = originFile.substr( trueStart, originFile.find_first_of(">/ ", trueStart) - trueStart );}
            auto vals = lineFinder( originFile.substr( trueStart, originFile.find("<", trueStart +  3 ) - trueStart + 3 ));
            header = evHead(originFile, vals->at(0) + trueStart, vals->at(1) + trueStart + 1 );
            prts.reserve(vals->size());
            for( int k = 1 ; k < std::stoi(std::string(header.getNprt())) + 1; ++k)
            {
                prts.push_back( std::make_shared<lhePrt>(originFile, vals->at(k) + trueStart + 1, vals->at(k+1) + trueStart + 1) );
            }
        }
        event( const xmlNode& originFile )
        : xmlNode( originFile ) {
            size_t trueStart = xmlFile.find_first_not_of(" ", start+1);
            auto vals = lineFinder( xmlFile.substr( trueStart, xmlFile.find("<", trueStart +  3 ) - trueStart + 3 ));
            header = evHead(xmlFile, vals->at(0) + trueStart, vals->at(1) + trueStart );
            prts.reserve(vals->size());
            for( int k = 1 ; k < std::stoi(std::string(header.getNprt())) + 1; ++k)
            {
                prts.push_back( std::make_shared<lhePrt>(xmlFile, vals->at(k) + trueStart + 1, vals->at(k+1) + trueStart) );
            }
        }
        bool prtsAreMod(){
            for( auto prt : prts ){ if( prt->isModded() ){ return true; } }
            return false;
        }
        bool headIsMod(){
            return header.isModded();
        }
    protected:
        std::vector<std::shared_ptr<bodyWgt>> rwgt;
        std::shared_ptr<xmlNode> childRwgt;
        bool hasRwgt(){
            if( rwgt.size() > 0 ){ return true; }
            return false;
        }
        bool rwgtChild(){
            if( childRwgt != nullptr ){ return true; }
            for( auto child : children ){ if( clStringComp(child->getName(), std::string("rwgt") ) ){ childRwgt = child; return true; } }
            return false;
        }
        bool bothRwgt(){ return (hasRwgt() && rwgtChild() ); }
        bool eitherRwgt(){ return (hasRwgt() || rwgtChild() ); }
        evHead header;
        std::vector<std::shared_ptr<lhePrt>> prts;
        bool inRwgtChild( std::string_view name ){ 
            for( auto child : childRwgt->getChildren() ){ 
                for( auto tag : child->getTags() ){ if(clStringComp(tag->getVal(), name)){ return true; } }
            }
            return false;
        }
        bool checkRwgtOverlap(){
            for( auto wgt : rwgt ){ 
                for( auto tag : wgt->getTags() ){ if( inRwgtChild( tag->getVal() ) ){ return true; } }
            }
            return false;
        }
        void childRwgtWriter(){
            if( rwgtChild() ){ nodeContent += *childRwgt->nodeWriter(); }
        }
        void vecRwgtWriter( bool midNode = false ){
            if( !midNode ){ nodeContent += "<rwgt>\n"; }
            for( auto wgt : rwgt ){ 
                nodeContent += *wgt->nodeWriter();
            }
            nodeContent += "</rwgt>\n";
        }
        void rwgtWriter(){
            if( bothRwgt() ){ if( checkRwgtOverlap() ){ childRwgtWriter(); return; } 
                childRwgtWriter();
                nodeContent.erase( nodeContent.size() - 8, 8 );
                vecRwgtWriter();
                return;
            } else {
                if( hasRwgt() ){ vecRwgtWriter(); return; }
                if( rwgtChild() ){ childRwgtWriter(); return; }
            }
        }
        void contWriter() override {
            nodeContent = "\n" + *header.getContent();
            for( auto prt : prts ){
                nodeContent += *prt->getContent();
            }
        }
        void childWriter() override {
            for( auto child : children ){
                if( clStringComp( child->getName(), std::string("wgt") ) ){ continue; }
                nodeContent += *child->nodeWriter();
            }
        }
        bool addedWgt = false;
        void fullWriter() override {
            if( isModded( false ) ){
                headWriter();
                contWriter();
                childWriter();
                rwgtWriter();
                endWriter();
                writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
                modded = false;
            } else if( !isWritten() ){
            writtenSelf = std::make_shared<std::string>( xmlFile.substr( start, end - start ) );
            written = true;
            }
        }
        void fullWriter( bool deep ){
            if( !deep ){ fullWriter(); return; }
            if( isModded( true ) ){
                headWriter();
                contWriter();
                childWriter();
                rwgtWriter();
                endWriter();
                writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
                modded = false;
                written = true;
            } else if( !isWritten() ){
            writtenSelf = std::make_shared<std::string>( xmlFile.substr( start, end - start ) );
            written = true;
            }
        }
        void appendWgts(){
            if( !addedWgt ){ return; }
            writtenSelf->erase( writtenSelf->size() - 17, 17 );
            for( auto wgt : rwgt ){
                if( !wgt->isWritten() ){ wgt->appendWgt( writtenSelf ); }
            }
            *writtenSelf += "</rwgt>\n</event>\n";
        }
    public:
        std::shared_ptr<std::string> nodeWriter() override {
            if( isModded(false) || !isWritten() ){ fullWriter(); return writtenSelf; }
            if( addedWgt ){ appendWgts(); }
            return writtenSelf;
        }
        std::shared_ptr<std::string> nodeWriter( bool recursive ){
            if( isModded( recursive ) || !isWritten() ){ fullWriter(); return writtenSelf; }
            if( addedWgt ){ appendWgts(); }
            return writtenSelf;
        }
    };

    // ZW: struct for handling the first line of
    // LHE format init tag
    struct lheInitHead{
    public:
        std::string_view idbmup[2];
        std::string_view ebmup[2];
        std::string_view pdfgup[2];
        std::string_view pdfsup[2];
        std::string_view idwtup;
        std::string_view nprup;
        bool isWritten(){ return written; }
        bool isModded(){ return modded; }
        std::shared_ptr<std::string> getContent(){ 
            if( isModded() || !isWritten() ){ writer(); }
            return content; }
        lheInitHead( std::string_view initHead ){
            auto vals = *nuBlankSplitter( initHead );
            if( vals.size() < 10 ){ return; }
            idbmup[0] = vals[0]; idbmup[1] = vals[1];
            ebmup[0] = vals[2]; ebmup[1] = vals[3];
            pdfgup[0] = vals[4]; pdfgup[1] = vals[5];
            pdfsup[0] = vals[6]; pdfsup[1] = vals[7];
            idwtup = vals[8]; nprup = vals[9];
        }
        lheInitHead( xmlNode& initNode )
        {
            if( initNode.getName() != "init" ){ return; }
            auto startPos = initNode.getFile().find( ">", initNode.getStart() ) + 1;
            auto endPos = initNode.getFile().find( "\n", startPos );
            auto vals = *nuBlankSplitter( initNode.getFile().substr( startPos, endPos - startPos ) );
            idbmup[0] = vals[0]; idbmup[1] = vals[1];
            ebmup[0] = vals[2]; ebmup[1] = vals[3];
            pdfgup[0] = vals[4]; pdfgup[1] = vals[5];
            pdfsup[0] = vals[6]; pdfsup[1] = vals[7];
            idwtup = vals[8]; nprup = vals[9];
        }
    protected:
        std::shared_ptr<std::string> content;
        bool written = false;
        bool modded = false;
        void writer(){
            *content = std::string(idbmup[0]) + " " + std::string(idbmup[1]) + " " + std::string(ebmup[0]) + " " + std::string(ebmup[1]) + " " + std::string(pdfgup[0]) 
    + " " + std::string(pdfgup[1]) + " " + std::string(pdfsup[0]) + " " + std::string(pdfsup[1]) + " " + std::string(idwtup) + " " + std::string(nprup) +"\n";
            written = true;
            modded = false;
        }
    };

    // ZW: struct for handling process lines
    // in LHE format init tag
    struct lheInitLine {
    public:
        std::string_view xsecup;
        std::string_view xerrup;
        std::string_view xmaxup;
        std::string_view lprup;
        bool isWritten(){ return written; }
        bool isModded(){ return modded; }
        std::shared_ptr<std::string> getContent(){ 
            if( isModded() || !isWritten() ){ writer(); }
            return content; }
        lheInitLine(){}
        lheInitLine( std::string_view procLine )
        {
            auto vals = *nuBlankSplitter( procLine );
            if( vals.size() < 4 ){ return; }
            xsecup = vals[0];
            xerrup = vals[1];
            xmaxup = vals[2];
            lprup = vals[3];
        }
    protected:
        std::shared_ptr<std::string> content;
        bool written = false;
        bool modded = false;
        void writer(){
            *content = std::string(xsecup) + " " + std::string(xerrup) + " " + std::string(xmaxup) + " " + std::string(lprup) + "\n";
            written = true;
            modded = false;
        }
    };

    // ZW: struct for handling single parameter line in
    // SLHA format parameter card
    struct paramVal{
    public:
        double value = 0;
        int id = 0;
        std::string_view realLine;
        std::string_view comment;
        std::string_view idStr;
        std::string_view valStr;
        virtual void parse(){
            id = std::stoi( std::string(idStr) );
            value = std::stod( std::string(valStr) );
        }
        paramVal(){ realLine = ""; idStr = ""; valStr = ""; }
        paramVal( std::string_view paramLine, bool parseOnline = false )
        {
            if( paramLine.find("\n") != std::string_view::npos ){
                auto startPos = paramLine.find_first_not_of(" \n", paramLine.find("\n"));
                if( startPos!= std::string_view::npos ){
                auto endPos = paramLine.find("\n", startPos);
                realLine = paramLine.substr(startPos, endPos - startPos - 1);
                } else{
                    realLine = paramLine.substr( 0, paramLine.find("\n") - 1 );
                }
            }
            realLine = paramLine;
            auto vals = *nuBlankSplitter( realLine );
            idStr = vals[0];
            valStr = vals[1];
            if( parseOnline ){ 
            if( vals.size() > 2 )
            {
                auto comStart = realLine.find("#");
                comStart = realLine.find_first_not_of( " #", comStart );
                comment = realLine.substr( comStart, realLine.find("\n", comStart) - comStart );
            }
            parse(); }
        }
        bool isMod(){ return modded; }
        bool modded = false;
        virtual std::shared_ptr<std::string> selfWrite(){
            auto writeVal = std::make_shared<std::string>("");
            if( isMod() )
            {
                for( int k = idStr.size() ; k < 5 ; ++k ){ *writeVal += " "; }
                *writeVal += std::string( idStr ) + " " + std::string( valStr );
                if( comment.size() != 0 ){
                    *writeVal += " # " + std::string( comment );
                }
                *writeVal += "\n";
            }
            else{ *writeVal = std::string( realLine ) + "\n"; }
            return writeVal;
        }
    };

    // ZW: struct for handling single DECAY line
    // in SLHA format parameter card
    struct decVal : paramVal{
    public:
        void parse() override {
            auto vals = *nuBlankSplitter( realLine );
            id = std::stoi( std::string(vals[1]) );
            value = std::stod( std::string(vals[2]) );
            if( vals.size() > 3 )
            {
                auto comStart = realLine.find("#");
                comment = realLine.substr( comStart, realLine.find("\n", comStart) - comStart );
            }
        }
        decVal( std::string_view paramLine = "", bool parseOnline = false ) : paramVal( paramLine, false )
        {
            if( parseOnline ){ parse(); }
        }
        std::shared_ptr<std::string> selfWrite() override {
            auto writeVal = std::make_shared<std::string>("");
            if( isMod() )
            {
                *writeVal += "DECAY " + std::string( idStr ) + " " + std::string( valStr );
                if( comment.size() != 0 ){
                    *writeVal += " # " + std::string( comment );
                }
                *writeVal += "\n";
            }
            else{ *writeVal = std::string( realLine ) + "\n"; }
            return writeVal;
        }
    };

    // ZW: struct for handling parameter block
    // in SLHA format parameter card
    struct paramBlock {
    public:
        std::string_view realBlock;
        size_t startPt;
        std::string_view comment;
        std::string_view initComm;
        std::string_view name;
        std::vector<paramVal> params;
        virtual void parse( bool parseOnline = false ){
            if( realBlock.size() == 0 ){ return; }
            if( !(clStringComp(realBlock.substr(startPt+1, 5), std::string("block"))) ){ startPt = clStringFind( realBlock, std::string("\nblock") ); }
            auto namePt = realBlock.find_first_not_of( " ", startPt + 7 );
            name = realBlock.substr( namePt, realBlock.find_first_of( " \n", namePt ) - namePt );
            if( realBlock.find( " ", namePt ) < realBlock.find( "\n", namePt ) )
            {comment = realBlock.substr( namePt + name.size(), realBlock.find( "\n", namePt ) - namePt - name.size() ); }
            auto paramLines = blockLineStractor( realBlock.substr( startPt ) );
            params.reserve( paramLines.size() );
            for( auto line : paramLines )
            {
                params.push_back( paramVal( line, parseOnline ) );
            }
        }
        paramBlock(){ return; }
        paramBlock( std::string_view paramSet, bool parseOnline = false )
        {
            realBlock = paramSet;
            startPt = clStringFind( realBlock, std::string("\nB") );
            if( parseOnline ){ parse(parseOnline);  }
        }
        bool isMod(){ return modded; }
        bool modded = false;
        virtual std::shared_ptr<std::string> selfWrite(){
            auto writeBlock = std::make_shared<std::string>("");
            if( isMod() )
            {
                *writeBlock += "\nBLOCK " + std::string(name);
                if( comment.size() > 0 ){
                    *writeBlock += " # " + std::string( comment );
                }
                *writeBlock += "\n";
                for ( auto val : params )
                {
                    *writeBlock += *val.selfWrite();
                } 
            }
            else{ if( startPt == size_t(-1) ){
                *writeBlock += realBlock;
            } else {
                *writeBlock = realBlock.substr( startPt );
            } }
            return writeBlock;
        }
    };

    // ZW: struct for handling DECAY lines
    // in SLHA format parameter card
    struct decBlock : paramBlock {
    public:
        std::vector<decVal> decays;
        void parse( bool parseOnline = false ) override{
        if( realBlock.size() == 0 ){ return; }
            auto decLines = clFindEach( realBlock, std::string("\ndecay") );
            decays.reserve(decLines->size());
            if( realBlock.size() > 5 ){  if( clStringComp( realBlock.substr(0,5), std::string("decay")) )
            { decays.push_back( decVal(realBlock.substr( 0, realBlock.find("\n") ), parseOnline) ); } }
            for( auto pts : *decLines )
            {
                auto lineBr = realBlock.find( "\n", pts + 1 );
                if( lineBr == size_t(-1) ){ decays.push_back( decVal( realBlock.substr( pts + 1), parseOnline ) ); continue; }
                decays.push_back( decVal( realBlock.substr( pts + 1, lineBr - pts - 1 ), parseOnline ) );
            }
        }
        void parse( std::shared_ptr<std::vector<size_t>> decLines, bool parseOnline = false ) {
            decays.reserve(decLines->size());
            if( realBlock.size() > 5 ){ if( clStringComp( realBlock.substr(0,5), std::string("decay")) )
            { decays.push_back( decVal(realBlock.substr( 0, realBlock.find("\n") ), parseOnline) ); } }
            for( auto pts : *decLines )
            {
                auto lineBr = realBlock.find( "\n", pts + 1 );
                if( lineBr == size_t(-1) ){ decays.push_back( decVal( realBlock.substr( pts + 1), parseOnline ) ); continue; }
                decays.push_back( decVal( realBlock.substr( pts + 1, lineBr - pts - 1 ), parseOnline ) );
            }
        }
        decBlock( std::string_view paramSet = "", bool parseOnline = false ) : paramBlock( paramSet, parseOnline )
        {
            realBlock = paramSet;
            if( parseOnline ){ parse(parseOnline); }
        }
        std::shared_ptr<std::string> selfWrite() override {
            auto writeBlock = std::make_shared<std::string>("");
            *writeBlock += "\n";
            for ( auto val : decays )
            {
                *writeBlock += *val.selfWrite(); 
            }
            return writeBlock;
        }
    };

    // ZW: struct for handling SLHA parameter cards
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
        void parse( bool parseOnline = false )
        {
            if( parsed ){ return; }
            if( xmlFile.substr(start,1).find_first_of("BbDd#") == std::string_view::npos ){ start = clStringFindIf( xmlFile, std::string("\n"), lambdaNu ); }
            auto blockPts = clFindEach( xmlFile, std::string("\nblock") );
            auto decLines = clFindEach( xmlFile, std::string("\ndecay") );
            header = xmlFile.substr( start, std::min( blockPts->at(0), decLines->at(0) ) - start );
            auto startPt = blockStart;
            for( int k  = 0 ; k < blockPts->size() - 1 ; ++k )
            {
                blocks.push_back( paramBlock( xmlFile.substr( blockPts->at(k), blockPts->at(k+1) - blockPts->at(k) ), parseOnline ) );
            }
            blocks.push_back(paramBlock(xmlFile.substr(blockPts->at(blockPts->size()-1), clStringFindIf( xmlFile, std::string("\n"),
             lambda, blockPts->at(blockPts->size()-1) + 1) - blockPts->at(blockPts->size()-1)), parseOnline));
            decays = decBlock( xmlFile );
            decays.parse( decLines, parseOnline );
            parsed = true;
        } 
        lesHouchesCard( const std::string_view originFile = "", const size_t& begin = 0, bool parseOnline = false ){ 
            xmlFile = originFile; start = begin; size_t trueStart = originFile.find_first_not_of("\n ", begin+1);
            modded = false; blockStart = clStringFindIf( xmlFile, std::string("\n"), lambda, start + 1); end = xmlFile.find("</", blockStart);
            parsed = false;
            if( parseOnline ){ parse( parseOnline ); }
        }
        bool isMod(){ return modded; }
        std::shared_ptr<std::string> selfWrite(){
            auto writeCard = std::make_shared<std::string>(header);
            if( isMod() )
            { for( auto block : blocks )
                { *writeCard += *block.selfWrite(); }
                *writeCard += *decays.selfWrite(); }
            else{
                if( end != size_t(-1) ){ *writeCard += std::string( xmlFile.substr( blockStart, end - blockStart ) );
                } else{ *writeCard += std::string( xmlFile.substr( blockStart ) ); }
            }
            return writeCard;
        }
    };

    struct slhaNode : xmlNode {
    public:
        std::shared_ptr<lesHouchesCard> getParameters(){
            modded = true;
            return parameterCard;
        }
        slhaNode() : xmlNode(){}
        slhaNode( lesHouchesCard parameters ) : xmlNode(){
            parameterCard = std::make_shared<lesHouchesCard>( parameters );
            pCardInit = true;
        } 
        slhaNode( std::shared_ptr<lesHouchesCard> parameters ) : xmlNode(){
            parameterCard = parameters;
            pCardInit = true;
        }
        slhaNode( xmlNode& node, bool parseOnline = false ) : xmlNode( node ){ 
            parameterCard = std::make_shared<lesHouchesCard>( node.getFile(), node.getStart(), parseOnline );
        }
        slhaNode( xmlNode* node, bool parseOnline = false ) : xmlNode( *node ){ 
            parameterCard = std::make_shared<lesHouchesCard>( node->getFile(), node->getStart(), parseOnline );
        }
        slhaNode( std::shared_ptr<xmlNode> node, bool parseOnline = false ) : xmlNode( *node ){ 
            parameterCard = std::make_shared<lesHouchesCard>( node->getFile(), node->getStart(), parseOnline );
        }
        slhaNode( const std::string_view originFile, const size_t& begin = 0, bool parseOnline = false )
        : xmlNode( originFile, begin ){
            if( parse() ){ parameterCard = std::make_shared<lesHouchesCard>( content, begin, parseOnline ); pCardInit = true; }
        }
    protected:
        std::shared_ptr<lesHouchesCard> parameterCard;
        bool pCardInit = false;
        void headWriter() override{
            nodeHeader = "<slha";
            for( auto tag : tags ){
                nodeHeader += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            nodeHeader += ">";
        }
        void endWriter() override{ nodeEnd += "</slha>\n"; }
        void contWriter() override{
            if( pCardInit ){
                nodeContent = *parameterCard->selfWrite();
            } else {
                nodeContent = content;
            }
        }
    };

    // ZW: struct for handling LHE init nodes
    struct initNode : xmlNode {
    public:
        std::shared_ptr<lheInitHead> getHead(){ return initHead; }
        std::vector<std::shared_ptr<lheInitLine>> getLines(){ return initLines; }
        void setHead( std::shared_ptr<lheInitHead> head ){ modded = true; initHead = head; }
        void setLines( std::vector<std::shared_ptr<lheInitLine>> lines ){ modded = true; initLines = lines; initHead->nprup = std::to_string( initLines.size() ); }
        void addLine( std::shared_ptr<lheInitLine> line ){ modded = true; initLines.push_back( line ); initHead->nprup = std::to_string( initLines.size() ); }
        initNode() : xmlNode(){ name = "init"; }
        initNode( const std::string_view originFile, const size_t& begin = 0, bool parseOnline = false )
        : xmlNode( originFile, begin ){
            auto strtPt = originFile.find_first_not_of(" \n", originFile.find(">", start+1));
            content = originFile.substr( strtPt, originFile.find("</", strtPt) - strtPt );
        }
    protected:
        std::shared_ptr<lheInitHead> initHead;
        std::vector<std::shared_ptr<lheInitLine>> initLines;
        bool parseContent() override{
            if( content.size() == 0 ){ return false; }
            auto linebreaks = lineFinder( content );
            if( linebreaks->size() == 0 ){ return false; }
            initHead = std::make_shared<lheInitHead>(content.substr( 0, linebreaks->at(0) ) );
            for( int k = 0 ; k < linebreaks->size() - 1 ; ++k ){
                initLines.push_back( std::make_shared<lheInitLine>( content.substr( linebreaks->at(k), linebreaks->at(k+1) - linebreaks->at(k) ) ) );
            }
            return true;
        }
        void contWriter() override{
            if( isModded() ){nodeContent = std::string( content ); return; }
            nodeContent = *initHead->getContent();
            for( auto line : initLines ){
                nodeContent += *line->getContent();
            }
        }
    };
    
    // ZW: struct for explicitly handling LHE header nodes
    struct lheHead : xmlNode {
    public:
        size_t addWgtGroup( std::shared_ptr<weightGroup>& wgtGroup ){
            hasRwgt = true; 
            modded = true; 
            if( wgtGrpInit( wgtGroup ) ){
                rwgtNodes->addGroup( wgtGroup );
            }
            return (rwgtNodes->noGrps() - 1);
        }
        size_t addWgtGroup( weightGroup wgtGroup ){
            hasRwgt = true;
            modded = true;
            auto wgtGrpPtr = std::make_shared<weightGroup>( wgtGroup );
            if( wgtGrpInit( wgtGrpPtr ) ){
                rwgtNodes->addGroup( std::make_shared<weightGroup>( wgtGroup ) );
            }
            return (rwgtNodes->noGrps() - 1);
        }
        void addWgt( unsigned int index, std::shared_ptr<headWeight> nuWgt ){
            if( index >= (int)rwgtNodes->getGroups().size() )
                throw std::runtime_error( "Appending weight to uninitialised weightgroup." );
            hasRwgt = true;
            modded = true;
            rwgtNodes->addWgt( index, nuWgt );
        }
        void addWgt( unsigned int index, headWeight nuWgt ){
            if( index >= (int)rwgtNodes->getGroups().size() )
                throw std::runtime_error( "Appending weight to uninitialised weightgroup." );
            hasRwgt = true;
            modded = true;
            rwgtNodes->addWgt( index, nuWgt );
        }
        void setInitRwgt( initRwgt initWgt ){  hasRwgt = true; modded = true; rwgtNodes = std::make_shared<initRwgt>(initWgt); }
        void setInitRwgt( std::shared_ptr<initRwgt> initWgt ){ hasRwgt = true; modded = true; rwgtNodes = initWgt; }
        std::vector<std::shared_ptr<weightGroup>> getWgtGroups(){ return rwgtNodes->getGroups(); }
        std::shared_ptr<initRwgt> getInitRwgt(){ return rwgtNodes; }
        std::shared_ptr<slhaNode> getParameters(){ return parameters; }
        void setParameters( std::shared_ptr<slhaNode> params ){ parameters = params; }
        bool rwgtInc(){ return hasRwgt; }
        lheHead(){ return; }
        lheHead( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} )
        : xmlNode(originFile, begin, childs){
            xmlFile = originFile; start = begin; children = childs; size_t trueStart = originFile.find_first_not_of(" ", begin+1);
            if( trueStart != std::string_view::npos ){name = originFile.substr( trueStart, originFile.find_first_of(">/ ", trueStart) - trueStart );}
        }
    protected:
        bool wgtGrpIsInit = false;
        bool wgtGrpInit( std::shared_ptr<weightGroup>& wgtGrp ){
            if( wgtGrpIsInit ){ return true; }
            if( rwgtNodes == nullptr ){
                rwgtNodes = std::make_shared<initRwgt>();
                wgtGrpIsInit = true;
                rwgtNodes->addGroup( wgtGrp );
                return false;
            } else throw std::runtime_error( "Error while initiating return LHE file header (initrwgt node is defined in an unrecognised manner)." );
        }
        std::shared_ptr<slhaNode> parameters;
        bool hasRwgt = false;
        std::shared_ptr<initRwgt> rwgtNodes;
        std::vector<std::shared_ptr<weightGroup>> initrwgt;
        bool relChildSet = false;
        std::vector<int> relChild;
        void setRelChild(){
            if( relChildSet ){ return; }
            relChild.reserve( children.size() );
            for( int k = 0 ; k < children.size() ; ++k ){
                auto child = &children[k];
                if( (*child)->getName() == "slha" ){ continue; }
                if( (*child)->getName() == "initrwgt" ){ continue; }
                relChild.push_back( k );
            }
            relChildSet = true;
        }
        bool parseChildren( bool recursive ){
            bool status = true;
            for( auto child : children ){
                if( child->getName() == "slha" || child->getName() == "initrwgt" ){ continue; }
                child->parser( recursive );
                status = (status && child->isParsed() );
                deepParsed = true;
            }
            return status;
        }
        void headWriter() override{
            nodeHeader =  "<header";
            for( auto tag : tags ){
                nodeHeader += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            nodeHeader += ">\n";
        }
        void childWriter() override{
            setRelChild();
            for( auto relKid : relChild ){
                nodeContent += *(children[relKid]->nodeWriter());
            }
            if( parameters != nullptr ){ nodeContent += *parameters->nodeWriter(); }
            if( hasRwgt ){ 
                nodeContent += *rwgtNodes->nodeWriter();
            }
        }
        void fullWriter() override{
            if( isModded() ){
            headWriter();
            contWriter();
            childWriter();
            endWriter();
            writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
            written = true;
            }
        }
    };

    // ZW: struct for keeping track of appended weights in LHE node,
    // since weight information is stored both in the header 
    // and in the individual events
    struct newWgt{
    protected:
        std::shared_ptr<headWeight> headWgt;
        std::vector<std::shared_ptr<bodyWgt>> bodyWgts;
    public:
        newWgt( std::shared_ptr<headWeight> heaWgt, std::vector<std::shared_ptr<bodyWgt>> bodWgts ){
            headWgt = heaWgt; bodyWgts = bodWgts;
        }
        newWgt( std::shared_ptr<headWeight> heaWgt, std::shared_ptr<std::vector<double>> wgts ){
            headWgt = heaWgt;
            bodyWgts = std::vector<std::shared_ptr<bodyWgt>>(wgts->size());
            auto idTag = std::string(headWgt->getTag());
            if( idTag != "" ){
                for( size_t i = 0 ; i < wgts->size() ; ++i ){
                    bodyWgts[i] = std::make_shared<bodyWgt>(wgts->at(i), idTag);
                }
            } else{
                for( size_t i = 0 ; i < wgts->size() ; ++i ){
                    bodyWgts[i] = std::make_shared<bodyWgt>(wgts->at(i));
                }
            }
        }
        newWgt( std::string_view parameters, std::shared_ptr<std::vector<double>> wgts, std::string idTag = "pep_rwgt" ){
            headWgt = std::make_shared<headWeight>(parameters, idTag);
            bodyWgts = std::vector<std::shared_ptr<bodyWgt>>(wgts->size());
            for( size_t i = 0 ; i < wgts->size() ; ++i ){
                bodyWgts[i] = std::make_shared<bodyWgt>(wgts->at(i), idTag);
            }
        }
        newWgt( std::string_view parameters, int idNum, std::shared_ptr<std::vector<double>> wgts, std::string idTag = "pep_rwgt" ){
            std::string newTag = std::string( idTag ) + "_" + std::to_string( idNum );
            headWgt = std::make_shared<headWeight>(parameters, newTag);
            bodyWgts = std::vector<std::shared_ptr<bodyWgt>>(wgts->size());
            for( size_t i = 0 ; i < wgts->size() ; ++i ){
                bodyWgts[i] = std::make_shared<bodyWgt>(wgts->at(i), newTag);
            }
        }
        newWgt( std::string& parameters ){
            headWgt = std::make_shared<headWeight>(parameters);
        }
        newWgt( std::string& parameters, std::string& idTag ){
            headWgt = std::make_shared<headWeight>(parameters, idTag);
        }
        std::shared_ptr<headWeight> getHeadWgt(){ return headWgt; }
        std::vector<std::shared_ptr<bodyWgt>> getBodyWgts(){ return bodyWgts; }
        void addBdyWgts( std::shared_ptr<std::vector<double>> wgts ){
            auto idTag = std::string(headWgt->getTag());
            if( idTag != "" ){
                for( size_t i = 0 ; i < wgts->size() ; ++i ){
                    bodyWgts[i] = std::make_shared<bodyWgt>(wgts->at(i), idTag);
                }
            } else{
                for( size_t i = 0 ; i < wgts->size() ; ++i ){
                    bodyWgts[i] = std::make_shared<bodyWgt>(wgts->at(i));
                }
            }
        }
    };

    // ZW: general struct for handling LHE files explicitly
    struct lheNode : xmlNode {
    public:
        std::vector<std::shared_ptr<event>> events = {};
        std::shared_ptr<lheHead> header =  std::make_shared<lheHead>(xmlFile, start);
        std::shared_ptr<initNode> init = std::make_shared<initNode>(xmlFile, start);
        lheNode() : xmlNode(){}
        lheNode( const std::string_view originFile, const size_t& begin = 0, const std::vector<std::shared_ptr<xmlNode>>& childs = {} )
        : xmlNode(originFile, begin, childs){
            xmlFile = originFile; start = begin; children = childs; size_t trueStart = originFile.find_first_not_of(" ", begin+1);
            if( trueStart != std::string_view::npos ){name = originFile.substr( trueStart, originFile.find_first_of(">/ ", trueStart) - trueStart );}
        }
        bool isModded() override{ return modded; }
        bool isModded( bool deep ) override{
            if( !deep ){ return isModded(); }
            bool modStat = isModded();
            for( auto child : children ){ modStat = ( modStat || child->isModded( deep ) ); }
            for( auto event : events ){ modStat = ( modStat || event->isModded( deep ) ); }
            return modStat;
        }
        void addWgt( size_t index, newWgt& addedWgt ){
            header->addWgt( index, addedWgt.getHeadWgt() );
            auto wgtsVec = addedWgt.getBodyWgts();
            for( int k = 0 ; k < wgtsVec.size() ; ++k ){
                events[k]->addWgt( wgtsVec[k] );
            }
        }
    protected:
        virtual void headerWriter(){
            nodeContent += "\n" + *header->nodeWriter();
        }
        virtual void initWriter(){
            nodeContent += *init->nodeWriter();
        }
        virtual void eventWriter(){
            for( auto event : events ){
                nodeContent += *event->nodeWriter();
            }
        }
        void contWriter() override{
            nodeContent = "";
            headerWriter();
            initWriter();
            eventWriter();
        }
        void fullWriter() override{
            if( isModded( true ) ){
            headWriter();
            contWriter();
            endWriter();
            writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
            written = true;
            modded = false;
            } else if( !isWritten() ){
                writtenSelf = std::make_shared<std::string>( xmlFile.substr(start, end - start ) );
                written = true;
            }
        }
    public:    
        virtual std::shared_ptr<std::string> nodeWriter() {
            if( isModded( true ) || !isWritten() ){ fullWriter(); }
            return writtenSelf;
        }
    };

    // ZW: function for extracting event information from
    // LHE files
    std::vector<std::shared_ptr<std::vector<double>>> valExtraction( const lheNode& lheFile )
    {
        bool getGs = true;
        auto momVec = std::make_shared<std::vector<double>>();
        auto wgtVec = std::make_shared<std::vector<double>>();
        auto gVec = std::make_shared<std::vector<double>>();
        momVec->reserve( lheFile.events.size() * 4 * std::stoi(std::string(lheFile.events[0]->getHead().getNprt())) );
        wgtVec->reserve( lheFile.events.size() );
        gVec->reserve( lheFile.events.size() );
        if( getGs ){
        for( auto event : lheFile.events )
        {
            wgtVec->push_back(std::stod(std::string( event->getHead().getWeight() )));
            gVec->push_back( std::sqrt( 4.0 * M_PI * std::stod(std::string( event->getHead().getAQCD() ))));
            for( auto prt : event->getPrts() )
            {
                momVec->push_back(std::stod(std::string(prt->getE())));
                for( int p = 0 ; p < 3 ; ++p )
                { momVec->push_back(std::stod(std::string(prt->getMom()[p]))); }
            }
        }
        } else{
        for( auto event : lheFile.events )
        {
            wgtVec->push_back(std::stod(std::string( event->getHead().getWeight() )));
            gVec->push_back( std::stod(std::string( event->getHead().getAQCD() )));
            for( auto prt : event->getPrts() )
            {
                momVec->push_back(std::stod(std::string(prt->getE())));
                for( int p = 0 ; p < 3 ; ++p )
                { momVec->push_back(std::stod(std::string(prt->getMom()[p]))); }
            }
            
        } }
        return {momVec, gVec, wgtVec};
    }

    // ZW: fcn for parsing an LHE format event block
    // and return a PEP format event object
    std::shared_ptr<event> evPtrParsor( std::string_view parseFile, size_t& initPos, size_t& endPos )
    {
        auto currNode = std::make_shared<event>(parseFile, initPos);
        initPos = *nodeStartFind( parseFile, initPos + 1 );
        while( initPos < endPos )
        {
            currNode->addChild(xmlPtrParser( parseFile, initPos, endPos ));
        }
        size_t equalSign = parseFile.find("=", initPos);
        size_t nodeInitEnd = parseFile.find(">", initPos);
        while( equalSign < nodeInitEnd ){
            currNode->addTag( xmlTagParser(parseFile, equalSign) );
        }
        initPos = *nodeStartFind( parseFile, endPos );
        endPos = *nodeEndFind( parseFile, endPos + 1 );
        return currNode;
    }

    // ZW: fcn for parsing an LHE format header
    // and return a PEP format lheHead object
    std::shared_ptr<lheHead> lheHeadParser( std::string_view parseFile, size_t& initPos, size_t& endPos )
    {
        auto currNode = std::make_shared<lheHead>(parseFile, initPos);
        initPos = *nodeStartFind( parseFile, initPos + 1 );
        while( initPos < endPos )
        {
            auto nuStrtPos =  *nodeStartFind( parseFile, initPos);
            currNode->addChild(xmlPtrParser( parseFile, initPos, endPos ));
            if( currNode->getChildren()[ currNode->getChildren().size() - 1 ]->getName() == "init" ){ continue; }
            if( currNode->getChildren()[ currNode->getChildren().size() - 1 ]->getName() == "slha" ){
                auto nuLine = parseFile.find("\n", parseFile.find("<", initPos));
                currNode->setParameters( std::make_shared<slhaNode>(currNode->getChildren()[ currNode->getChildren().size() - 1 ]) );
            }
            if( currNode->getChildren()[ currNode->getChildren().size() - 1 ]->getName() == "initrwgt" ){
                currNode->setInitRwgt( std::make_shared<initRwgt>( currNode->getChildren()[ currNode->getChildren().size() - 1 ] ) );
            }
        }
        size_t equalSign = parseFile.find("=", initPos);
        size_t nodeInitEnd = parseFile.find(">", initPos);
        while( equalSign < nodeInitEnd ){
            currNode->addTag( xmlTagParser(parseFile, equalSign) );
        }
        initPos = *nodeStartFind( parseFile, endPos );
        endPos = *nodeEndFind( parseFile, endPos + 1 );
        return currNode;
    }

    // ZW: fcn for parsing an LHE format file
    // and return a PEP format LHE node object
    std::shared_ptr<lheNode> lheParser( std::string_view parseFile, size_t& initPos, size_t& endPos )
    {
        auto currNode = std::make_shared<lheNode>(parseFile, initPos);
        initPos = *nodeStartFind( parseFile, initPos + 1 );
        while( initPos < endPos )
        {
            auto nuStrtPos =  *nodeStartFind( parseFile, initPos);
            if( nuStrtPos == parseFile.find("<event", initPos) ){
                currNode->events.push_back( evPtrParsor( parseFile, initPos, endPos ) );
                continue;
            } else if( nuStrtPos == parseFile.find("<header", initPos) ){
                currNode->header = lheHeadParser( parseFile, initPos, endPos );
                continue;
            } else if( nuStrtPos == parseFile.find("<init", initPos) ){
                currNode->init = std::make_shared<initNode>( parseFile, initPos );
                initPos = *nodeStartFind( parseFile, endPos );
                endPos = *nodeEndFind( parseFile, *nodeEndFind( parseFile, endPos + 1 ) + 1);
                continue;
            } else {
            currNode->addChild(xmlPtrParser( parseFile, initPos, endPos ));
            }
        }
        size_t equalSign = parseFile.find("=", initPos);
        size_t nodeInitEnd = parseFile.find(">", initPos);
        while( equalSign < nodeInitEnd ){
            currNode->addTag( xmlTagParser(parseFile, equalSign) );
        }
        initPos = *nodeStartFind( parseFile, endPos );
        endPos = *nodeEndFind( parseFile, endPos + 1 );
        return currNode;
    }

    // ZW: struct for treating individual HEP
    // processes, formatted based on PDG codes
    // and the LHE particle status standard 
    struct lheProc {
    public:
        std::vector<std::string_view> minusOne;
        std::vector<std::string_view> plusOne;
        std::vector<std::string_view> minusTwo;
        std::vector<std::string_view> plusTwo;
        std::vector<std::string_view> plusThree;
        std::vector<std::string_view> minusNine;
        std::map<std::string_view,std::vector<std::string_view>> valVecs{{"-1", minusOne}, {"1", plusOne}, {"-2", minusTwo}, {"2", plusTwo}, {"3", plusThree}, {"-9", minusNine}};
        lheProc( event& eventNode )
        {
            for( auto prt : eventNode.getPrts() )
            {
                valVecs[prt->getStatus()].push_back(prt->getPDG());
            }
        }
    };

    // ZW: fcn for uploading text files
    // to the program, pushing all characters to lowercase
    std::shared_ptr<std::string> filePuller( const std::string& fileLoc )
    {
        std::ifstream fileLoad( fileLoc );
        std::stringstream buffer;
        buffer << fileLoad.rdbuf();
        auto fileContent = std::make_shared<std::string>(buffer.str());
        //std::transform( fileContent->begin(), fileContent->end(), fileContent->begin(), ::tolower );
        buffer.str(std::string());
        fileLoad.close();
        return fileContent;
    }

    // ZW: fcn for saving std::string to disk
    bool filePusher( std::string fileLoc, std::string fileCont )
    {
        std::ofstream fileWrite( fileLoc );
        if(!fileWrite){return false;}
        fileWrite << fileCont;
        fileWrite.close();
        return true;
    }

    // ZW: fcn for extracting the fill
    // process information from an LHE event
    std::shared_ptr<std::map<std::string_view, std::vector<std::string_view>>> pgdXtract( event& currEv, const std::vector<std::string>& pdgVec )
    {
        auto currProc = std::make_shared<std::map<std::string_view, std::vector<std::string_view>>>();
        auto &useProc = *currProc;
        for( auto prt : currEv.getPrts() )
        {
            useProc[ prt->getStatus() ].push_back(prt->getPDG());
        }
        return currProc;
    }

    // ZW: fcn for comparing two processes it the
    // format output by pgdXtract
    bool sameProcString( std::map<std::string_view, std::vector<std::string_view>>& firstVec, std::map<std::string_view,
     std::vector<std::string_view>>& secVec, const std::vector<std::string>& pdgVec )
    {
        if( firstVec.size() != secVec.size() ){return false;}
        for(auto code : pdgVec )
        {
            if( firstVec[code] != secVec[code] ){ return false; }
        }
        return true;
    }

    // ZW: fcn for processes in the lheProc struct format
    bool procComp( const lheProc& firstProc, const lheProc& secProc, const std::vector<std::string>& pdgVec )
    {
        for( auto stat : pdgVec )
        {
            if( firstProc.valVecs.at(stat).size() != secProc.valVecs.at(stat).size() ){ return false; }
            if( firstProc.valVecs.at(stat) != secProc.valVecs.at(stat) ){ return false; }
        }
        return true;
    }

    // ZW: fcn for checking whether a list of pdgKtract format
    // processes sourceProcList contains a given process newProc
    bool procVecContains( std::vector<std::shared_ptr<std::map<std::string_view, std::vector<std::string_view>>>>& sourceProcList, 
    std::map<std::string_view, std::vector<std::string_view>>& newProc, const std::vector<std::string>& pdgVec  )
    {
        int noProcs = sourceProcList.size();
        for( auto proc : sourceProcList )
        {
            if( sameProcString( *proc, newProc, pdgVec ) ){ return true; }
        }
        return false;
    }

    // ZW: fcn for checking whether a vector of lheProc structs
    // procList contains a given lheProc nuProc
    bool procListComp( const std::vector<std::shared_ptr<lheProc>>& procList, const lheProc& nuProc, const std::vector<std::string>& pdgVec )
    {
        if( procList.size() != 0 ){
            for(auto proc : procList )
            {
                if( procComp( *proc, nuProc, pdgVec ) ){ return true; }
            }
        }
        return false;
    }

    // ZW: fcn for extracting the different processes
    // in a given PEP format LHE file in the pdgXtract format
    std::vector<std::shared_ptr<std::map<std::string_view, std::vector<std::string_view>>>> procExtractor( const lheNode& lheFile )
    {
        std::vector<std::shared_ptr<std::map<std::string_view, std::vector<std::string_view>>>> procList;
        const static std::vector<std::string> pdgVec = { "-1", "1", "-2", "2", "3", "-9" };
        for( auto event : lheFile.events )
        {
            auto currProc = pgdXtract( *event, pdgVec );
            if( procVecContains( procList, *currProc, pdgVec ) ){ continue; }
            procList.push_back(currProc);
        }
        return procList;
    }

    // ZW: fcn for extracting the differenty processes
    // in a given PEP format LHE file in the lheProc format
    std::vector<std::shared_ptr<lheProc>> processPull( const lheNode& lheFile )
    {
        const static std::vector<std::string> pdgVec = { "-1", "1", "-2", "2", "3", "-9" };
        std::vector<std::shared_ptr<lheProc>> procsList{};
        for( auto event : lheFile.events )
        {
            auto currProc =  std::make_shared<lheProc>( *event );
            if( procListComp( procsList, *currProc, pdgVec ) ){ continue; }
            procsList.push_back( currProc );
        }
        return procsList;
    }

    // ZW: fcn for keeping track of subprocess ordering
    // in LHE file
    int procPos( const std::vector<std::shared_ptr<lheProc>>& evtSet, lheProc& currProc, 
        const std::vector<std::string>& pdgVec )
    {
        int truPos = 0;
        for( auto k = 0 ; k < evtSet.size() ; ++k )
        {
            for( auto stat : pdgVec )
            {
                if( evtSet[k]->valVecs[stat] != currProc.valVecs[stat] ){ break; }   
            }
            return k;
        }
        return evtSet.size();
    }

    // ZW: fcn for extracting the subprocess ordering
    // of LHE file
    std::vector<std::shared_ptr<std::vector<bool>>> procOrder( const lheNode& lheFile, const std::vector<std::shared_ptr<lheProc>>& evtSet )
    {
        const static std::vector<std::string> pdgVec = { "-1", "1", "-2", "2", "3", "-9" };
        std::vector<std::shared_ptr<std::vector<bool>>> eventBools( evtSet.size());
        std::vector<std::vector<bool>> pracBools( evtSet.size(), std::vector<bool> ( lheFile.events.size() ));
        for( auto boolSets : pracBools ){
            std::fill( boolSets.begin(), boolSets.end(), false );
        }
        for( auto k = 0 ; k < lheFile.events.size() ; ++k )
        {
            auto currProc = lheProc(*lheFile.events[k]);
            pracBools[ procPos(evtSet, currProc, pdgVec) ][ k ] = true;
        }
        for( int k = 0 ; k < eventBools.size() ; ++k )
        {
            eventBools[k] = std::make_shared<std::vector<bool>>( pracBools[k] );
        }
        return eventBools; 
    }

    // ZW: fcn for reordering LHE file based on subprocess
    std::shared_ptr<std::vector<std::shared_ptr<event>>> eventReOrder( const lheNode& lheFile, std::vector<bool> relProc )
    {
        auto reOrdered = std::make_shared<std::vector<std::shared_ptr<event>>>();
        reOrdered->reserve( std::count( relProc.begin(), relProc.end(), true ) );
        for( int k = 0 ; k < relProc.size() ; ++k )
        {
            if(!relProc[k]){continue;}
            reOrdered->push_back( lheFile.events[k] );
        }
        return reOrdered;
    }

    // ZW: wrapper for eventReOrder
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> lheReOrder( const lheNode& lheFile )
    {
        auto procSets = processPull( lheFile ); 
        auto relProcs = procOrder( lheFile, procSets ); 
        std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> ordProcs(procSets.size());
        for( int k = 0 ; k < relProcs.size() ; ++k )
        { 
            ordProcs[k] = eventReOrder( lheFile, *relProcs[k] );
        }
        return ordProcs;
    }

    // ZW: transposed event information struct
    struct evtInfo {
    public:
        std::vector<std::string_view> wgts;
        std::vector<std::string_view> scales;
        std::vector<std::string_view> aQEDs;
        std::vector<std::string_view> aQCDs;
        std::vector<std::string_view> nprts;
        std::vector<std::string_view> procIDs;
        evtInfo( const std::vector<std::shared_ptr<PEP::event>>& lheFile = {} ){
            int nEvt = lheFile.size();
            wgts.reserve(nEvt); scales.reserve(nEvt); aQEDs.reserve(nEvt); aQCDs.reserve(nEvt); procIDs.reserve(nEvt);
            for( auto evt : lheFile )
            {
                wgts.push_back(evt->getHead().getWeight());
                scales.push_back(evt->getHead().getScale());
                aQEDs.push_back(evt->getHead().getAQED());
                aQCDs.push_back(evt->getHead().getAQCD());
                nprts.push_back(evt->getHead().getNprt());
                procIDs.push_back(evt->getHead().getProcID());
            }
        }
    };

    // ZW: transposed particle information struct
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
        prtInfo( const std::vector<std::shared_ptr<PEP::event>>& lheFile = {}, const int nPrt = 8 ){
            int nEvt = lheFile.size(); 
            moms.reserve(4*nPrt*nEvt); vtims.reserve(nPrt*nEvt); masses.reserve(nPrt*nEvt); pdgs.reserve(nPrt*nEvt); 
            spins.reserve(nPrt*nEvt); statuses.reserve(nPrt*nEvt); mothers.reserve(2*nPrt*nEvt); icols.reserve(2*nPrt*nEvt);
            for( auto evt : lheFile )
            {
                for( auto prt : evt->getPrts() )
                {
                    moms.push_back( prt->getE() );
                    masses.push_back( prt->getMass() );
                    vtims.push_back( prt->getVTim() );
                    spins.push_back( prt->getSpin() );
                    statuses.push_back( prt->getStatus() );
                    pdgs.push_back( prt->getPDG() );
                    for( int k = 0 ; k < 2 ; ++k )
                    {
                        moms.push_back( prt->getMom()[k] );
                        mothers.push_back( prt->getMothers()[k] );
                        icols.push_back( prt->getColor()[k] );
                    }
                    moms.push_back( prt->getMom()[2] );
                }
            }
        }
    };

    // ZW: transposed LHE file with a single process type
    struct transMonoLHE {
    public:
        evtInfo evtsHead;
        prtInfo evtsData;
        transMonoLHE( const std::vector<std::shared_ptr<PEP::event>>& lheFile = {}, const int nPrt = 8 ){
            evtsHead = evtInfo(lheFile);
            evtsData = prtInfo(lheFile, nPrt);
        }
    };

    // ZW: transposed LHE file ordered by subprocess
    struct transLHE {
    public:
        std::string_view xmlFile;
        std::vector<std::shared_ptr<transMonoLHE>> subProcs;
        transLHE( lheNode& lheFile )
        {
            xmlFile = lheFile.getFile();
            auto procsOrdered = lheReOrder( lheFile ); 
            subProcs = std::vector<std::shared_ptr<transMonoLHE>>( procsOrdered.size() ); 
            for( int k = 0 ; k < procsOrdered.size() ; ++k )
            { 
                subProcs[k] = std::make_shared<transMonoLHE>( *procsOrdered[k], procsOrdered[k]->at(0)->getNprt() );
            }
        }
    };

    // ZW: vector transformation string_to_double
    std::shared_ptr<std::vector<double>> vecStoD( const std::vector<std::string_view> dataVec )
    {
        auto valVec = std::make_shared<std::vector<double>>( dataVec.size() );
        std::transform( dataVec.begin(), dataVec.end(), valVec->begin(), []( const std::string_view& stv ){
            return std::stod(std::string(stv));
        } );
        return valVec;
    }
    
    // ZW: vector transformation string_to_int
    std::shared_ptr<std::vector<int>> vecStoI( const std::vector<std::string_view> dataVec )
    {
        auto valVec = std::make_shared<std::vector<int>>( dataVec.size() );
        std::transform( dataVec.begin(), dataVec.end(), valVec->begin(), []( const std::string_view& stv ){
            return std::stoi(std::string(stv));
        } );
        return valVec;
    }
    
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

    // ZW: bool struct to define which double values
    // to extract transposed from LHE file
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
        std::vector<bool> getBools(){
            return { ebmup, xsecup, xerrup, xmaxup, xwgtup, scalup, aqedup, aqcdup, 
        pup, mass, vtimup, spinup };
        }
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
        std::vector<bool> getBools(){
            return { idbmup, pdfgup, pdfsup, idwtup, nprup, lprup,
            nup, idprup, idup, istup, mothup, icolup };
        }
    };

    // ZW: function for extracting transposed double values
    // from LHE file
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> lheValDoubles( lheNode& lheFile, lheRetDs vals = lheRetDs() )
    {
        // ZW: hard-setting returning g_S instead of a_S for now
        bool aStogS = true;
        auto boolVec = vals.getBools();
        const int noVals = std::count(boolVec.begin(), boolVec.end(), true);
        auto lheAOS = transLHE( lheFile );
        auto lheDos = std::make_shared<std::vector<std::shared_ptr<std::vector<double>>>>(noVals * lheAOS.subProcs.size() ); 
        std::vector<std::shared_ptr<std::vector<double>>> &lheDs = *lheDos;
        int currInd = 0;
        if( boolVec[0] ){ lheDs[currInd] = vecStoD( { lheFile.init->getHead()->ebmup[0], lheFile.init->getHead()->ebmup[1] } ); ++currInd; }
        if( boolVec[1] ){ 
            std::vector<std::string_view> xsecVec( lheFile.init->getLines().size() );
            for( auto line : lheFile.init->getLines() )
            {
                xsecVec.push_back(line->xsecup);
            }
            lheDs[currInd] = vecStoD( xsecVec );
             ++currInd; }
        if( boolVec[2] ){ 
            std::vector<std::string_view> xerrVec( lheFile.init->getLines().size() );
            for( auto line : lheFile.init->getLines() )
            {
                xerrVec.push_back(line->xerrup);
            }
            lheDs[currInd] = vecStoD( xerrVec );
             ++currInd; }
        if( boolVec[3] ){ 
            std::vector<std::string_view> xmaxVec( lheFile.init->getLines().size() );
            for( auto line : lheFile.init->getLines() )
            {
                xmaxVec.push_back(line->xmaxup);
            }
            lheDs[currInd] = vecStoD( xmaxVec );
             ++currInd; } 
        for( int k = 0 ; k < lheAOS.subProcs.size() ; ++k )
        {
            if( boolVec[4] ){ lheDs[currInd] = vecStoD( lheAOS.subProcs[k]->evtsHead.wgts ); ++currInd; }
            if( boolVec[5] ){ lheDs[currInd] = vecStoD( lheAOS.subProcs[k]->evtsHead.scales ); ++currInd; }
            if( boolVec[6] ){ lheDs[currInd] = vecStoD( lheAOS.subProcs[k]->evtsHead.aQEDs ); ++currInd; }
            if( boolVec[7] ){ lheDs[currInd] = vecStoD( lheAOS.subProcs[k]->evtsHead.aQCDs ); 
                if( aStogS ){
                    std::transform( lheDs[currInd]->begin(), lheDs[currInd]->end(), lheDs[currInd]->begin(), 
                    []( double alphaS ){
                        auto gS = std::sqrt( 4. * M_PI * alphaS );
                        return gS;
                    } );
                }
                ++currInd;
            }
            if( boolVec[8] ){ lheDs[currInd] = vecStoD( lheAOS.subProcs[k]->evtsData.moms ); ++currInd; }
            if( boolVec[9] ){ lheDs[currInd] = vecStoD( lheAOS.subProcs[k]->evtsData.masses ); ++currInd; }
            if( boolVec[10] ){ lheDs[currInd] = vecStoD( lheAOS.subProcs[k]->evtsData.vtims ); ++currInd; }
            if( boolVec[11] ){ lheDs[currInd] = vecStoD( lheAOS.subProcs[k]->evtsData.spins ); ++currInd; }
        }

        return lheDos;
    }

    // ZW: function for extracting transposed int values
    // from LHE file
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> lheValInts( lheNode& lheFile, lheRetInts vals = lheRetInts() )
    {
        auto boolVec = vals.getBools();
        const int noVals = std::count(boolVec.begin(), boolVec.end(), true);
        auto lheAOS = transLHE( lheFile );
        auto lheIs = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(noVals * lheAOS.subProcs.size() );
        std::vector<std::shared_ptr<std::vector<int>>> &lheDs = *lheIs;
        int currInd = 0;
        if( boolVec[0] ){ lheDs[currInd] = vecStoI( { lheFile.init->getHead()->idbmup[0], lheFile.init->getHead()->idbmup[1] } ); ++currInd; }
        if( boolVec[1] ){ lheDs[currInd] = vecStoI( { lheFile.init->getHead()->pdfgup[0], lheFile.init->getHead()->pdfgup[1] } ); ++currInd; }
        if( boolVec[2] ){ lheDs[currInd] = vecStoI( { lheFile.init->getHead()->pdfsup[0], lheFile.init->getHead()->pdfsup[1] } ); ++currInd; }
        if( boolVec[3] ){ lheDs[currInd] = vecStoI( { lheFile.init->getHead()->idwtup } ); ++currInd; }
        if( boolVec[4] ){ lheDs[currInd] = vecStoI( { lheFile.init->getHead()->nprup } ); ++currInd; }
        if( boolVec[5] ){ 
            std::vector<std::string_view> lprVec( lheFile.init->getLines().size() );
            for( auto line : lheFile.init->getLines() )
            {
                lprVec.push_back(line->lprup);
            }
            lheDs[currInd] = vecStoI( lprVec );
             ++currInd; }
        for( int k = 0 ; k < lheAOS.subProcs.size() ; ++k )
        {
            if( boolVec[6] ){ lheDs[currInd] = vecStoI( lheAOS.subProcs[k]->evtsHead.nprts ); ++currInd; }
            if( boolVec[7] ){ lheDs[currInd] = vecStoI( lheAOS.subProcs[k]->evtsHead.procIDs ); ++currInd; }
            if( boolVec[8] ){ lheDs[currInd] = vecStoI( lheAOS.subProcs[k]->evtsData.pdgs ); ++currInd; }
            if( boolVec[9] ){ lheDs[currInd] = vecStoI( lheAOS.subProcs[k]->evtsData.statuses ); ++currInd; }
            if( boolVec[10] ){ lheDs[currInd] = vecStoI( lheAOS.subProcs[k]->evtsData.mothers ); ++currInd; }
            if( boolVec[11] ){ lheDs[currInd] = vecStoI( lheAOS.subProcs[k]->evtsData.icols ); ++currInd; }
        }
        return lheIs;
    }
}