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

#ifndef _REX_CC_
#define _REX_CC_

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
#include "REX.h"
#include <typeinfo>

// ZW: all fcns within the REX standard sit in the
// namespace REX
// Note that as a convention, std::string_view objects will be
// referred to as strings unless the difference is relevant
namespace REX
{

    using sortFcn = std::function<std::shared_ptr<std::vector<size_t>>(std::vector<std::string_view>)>;
    using statSort = std::function<std::shared_ptr<std::vector<size_t>>(std::string_view, std::vector<std::string_view>)>;

    // ZW: index sorting function, which returs vector
    // of the indices of the original vector sorted 
    // by default in ascending order
    // ie, for [5.0, 0.25, 2.0, 9.2] returns [1, 2, 0, 3]
    template <typename T>
    std::shared_ptr<std::vector<size_t>> indSort(const std::vector<T> &vector, std::function<bool(const T&, const T&)> comp = std::less<T>())
    {
        auto sorted = std::make_shared<std::vector<size_t>>(vector.size());
        std::iota(sorted->begin(), sorted->end(), 0);
        std::stable_sort(sorted->begin(), sorted->end(), [&](size_t i, size_t j) { return comp(vector[i], vector[j]); });
        return sorted;
    }

    // ZW: wrapper for indSort for comparing string-type arguments representing integers
    template <typename T>
    std::shared_ptr<std::vector<size_t>> stoiSort(const std::vector<T> &vector)
    {
        std::function<bool(const T&, const T&)> stoicomp = [](const T& i, const T& j) { 
        return std::stoi(std::string(i)) < std::stoi(std::string(j)); };
        return indSort(vector, stoicomp);
    }
    template std::shared_ptr<std::vector<size_t>> stoiSort<std::string_view>(const std::vector<std::string_view> &vector);

    // ZW: wrapper for indSort for comparing string-type arguments representing doubles
    template <typename T>
    std::shared_ptr<std::vector<size_t>> stodSort(const std::vector<T> &vector)
    {
        std::function<bool(const T&, const T&)> stodcomp = [](const T& i, const T& j) { return std::stod(std::string(i)) < std::stod(std::string(j)); };
        return indSort(vector, stodcomp);
    }

    // ZW: templated fcn for finding the order of elements in a vector to_sort
    // based on their order in a reference vector reference
    // Elements not found in reference are represented by npos,
    // including if to_sort is longer than reference
    template <typename T>
    std::shared_ptr<std::vector<size_t>> getRefOrder(const std::vector<T>& reference, const std::vector<T>& to_sort) {
        std::unordered_map<T, std::queue<size_t>> indexMap;

        // Populate indexMap with indices from vec1
        for (size_t i = 0; i < reference.size(); ++i) {
            indexMap[reference[i]].push(i);
        }

        auto order = std::make_shared<std::vector<size_t>>(std::vector<size_t>(to_sort.size(), npos));
        //order->reserve(to_sort.size()); // Pre-allocate memory
        size_t pos = 0;
        for (const auto& elem : to_sort) {
            auto it = indexMap.find(elem);
            if (it != indexMap.end() && !it->second.empty()) {
                order->at(pos) = (it->second.front());
                it->second.pop();
            } //else {
                // Element in vec2 not found in vec1
            //    order->at(pos) = npos;
            //}
            ++pos;
        }

        return order;
    }
    template std::shared_ptr<std::vector<size_t>> getRefOrder<std::string_view>(const std::vector<std::string_view>& reference, const std::vector<std::string_view>& to_sort);

    // ZW: minimal fcn for counting the amount of times
    // a given search term appears in a string
    int nuStrCount( std::string_view searchString, std::string_view searchTerm )
    {
        int count = 0;
        size_t pos = 0;
        while((pos = searchString.find(searchTerm, pos)) != npos ){
            ++count;
            ++pos;
        }
        return count;
    }

    // ZW: fcn for finding the location of each
    // entry of seachTerm in the given string textFile
    // Pre-allocates vector memory using nuStrCount
    std::shared_ptr<std::vector<size_t>> nuFindEach( std::string_view textFile, std::string_view searchTerm )
    {
        auto eachPos = std::make_shared<std::vector<size_t>>();
        eachPos->reserve( nuStrCount(textFile, searchTerm) );
        eachPos->push_back( textFile.find( searchTerm ) );
        size_t currPos = textFile.find( searchTerm, eachPos->at(0) + 1 );
        while( currPos != npos )
        {
            eachPos->push_back( currPos );
            currPos = textFile.find( searchTerm, currPos + 1 );
        }
        return eachPos;
    }

    // ZW: fcn for splitting a string into a vector of strings,
    // each element differentiated by linebreaks in the original string
    // Removes sequential linebreaks, ie "\n\n\n" would
    // only result in a single element separation
    std::shared_ptr<std::vector<std::string_view>> nuLineSplitter( std::string_view currEvt )
    {
        auto lineBreaks = nuFindEach( currEvt, "\n" );
        std::vector<size_t> trueBreaks;
        trueBreaks.reserve( lineBreaks->size() );
        for( size_t k = 0 ; k < lineBreaks->size() - 1 ; ++k )
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
        if( currEvt.substr( startPos ).size() > 1 ){ splitLines->push_back( currEvt.substr( startPos ) ); }
        return splitLines;
    }

    // ZW: fcn for finding each linebreak in a string,
    // returning a vector of the positions of "\n" characters
    // Ignores sequential linebreaks, ie would only return { }
    // for the string "\n\n\n\n"
    std::shared_ptr<std::vector<size_t>> lineFinder( std::string_view currEvt, size_t startPos = 0, size_t endPos = npos )
    {
        auto lineBreaks = nuFindEach( currEvt.substr( startPos, endPos - startPos), "\n" );
        auto truBreaks = std::make_shared<std::vector<size_t>>();
        truBreaks->reserve( lineBreaks->size() );
        for( size_t k = 0 ; k < lineBreaks->size() ; ++k )
        {
            if( int( (*lineBreaks)[k+1] - (*lineBreaks)[k]) == 1){continue;}
            truBreaks->push_back( (*lineBreaks)[k] );
        }
        return truBreaks;
    }
    
    // ZW: fcn for splitting a string into a vector of strings,
    // each element separated by blankspace (" ") in the original string
    // Ignores sequential blankspaces, as well as linebreaks
    // ie "hello     \n\n\n     world" would return {"hello", "world"}
    // Does not ignore linebreaks that are not separated from words
    // by anything other than blankspace,
    // ie "hello     \n\n\nworld   \n\n" would return {"hello", "\n\nworld"}
    std::shared_ptr<std::vector<std::string_view>> nuWordSplitter( std::string_view currEvt )
    {
        std::vector<size_t> noSpace;
        size_t nuStart = currEvt.find_first_not_of( " " );
        size_t nuEnd = currEvt.find(" ", nuStart+1 );
        auto splitWords = std::make_shared<std::vector<std::string_view>>();
        splitWords->reserve(13);
        while( nuStart != npos )
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
    // elements separated by any form of blankspace in the original string
    // Ignores sequential blankspaces of all forms
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
    bool clStringComp( std::string_view org, std::string comp ){
        return std::equal( org.begin(), org.end(), comp.begin(), comp.end(), 
        []( const char& x, char y ){ return (std::toupper(x) == std::toupper(y)); } );
    }
    bool clStringComp( std::string_view org, std::string_view comp ){
        return std::equal( org.begin(), org.end(), comp.begin(), comp.end(), 
        []( const char& x, char y ){ return (std::toupper(x) == std::toupper(y)); } );
    }
    bool clStringComp( std::string org, std::string_view comp ){
        return std::equal( org.begin(), org.end(), comp.begin(), comp.end(), 
        []( const char& x, char y ){ return (std::toupper(x) == std::toupper(y)); } );
    }
    bool clStringComp( std::string org, std::string comp ){
        return std::equal( org.begin(), org.end(), comp.begin(), comp.end(), 
        []( const char& x, char y ){ return (std::toupper(x) == std::toupper(y)); } );
    }
    // template<typename Str1, typename Str2>
    // bool clStringComp( const Str1& org, const Str2& comp ){
    //     return std::equal( org.begin(), org.end(), comp.begin(), comp.end(), 
    //     []( const char& x, char y ){ return (std::toupper(x) == std::toupper(y)); } );
    // }
    // template<typename Str1Pt, typename Str2>
    // bool clStringComp( const Str1Pt& orgStrt, const Str1Pt& orgEnd, const Str2& comp ){
    //     return std::equal( orgStrt, orgEnd, comp.begin(), comp.end(), 
    //     []( const char& x, char y ){ return (std::toupper(x) == std::toupper(y)); } );
    // }

    // ZW: templated fcn for finding a caseless substring searchTerm in srcFile
    // On failure to find searchTerm, returns REX::npos
    template<typename Str1, typename Str2>
    size_t clStringFind( const Str1& srcFile, const Str2& searchTerm, size_t strtPt = 0 ){
        size_t strLen = searchTerm.size();
        if( srcFile.size() == 0  || srcFile.size() < strLen ){ return npos; }
        for( size_t k = strtPt ; k < srcFile.size()  - strLen; ++k )
        {
            if( clStringComp( srcFile.substr(k, strLen), searchTerm ) ){ return k; }
        }
        return npos;
    }

    // ZW: templated fcn for finding a caseless substring searchTerm of srcFile 
    // fulfilling a particular predicate cond( size_t, string )
    template<typename Str1, typename Str2>
    size_t clStringFindIf( const Str1& srcFile, const Str2& searchTerm, std::function<bool(size_t&, const Str1&)>& cond, size_t strtPt = 0 )
    {
        auto currPt = clStringFind( srcFile, searchTerm, strtPt ); 
        bool condStat = cond( currPt, srcFile );
        while( !( condStat ) && currPt != npos)
        {
            currPt = clStringFind( srcFile, searchTerm, currPt + 1 );
            condStat = cond( currPt, srcFile );
        } 
        return currPt;
    }

    // ZW: templated fcn for counting the number of occurances of 
    // caseless substring searchTerm in string-like object srcFile
    template<typename Str1, typename Str2>
    int clStrCount( Str1 srcFile, Str2 searchTerm )
    {
        int count = 0;
        size_t pos = 0;
        while((pos = clStringFind( srcFile, searchTerm, pos ) ) != npos ){
            ++count;
            ++pos;
        }
        return count;
    }

    // ZW: templated fcn for finding each instance of 
    // of substring searchTerm of string-like object srcFile
    template<typename Str1, typename Str2>
    std::shared_ptr<std::vector<size_t>> clFindEach( Str1 srcFile, Str2 searchTerm )
    {
        auto eachPos = std::make_shared<std::vector<size_t>>();
        auto nos = clStrCount(srcFile, searchTerm);
        if( nos == 0 ){ return eachPos; }
        eachPos->reserve( nos );
        eachPos->push_back( clStringFind( srcFile, searchTerm ) );
        size_t currPos = clStringFind( srcFile, searchTerm, eachPos->at(0) + 1);
        while( currPos != npos )
        {
            eachPos->push_back( currPos );
            currPos = clStringFind( srcFile, searchTerm, currPos + 1 );
        }
        return eachPos;
    }

    // ZW: fcn for finding left angle bracket
    // indicating the start of a new node in an XML file
    size_t nodeStartFind( std::string_view parseFile, size_t strtPos )
    {
        auto retPtr = parseFile.find("<", strtPos);
        while( parseFile[retPtr + 1] == '!' || parseFile[retPtr +1] == '/' || parseFile[retPtr +1] == '?' ){
            retPtr = parseFile.find("<", retPtr +1);
        }
        return retPtr;
    }

    size_t endNodeStartFind( std::string_view parseFile, size_t strtPos )
    {
        return parseFile.find(">", nodeStartFind( parseFile, strtPos ));
    }

    std::pair<size_t, size_t> startNodePts( std::string_view parseFile, size_t strtPos )
    {
        return { nodeStartFind( parseFile, strtPos ), endNodeStartFind( parseFile, strtPos ) };
    }

    // ZW: fcn for finding left angle bracket
    // indicating an end of a node in an XML file
    size_t nodeEndFind( std::string_view parseFile, size_t strtPos )
    { 
        auto retPtr = parseFile.find("<", strtPos); 
        while( parseFile[retPtr + 1] != '/' ){ 
            retPtr = parseFile.find("<", retPtr +1);
        } 
        return retPtr;
    }

    size_t endNodeEndFind( std::string_view parseFile, size_t strtPos )
    {
        return parseFile.find(">", nodeEndFind( parseFile, strtPos ));
    }

    std::pair<size_t, size_t> endNodePts( std::string_view parseFile, size_t strtPos )
    {
        return { nodeEndFind( parseFile, strtPos ), endNodeEndFind( parseFile, strtPos ) };
    }

    // ZW: struct for handling tags in XML node opening tags
        void xmlTag::setVal( std::string_view valSet ){ modded = true; val = valSet; }
        void xmlTag::setId( std::string_view idSet ){ modded = true; id = idSet; }
        std::string_view xmlTag::getVal(){ return val; }
        std::string_view xmlTag::getId(){ return id; }
        bool xmlTag::isModded(){ return modded; }
        xmlTag::xmlTag(){ modded = false; return; }
        xmlTag::xmlTag( xmlTag& oldTag ){
            modded = false; val = oldTag.getVal(); id = oldTag.getId();
        }
        xmlTag::xmlTag( std::string_view initId, std::string_view initVal){
            modded = false; val = initVal; id = initId;
        }

    // ZW: function for parsing XML opening
    // tags and returning the next header tag
    std::shared_ptr<xmlTag> xmlTagParser( std::string_view tagLine, size_t& equPt )
    {
        auto tagBreaker = tagLine.find_first_not_of(" ", equPt+1); // ZW: need to determine what type of quotation marks are used
        auto tagEnder = tagLine.find( tagLine[tagBreaker], tagBreaker+1);
        auto attrEnd = tagLine.find_last_not_of(" ", equPt - 1) ;
        auto attrStart = tagLine.find_last_of(" ", attrEnd) + 1;
        auto tagPtr = std::make_shared<xmlTag>(tagLine.substr(attrStart, attrEnd - attrStart + 1), tagLine.substr(tagBreaker + 1, tagEnder - tagBreaker - 1));
        equPt = tagLine.find("=", equPt + 1); // ZW: modifies input equPt to point to the next equality sign in tagLine
        return tagPtr;
    }

    // ZW: struct for handling the tree structure of XML files,
    // essentially just giving the positions of the beginning and
    // end of each node s.t. the proper node structures can accurately
    // detail where children begin and end while allowing for personal
    // content between child nodes 
        xmlTree::xmlTree(){ return; }
        xmlTree::xmlTree( std::string_view file ){
            origin = file;
            children = std::make_shared<std::vector<std::shared_ptr<xmlTree>>>();
            start = file.find_first_not_of(" \n\r\f\t\v");
            if( file.compare(start, 1, "<") != 0 ) { 
                faux = true; 
                contSt = start; 
                end = std::min( nodeStartFind(file, start), nodeEndFind(file, start) ); 
                contEnd = end; 
                initialised = true; 
                return; 
            }
            if( file.compare(start + 1, 1, "!") == 0 || file.compare(start + 1, 1, "?")  == 0 ) {
                faux = true;
                contSt = start;
                contEnd = file.find(">", start + 1);
                end = std::min( nodeStartFind(file, contEnd), nodeEndFind(file, contEnd) );
                initialised = true;
                return;
            }
            auto stEnd = file.find(">", start);
            if( file.compare(stEnd - 1, 1, "/" ) == 0 ) {
                end = file.find_first_not_of(" \n\r\f\t\v", stEnd + 1); 
                contSt = npos;
                contEnd = npos;
                initialised = true;
                return; 
            }
            contSt = stEnd + 1;
            auto stPos = nodeStartFind(file, start + 1);
            stEnd = nodeEndFind(file, start + 1);
            contEnd = std::min(stPos, stEnd);
            while( stPos < stEnd )
            {
                children->push_back( std::make_shared<xmlTree>( file, stPos, stEnd ) ); 
            }
            stEnd = endNodeEndFind(file, stEnd); 
            end = file.find_first_not_of(" \n\r\f\t\v", stEnd + 1); 
            initialised = true;
        }
        xmlTree::xmlTree( std::string_view file, size_t& strt, size_t& nd ){ 
            origin = file;
            children = std::make_shared<std::vector<std::shared_ptr<xmlTree>>>();
            start = file.find_first_not_of(" \n\r\f\t\v", strt);
            if( file.compare(start, 1, "<") != 0)  { 
                faux = true;
                contSt = start;
                strt = nodeStartFind(file, start);
                nd = nodeEndFind(file, start);
                end = std::min( strt, nd ); 
                contEnd = end;  
                initialised = true; 
                return; 
            }
            if( file.compare(start + 1, 1, "!") == 0 ) {
                faux = true;
                contSt = start;
                contEnd = file.find(">", start + 1);
                strt = nodeStartFind(file, contEnd);
                nd = nodeEndFind(file, contEnd);
                end = std::min( strt, nd );
                initialised = true;
                return;
            }
            auto stEnd = file.find(">", start); 
            if( file.compare(stEnd - 1, 1, "/" ) == 0 ) {  
                end = file.find_first_not_of(" \n\r\f\t\v", stEnd + 1); 
                contSt = npos;
                contEnd = npos; 
                strt = nodeStartFind(file, start);
                nd = nodeEndFind(file, start); 
                initialised = true; 
                return;
            }
            contSt = stEnd + 1; 
            strt = nodeStartFind(file, start + 1); 
            nd = nodeEndFind(file, start + 1); 
            contEnd = std::min(strt, nd); 
            while( strt < nd )
            { 
                children->push_back( std::make_shared<xmlTree>( file, strt, nd ) ); 
            }
            end = file.find_first_not_of(" \n\r\f\t\v", endNodeEndFind(file, nd) + 1);
            initialised = true; 
            strt = end; 
            nd = nodeEndFind(file, strt);
        }

    // ZW: struct for handling nodes in generic XML files
        xmlNode::xmlNode(){ modded = false; return; }
        xmlNode::xmlNode( const std::string_view originFile, const size_t& begin, const std::vector<std::shared_ptr<xmlNode>>& childs ){
            modded = false;
            xmlFile = originFile.substr( begin );
            structure = xmlTree( originFile );
            faux = structure.isFaux();
            start = structure.getStart();
            end = structure.getEnd();
            size_t trueStart = xmlFile.find_first_not_of("< \n\r\f\t\v", start+1);
            name = xmlFile.substr( trueStart, xmlFile.find_first_of(">/ \n\r\f\t\v", trueStart) - trueStart );
            content = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            for( auto child : childs ){
                children.push_back( child );
            }
        }
        xmlNode::xmlNode( xmlTree &tree ){ 
            modded = false; 
            structure = tree; 
            if( !structure.isInit() ){ return; } 
            xmlFile = structure.getOrigin(); 
            faux = structure.isFaux(); 
            start = structure.getStart(); 
            end = structure.getEnd(); 
            size_t trueStart = xmlFile.find_first_not_of("< \n\r\f\t\v", start); 
            name = xmlFile.substr( trueStart, xmlFile.find_first_of(">/ \n\r\f\t\v", trueStart) - trueStart );
            content = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() ); 
            for( auto& child : *(structure.getChildren()) ){ 
                children.push_back( std::make_shared<xmlNode>( *child ) ); 
            } 
        }
        std::vector<std::shared_ptr<xmlNode>> xmlNode::getChildren(){ return children; }
        std::vector<std::shared_ptr<xmlTag>> xmlNode::getTags(){ return tags; }
        std::string_view xmlNode::getFile(){ return xmlFile; }
        std::string_view xmlNode::getName(){ return name; }
        std::string_view xmlNode::getContent(){ return content; }
        size_t xmlNode::getStart(){ return start; }
        size_t xmlNode::getEnd(){ return end; }
        xmlTree xmlNode::getTree(){ return structure; }
        bool xmlNode::isModded(){ return modded; }
        bool xmlNode::isModded( bool deep ){
            bool modStat = isModded();
            if( !deep ){ return modStat; }
            for( auto child : children ){ modStat = (modStat || child->isModded( deep )); }
            return modStat;
        }
        bool xmlNode::isWritten(){ return written; }
        bool xmlNode::isParsed(){ return parsed; }
        bool xmlNode::isFaux(){ return faux; }
        bool xmlNode::hasChildren(){ return children.size() > 0; }
        void xmlNode::setModded( bool mod ){ modded = mod; }
        bool xmlNode::deepModded(){ return deepMod; }
        bool xmlNode::deepParse(){ return deepParsed; }
        void xmlNode::parser( bool recursive ){
            parsed = parse( recursive );
        }
        void xmlNode::addChild( std::shared_ptr<xmlNode> child ){ modded = true; children.push_back(child); }
        void xmlNode::addTag( std::shared_ptr<xmlTag> tag ){ modded = true; tags.push_back(tag); }
        void xmlNode::setFile( std::string_view file ){ modded = true; xmlFile = file; }
        void xmlNode::setName( std::string_view newName ){ modded = true; name = newName; }
        void xmlNode::setCont( std::string_view cont ){ modded = true; content = cont; }

        bool xmlNode::parse(){ 
            auto topStat = parseTop();
            auto contStat = parseContent();
            return ( topStat && contStat );
        }
        bool xmlNode::parse( bool recurs )
        {
            bool parseSt = parse();
            if( !recurs ){ return parseSt; }
            bool childSt = parseChildren( recurs );
            deepMod = true;
            return (parseSt && childSt );
        }
        bool xmlNode::parseTop(){
            if( xmlFile == "" ){ return false; }
            if( isFaux() ){ return true; }
            size_t eqSgn = xmlFile.find( "=", start ); size_t nodeInitEnd = xmlFile.find( ">", start );
            while( eqSgn < nodeInitEnd ){ tags.push_back( xmlTagParser( xmlFile, eqSgn ) ); }
            return true;
        }
        bool xmlNode::parseContent(){
            if( xmlFile == "" ){ return false; }
            end = structure.getContEnd();
            for( auto branch : *(structure.getChildren()) ){
                children.push_back( std::make_shared<xmlNode>( *branch ) ); 
            }
            return true;
        }
        bool xmlNode::parseChildren( bool recursive ){
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
        void xmlNode::headWriter() {
            if( isFaux() ){ return; }
            nodeHeader =  "<" + std::string(name) ;
            for( auto tag : tags ){
                nodeHeader += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            nodeHeader += ">";
        }
        void xmlNode::endWriter() {
            if( isFaux() ){ return; }
            auto endSt = xmlFile.find_last_of("<", end);
            nodeEnd = xmlFile.substr( endSt, end - endSt );
        }
        void xmlNode::contWriter() {
            if( hasChildren() ){
            nodeContent = std::string(content.substr(0, children[0]->start - 1 ));
            } else {
            nodeContent = std::string(content);
            }
        }
        void xmlNode::childWriter() {
            for(auto child : children){
                nodeContent += (*child->nodeWriter());
            }
        }
        void xmlNode::endFinder(){
            auto headEnd = xmlFile.find(">", start);
            auto slashPos = xmlFile.find("/", start);
            if( headEnd > slashPos ){ end = headEnd; }
            else{ end = xmlFile.find( ">", xmlFile.find( "</" + std::string(name), start )); }
            if( end == npos ){ end = xmlFile.size(); return; }
            end += 2;
        }
        void xmlNode::fullWriter(){
            if( isModded() ){
            headWriter();
            contWriter();
            childWriter();
            endWriter();
            writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
            written = true;
            modded = false;
            } else if( !isWritten() ){
            writtenSelf = std::make_shared<std::string>( xmlFile.substr( start, end - start ) );
            written = true;
            }
        }
        
        void xmlNode::childCounter( int& noChilds )
        {
            for( auto child : children )
            {
                child->childCounter( noChilds );
                if( child->end == 0 || child->isFaux() ){ --noChilds; }
            }
            noChilds += children.size();
        }  
        int xmlNode::childCounter() {
            int noChilds = 0;
            childCounter( noChilds );
            return noChilds;
        }  
        std::shared_ptr<std::string> xmlNode::nodeWriter() {
            if( isModded( true ) || !isWritten() ){ fullWriter(); }
            return writtenSelf;
        }


    // ZW: function for large scale parsing of XML files
    // sequentially goes through the document and
    // recursively calls itself while the next node
    // beginning is closer than the next node ending
    std::shared_ptr<xmlNode> xmlPtrParser( std::string_view parseFile, size_t& initPos, size_t& endPos )
    {
        auto currNode = std::make_shared<xmlNode>(parseFile, initPos);
        size_t equalSign = parseFile.find("=", initPos);
        size_t nodeInitEnd = parseFile.find(">", initPos);
        initPos = nodeStartFind( parseFile, initPos + 1 );
        while( equalSign < nodeInitEnd ){
            currNode->addTag( xmlTagParser(parseFile, equalSign) );
        }
        while( initPos < endPos )
        {
            currNode->addChild(xmlPtrParser( parseFile, initPos, endPos ));
        }
        
        initPos = nodeStartFind( parseFile, endPos );
        endPos = nodeEndFind( parseFile, endPos + 1 );
        return currNode;
    }

    // ZW: struct for handling rwgt parameter sets
    // in the LHE header initrwgt node
        int headWeight::headWeight::getId(){ return id; }
        std::string_view headWeight::getTag(){ return idTag; }
        bool headWeight::hasTag(){ return (idTag.size() > 0); }
        headWeight::headWeight(){ name = "weight"; return; }
        headWeight::headWeight( std::string_view paramSet, const size_t& begin ) : xmlNode(){ name = "weight"; xmlFile = paramSet; content = paramSet; return; }
        headWeight::headWeight( std::string_view paramSet, std::string_view idText, int idNo, const size_t& begin ) : xmlNode(){
            name = "weight"; xmlFile = paramSet; content = paramSet; idTag = idText; id = idNo;
        }
        headWeight::headWeight( xmlNode& node ) : xmlNode( node ){
            parser( false );
            name = "weight";
            for (auto tag : tags ){
                if( tag->getId() == "id" ){
                    idTag = tag->getVal().substr(0, tag->getVal().find_last_of("_") - 1 );
                    id = std::stoi( std::string( tag->getVal().substr( idTag.size() + 1 ) ) );
                }
            }
        }
        headWeight::headWeight( xmlNode* node ) : xmlNode( *node ){
            parser( false );
            name = "weight";
            for (auto tag : tags ){
                if( tag->getId() == "id" ){
                    idTag = tag->getVal().substr(0, tag->getVal().find_last_of("_") - 1 );
                    id = std::stoi( std::string( tag->getVal().substr( idTag.size() + 1 ) ) );
                }
            }
        }
        headWeight::headWeight( std::shared_ptr<xmlNode> node ) : xmlNode( *node ){
            parser( false );
            name = "weight";
            for (auto tag : tags ){
                if( tag->getId() == "id" ){
                    idTag = tag->getVal().substr(0, tag->getVal().find_last_of("_") - 1 );
                    id = std::stoi( std::string( tag->getVal().substr( idTag.size() + 1 ) ) );
                }
            }
        }
        headWeight::headWeight( xmlTree& tree ) : xmlNode( tree ){
            parser( false );
            name = "weight";
            for (auto tag : tags ){
                if( tag->getId() == "id" ){
                    idTag = tag->getVal().substr(0, tag->getVal().find_last_of("_") - 1 );
                    id = std::stoi( std::string( tag->getVal().substr( idTag.size() + 1 ) ) );
                }
            }
        }
        headWeight::headWeight( xmlTree* tree ) : xmlNode( *tree ){
            parser( false );
            name = "weight";
            for (auto tag : tags ){
                if( tag->getId() == "id" ){
                    idTag = tag->getVal().substr(0, tag->getVal().find_last_of("_") - 1 );
                    id = std::stoi( std::string( tag->getVal().substr( idTag.size() + 1 ) ) );
                }
            }
        }
        headWeight::headWeight( std::shared_ptr<xmlTree> tree ) : xmlNode( *tree ){
            parser( false );
            name = "weight";
            for (auto tag : tags ){
                if( tag->getId() == "id" ){
                    idTag = tag->getVal().substr(0, tag->getVal().find_last_of("_") - 1 );
                    id = std::stoi( std::string( tag->getVal().substr( idTag.size() + 1 ) ) );
                }
            }
        }
        headWeight::headWeight( std::string_view paramSet, std::string& idText, unsigned int idNo, const size_t& begin ) : xmlNode(){
            name = "weight"; xmlFile = paramSet; content = paramSet; idTag = idText; id = idNo;
        }
        headWeight::headWeight( std::string_view paramSet, std::string& idText){
            name = "weight"; xmlFile = paramSet; content = paramSet; idTag = idText;
        }
        void headWeight::setId( std::string identity ){ modded = true; idTag = identity; }
        void headWeight::headWriter(){
            if( tags.size() == 0 ){
                if( idTag == "" ){ nodeHeader = "<weight>"; return; }
                if( id == npos ){ nodeHeader = "<weight id=\"" + std::string(idTag) + "\">"; return; }
                nodeHeader = "<weight id=\"" + std::string(idTag) + std::to_string(id) + "\">";
                return;
            }
            nodeHeader = "<weight";
            for( auto tag : tags ){
                nodeHeader += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            nodeHeader += ">";
        }
        void headWeight::headWriter( bool incId ){
            if( !incId ){ headWriter(); return; }
            if( idTag == "" ){ headWriter(); return; }
            if( id == npos ){ nodeHeader = "<weight id=\"" + std::string( idTag ) + "\""; }
            else{ nodeHeader = "<weight id=\"" + std::string( idTag ) + "_" + std::to_string(id) + "\""; }
            for( auto tag : tags ){
                if( tag->getId() == "id" ){ continue; }
                nodeHeader += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            nodeHeader += ">";
        }
        void headWeight::endWriter() {
            nodeEnd = "</weight>\n";
        }
        void headWeight::contWriter() { 
            nodeContent = std::string( content );
        }
        void headWeight::childWriter() {
            for( auto child : children){
                if( child->getName() == "weight" ){ continue; }
                nodeContent += *(child->nodeWriter());
            }
        }
        void headWeight::childWriter( bool hasChildren ){
            if( hasChildren ){ childWriter(); }
        }
        void headWeight::fullWriter(){
            if( isModded() || !isWritten() ){
                headWriter();
                contWriter();
                childWriter();
                endWriter();
                writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
                written = true;
                modded = false;
            }
        }
        void headWeight::fullWriter( bool incId, bool hasChildren ){
            if( isModded() || !isWritten() ){
            headWriter( incId );
            contWriter();
            childWriter( hasChildren );
            endWriter();
            writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
            modded = false;
            written = true;
            }
        }

    // ZW: struct for handling rwgt groups
    // in the LHE header initrwgt node
        bool weightGroup::getIncId(){ return includeId; }
        void weightGroup::setIncId( bool nuIncId ){ includeId = nuIncId; }
        std::vector<std::shared_ptr<headWeight>> weightGroup::getWgts(){ return paramSets; }
        void weightGroup::addWgt( headWeight nuWgt ){ modded = true; paramSets.push_back( std::make_shared<headWeight>( nuWgt ) ); if( nuWgt.hasTag() ){ includeId = true; } }
        void weightGroup::addWgt( std::shared_ptr<headWeight> nuWgt ){ modded = true; paramSets.push_back( nuWgt); if( nuWgt->hasTag() ){ includeId = true; }}
        weightGroup::weightGroup() : xmlNode(){ name = "weightgroup"; return; }
        weightGroup::weightGroup( std::vector<std::shared_ptr<headWeight>> nuWgts ) : xmlNode(){ name = "weightgroup"; paramSets = nuWgts; for( auto wgt : nuWgts ){ if( wgt->hasTag() ){ includeId = true; } } }
        weightGroup::weightGroup( std::vector<std::string> nuWgts ) : xmlNode(){
            name = "weightgroup";
            for( auto wgt : nuWgts ){
                paramSets.push_back( std::make_shared<headWeight>( wgt ) );
            }
            for( auto wgt : paramSets ){ if( wgt->hasTag() ){ includeId = true; } }
        }
        weightGroup::weightGroup( xmlNode& wgtNode ) : xmlNode( wgtNode ){
            parser( true );
            name = "weightgroup";
            paramSets.reserve( children.size() );
            for(  auto child : children ){
                if( child->getName() == "weight" ){ paramSets.push_back( std::make_shared<headWeight>( *child ) ); }
            }
            for( auto wgt : paramSets ){ if( wgt->hasTag() ){ includeId = true; } }
        }
        weightGroup::weightGroup( xmlNode* wgtNode ) : xmlNode( *wgtNode ){
            parser( true );
            name = "weightgroup";
            paramSets.reserve( children.size() );
            for(  auto child : children ){
                if( child->getName() == "weight" ){ paramSets.push_back( std::make_shared<headWeight>( *child ) ); }
            }
            for( auto wgt : paramSets ){ if( wgt->hasTag() ){ includeId = true; } }
        }
        weightGroup::weightGroup( xmlTree& wgtTree ) : xmlNode( wgtTree ){ 
            parser( true );
            name = "weightgroup";
            paramSets.reserve( children.size() );
            for(  auto child : children ){
                if( child->getName() == "weight" ){ paramSets.push_back( std::make_shared<headWeight>( *child ) ); }
            }
            for( auto wgt : paramSets ){ if( wgt->hasTag() ){ includeId = true; } }
        }
        weightGroup::weightGroup( xmlTree* wgtTree ) : xmlNode( *wgtTree ){
            parser( true );
            name = "weightgroup";
            paramSets.reserve( children.size() );
            for(  auto child : children ){
                if( child->getName() == "weight" ){ paramSets.push_back( std::make_shared<headWeight>( *child ) ); }
            }
            for( auto wgt : paramSets ){ if( wgt->hasTag() ){ includeId = true; } }
        }
        weightGroup::weightGroup( std::shared_ptr<xmlTree> wgtTree ) : xmlNode( *wgtTree ){
            parser( true );
            name = "weightgroup";
            paramSets.reserve( children.size() );
            for(  auto child : children ){
                if( child->getName() == "weight" ){ paramSets.push_back( std::make_shared<headWeight>( *child ) ); }
            }
            for( auto wgt : paramSets ){ if( wgt->hasTag() ){ includeId = true; } }
        }
        weightGroup::weightGroup( const std::string_view originFile, const size_t& begin, const std::vector<std::shared_ptr<xmlNode>>& childs )
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
        void weightGroup::headWriter() {
            nodeHeader = "<weightgroup";
            if( rwgtName != "" ){ nodeHeader += " name=\"" + std::string( rwgtName ) +"\""; }else if( isModded() ){ nodeHeader += " name=\"rex_reweighting\""; }
            if( wgtNamStrat != "" ){ nodeHeader += " weight_name_strategy=\"" + std::string( wgtNamStrat ) +"\""; }
            else if( wgtNamStrat == "" && includeId ){ nodeHeader += " weight_name_strategy=\"includeIdInWeightName\"";}
            nodeHeader += ">";
        }
        void weightGroup::contWriter() {
            nodeContent = "\n";
            for( auto wgt : paramSets ){
                nodeContent += (*wgt->nodeWriter());
            }
        }
        void weightGroup::childWriter() {
            for(auto child : children){
                if( child->getName() == "weight" ){ continue; }
                nodeContent += (*child->nodeWriter());
            }
        }
        void weightGroup::childWriter( bool hasChildren ){
            if( hasChildren ){ childWriter(); }
            return;
        }
        void weightGroup::endWriter() { nodeEnd = "</weightgroup>\n"; }

        std::vector<std::shared_ptr<weightGroup>> initRwgt::getGroups(){ return groups; }
        size_t initRwgt::noGrps(){ return groups.size(); }
        void initRwgt::addGroup( weightGroup nuGroup ){ 
            modded = true;
            auto nuGrpPtr = std::make_shared<weightGroup>( nuGroup );
            if( grpInit( nuGrpPtr ) ){ groups.push_back( std::make_shared<weightGroup>( nuGroup ) );  }
        }
        void initRwgt::addGroup( std::shared_ptr<weightGroup> nuGroup ){ 
            modded = true;
            if( grpInit( nuGroup ) ){ groups.push_back( nuGroup ); }
        }
        void initRwgt::addWgt( unsigned int index, std::shared_ptr<headWeight> nuWgt ){
            if( index < groups.size() ){ modded = true; groups[index]->addWgt( nuWgt ); }
            else throw std::range_error( "Appending weight to uninitialised weightgroup." );
        }
        void initRwgt::addWgt( unsigned int index, headWeight nuWgt ){
            if( index < groups.size() ){ modded = true; groups[index]->addWgt( nuWgt ); }
            else throw std::range_error( "Appending weight to uninitialised weightgroup." );
        }
        initRwgt::initRwgt() : xmlNode(){ name = "initrwgt"; return; }
        initRwgt::initRwgt( std::vector<std::shared_ptr<xmlNode>> nuGroups ) : xmlNode(){
            name = "initrwgt";
            for( auto group : nuGroups ){
                groups.push_back( std::make_shared<weightGroup>( *group ) );
            }
        }
        initRwgt::initRwgt( xmlNode& wgtNode ) : xmlNode( wgtNode ){
            parser( true );
            name = "initrwgt";
            groups.reserve( children.size() );
            for(  auto child : children ){
                groups.push_back( std::make_shared<weightGroup>( *child ) );
            }
        }
        initRwgt::initRwgt( xmlNode* wgtNode ) : xmlNode( *wgtNode ){
            parser( true );
            name = "initrwgt";
            groups.reserve( children.size() );
            for(  auto child : children ){
                groups.push_back( std::make_shared<weightGroup>( *child ) );
            }
        }
        initRwgt::initRwgt( std::shared_ptr<xmlNode> wgtNode ) : xmlNode( *wgtNode ){
            parser( true );
            name = "initrwgt";
            groups.reserve( children.size() );
            for(  auto child : children ){
                groups.push_back( std::make_shared<weightGroup>( *child ) );
            }
        }
        initRwgt::initRwgt( xmlTree& wgtTree ) : xmlNode( wgtTree ){
            parser( true );
            name = "initrwgt";
            groups.reserve( children.size() );
            for(  auto child : children ){
                groups.push_back( std::make_shared<weightGroup>( *child ) );
            }
        }
        bool initRwgt::grpInit( std::shared_ptr<weightGroup>& wgt ){
            if( grpIsInit ){ return true; }
            else{
                groups = std::vector<std::shared_ptr<weightGroup>>( 1, wgt );
                grpIsInit = true;
                return false;
            }
        }
        void initRwgt::contWriter(){
            nodeContent = "\n";
            for( auto group : groups ){
                nodeContent += (*group->nodeWriter());
            }
        }
        void initRwgt::childWriter(){
            for( auto child : children ){
                if( child->getName() == "weightgroup" ){ continue; }
                nodeContent += (*child->nodeWriter());
            }
        }
        void initRwgt::childWriter( bool hasChildren ){
            if( hasChildren ){ childWriter(); }
            return;
        }

    // ZW: struct for handling weights
    // in event blocks of LHE files
        void bodyWgt::setComment( std::string_view nuComment ){ modded = true; comment = nuComment; }
        void bodyWgt::setVal( std::string nuVal ){ modded = true; valS = nuVal; valD = std::stod(valS);}
        void bodyWgt::setVal( std::string_view nuVal ){ modded = true; valS = std::string(nuVal); valD = std::stod(valS);}
        void bodyWgt::setVal( double nuVal ){ modded = true; valD = nuVal; valS = std::to_string(valD);}
        void bodyWgt::setId( std::string nuId ){ 
            modded = true; id = nuId;
            for( auto tag : tags ){
                if( tag->getId() == "id" ){ tag->setVal( id ); return; }
            }
            addTag( std::make_shared<xmlTag>( "id", id ) );
        }
        void bodyWgt::setModded( bool nuModded ){ modded = nuModded; }
        std::string_view bodyWgt::getComment(){ return comment; }
        std::string_view bodyWgt::getValS(){ return valS; }
        double bodyWgt::getValD(){ return valD; }
        bodyWgt::bodyWgt() : xmlNode(){ return; }
        bodyWgt::bodyWgt( std::string_view value ) : xmlNode() { setVal( value ); modded = false; }
        bodyWgt::bodyWgt( double value ) : xmlNode() { setVal( value ); modded = false; }
        bodyWgt::bodyWgt( std::string_view value, xmlTag rwgtId ) : xmlNode() { setVal( value ); addTag( std::make_shared<xmlTag>(rwgtId) ); modded = false; }
        bodyWgt::bodyWgt( double value, xmlTag rwgtId ) : xmlNode() { setVal( value ); addTag( std::make_shared<xmlTag>(rwgtId) ); modded = false; }
        bodyWgt::bodyWgt( std::string_view value, std::shared_ptr<xmlTag> rwgtId ) : xmlNode() { setVal( value ); addTag( rwgtId ); modded = false; }
        bodyWgt::bodyWgt( double value, std::shared_ptr<xmlTag> rwgtId ) : xmlNode() { setVal( value ); addTag( rwgtId ); modded = false; }
        bodyWgt::bodyWgt( const std::string_view originFile, const size_t& begin, const std::vector<std::shared_ptr<xmlNode>>& childs )
        : xmlNode( originFile, begin, childs ){
            auto strtPt = originFile.find_first_not_of(" >+", originFile.find(">", begin)+1);
            valS = originFile.substr( strtPt, originFile.find(" ", strtPt) - strtPt );
            valD = std::stod( valS );
        }
        bodyWgt::bodyWgt( xmlNode& wgtNode ) : xmlNode( wgtNode ){
            parser( true );
            valS = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            valD = std::stod( valS );
        }
        bodyWgt::bodyWgt( xmlNode* wgtNode ) : xmlNode( *wgtNode ){
            parser( true );
            valS = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            valD = std::stod( valS );
        }
        bodyWgt::bodyWgt( std::shared_ptr<xmlNode> wgtNode ) : xmlNode( *wgtNode ){
            parser( true );
            valS = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            valD = std::stod( valS );
        }
        bodyWgt::bodyWgt( xmlTree& wgtTree ) : xmlNode( wgtTree ){
            parser( true );
            valS = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            valD = std::stod( valS );
        }
        bodyWgt::bodyWgt( xmlTree* wgtTree ) : xmlNode( *wgtTree ){
            parser( true );
            valS = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            valD = std::stod( valS );
        }
        bodyWgt::bodyWgt( std::shared_ptr<xmlTree> wgtTree ) : xmlNode( *wgtTree ){
            parser( true );
            valS = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            valD = std::stod( valS );
        }
        bodyWgt::bodyWgt( double value, std::string& idTag ){
            setVal( value );
            id = idTag;
            addTag( std::make_shared<xmlTag>("id",id) );
        }
        void bodyWgt::appendWgt( std::shared_ptr<std::string> document ){
            if( !isWritten() ){ fullWriter(); }
            *document += *writtenSelf;
        }
        void bodyWgt::appendWgt( std::string* document ){
            if( !isWritten() ){ fullWriter(); }
            *document += *writtenSelf;
        }
        std::shared_ptr<std::string> bodyWgt::appendWgt( std::string_view document ){
            if(!isWritten() ){ fullWriter(); }
            auto retDoc = std::make_shared<std::string>( document );
            *retDoc += *writtenSelf;
            return retDoc;
        }
        void bodyWgt::fullWriter() {
            writtenSelf = std::make_shared<std::string>( "<wgt" );
            for( auto tag : tags ){
                *writtenSelf += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            *writtenSelf += ">" + std::string(valS) + "</wgt>\n";
            modded = false;
            written = true;
        }

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
        if( newBlock == npos ){ parLines = nuLineSplitter( parseFile.substr( blockStrt ) ); }
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
        std::string_view evHead::getComment(){ return comment; }
        std::string_view evHead::getWeight(){ return weight; }
        std::string_view evHead::getScale(){ return scale; }
        std::string_view evHead::getAQED(){ return aqed; }
        std::string_view evHead::getAQCD(){ return aqcd; }
        std::string_view evHead::getNprt(){ return nprt; }
        std::string_view evHead::getProcID(){ return procid; }
        bool evHead::isModded(){ return modded; }
        bool evHead::isWritten(){ return written; }
        void evHead::setComment( std::string_view nuCom ){ modded = true; comment = nuCom; }
        void evHead::setWeight( std::string_view nuWgt ){ modded = true; weight = nuWgt; }
        void evHead::setScale( std::string_view nuScale ){ modded = true; scale = nuScale; }
        void evHead::setAQED( std::string_view nuAQED ){ modded = true; aqed = nuAQED; }
        void evHead::setAQCD( std::string_view nuAQCD ){ modded = true; aqcd = nuAQCD; }
        void evHead::setNprt( std::string_view nuNprt ){ modded = true; nprt = nuNprt; }
        void evHead::setNprt( int nuNprt ){ modded = true; nprtint = nuNprt; nprtstr = std::to_string(nuNprt); nprt = nprtstr;}
        void evHead::setProcID( std::string_view nuProcID ){ modded = true; procid = nuProcID; }
        std::shared_ptr<std::string> evHead::getContent(){
            if( !isWritten() || isModded() ){ writer(); }
            return content;
        }
        evHead::evHead(){ return; }
        evHead::evHead( const std::string_view originFile, size_t beginLine, size_t endLine )
        {
            if( originFile.size() == 0){ return; }
            beginLine = originFile.find_first_not_of("\n \r\f\t\v", beginLine);
            if( endLine == npos ){ endLine = originFile.find("\n", beginLine ) + 1; }
            sourceFile = originFile.substr( beginLine, endLine - beginLine );
            auto evLine = nuWordSplitter( sourceFile );
            nprt = evLine->at(0) ;
            procid = evLine->at(1);
            weight = evLine->at(2);
            scale = evLine->at(3);
            aqed = evLine->at(4);
            aqcd = evLine->at(5);
         }
        void evHead::writer(){
            if( isWritten() && !isModded() ){ return; }
            if( !isModded() ){ content = std::make_shared<std::string>( sourceFile ); return; }
            auto retText = std::make_shared<std::string>( " " );
            *content = " " + std::string( nprt );
            for( size_t k = 0 ; k < 8 - procid.length() ; ++k ){ *content += " "; }
            *content +=  std::string( procid ) + " " + std::string( weight ) + " " + std::string( scale ) + " " + std::string( aqed ) + " " + std::string( aqcd );
            if( comment != "" ){ *content += " # " + std::string( comment ); }
            *content += "\n";
            modded = false;
            written = true;
        }

    // ZW: struct for handling particle lines
    // in LHE format event block
        std::string_view lhePrt::getLine(){ return sourceFile; }
        std::string_view lhePrt::getComment(){ return comment; }
        std::vector<std::string_view> lhePrt::getMom(){ return std::vector<std::string_view>( std::begin( mom ), std::end( mom ) ); }
        std::string_view lhePrt::getE(){ return energy; }
        std::string_view lhePrt::getMass(){ return mass; }
        std::string_view lhePrt::getVTim(){ return vtim; }
        std::string_view lhePrt::getSpin(){ return spin; }
        std::string_view lhePrt::getPDG(){ return pdg; }
        std::string_view lhePrt::getStatus(){ return status; }
        std::vector<std::string_view> lhePrt::getMothers(){ return std::vector<std::string_view>( std::begin( mothers ), std::end( mothers ) ); }
        std::vector<std::string_view> lhePrt::getColor(){ return std::vector<std::string_view>( std::begin( icol ), std::end( icol ) ); }
        void lhePrt::setComment( std::string_view nuCom ){ modded = true; comment = nuCom; }
        void lhePrt::setMom( std::vector<std::string_view> nuMom ){ modded = true; mom[0] = nuMom[0]; mom[1] = nuMom[1]; mom[2] = nuMom[2]; }
        void lhePrt::setEnergy( std::string_view nuE ){ modded = true; energy = nuE; }
        void lhePrt::setMass( std::string_view nuM ){ modded = true; mass = nuM; }
        void lhePrt::setVTim( std::string_view nuVTim ){ modded = true; vtim = nuVTim; }
        void lhePrt::setSpin( std::string_view nuSpin ){ modded = true; spin = nuSpin; }
        void lhePrt::setPDG( std::string_view nuPDG ){ modded = true; pdg = nuPDG; }
        void lhePrt::setStatus( std::string_view nuSt ){ modded = true; status = nuSt; }
        void lhePrt::setMothers( std::vector<std::string_view> nuMum ){ modded = true; mothers[0] = nuMum[0]; mothers[1] = nuMum[1]; }
        void lhePrt::setColors( std::vector<std::string_view> nuCol ){ modded = true; icol[0] = nuCol[0]; icol[1] = nuCol[1]; }
        bool lhePrt::isModded(){ return modded; }
        bool lhePrt::isWritten(){ return written; }
        std::shared_ptr<std::string> lhePrt::getContent(){
            if( !isWritten() || isModded() ){ writer(); }
            return content;
        }
        lhePrt::lhePrt(){ return; }
        lhePrt::lhePrt( std::pair<int,int>& prtInfo ){
            status = std::to_string( prtInfo.first );
            pdg = std::to_string( prtInfo.second );
        }
        lhePrt::lhePrt( std::pair<std::string,std::string>& prtInfo ){
            status = std::string_view( prtInfo.first );
            pdg = std::string_view( prtInfo.second );
        }
        lhePrt::lhePrt( const std::string_view originFile, const size_t& beginLine, const size_t& endLine )
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
        void lhePrt::writer(){
            if( isWritten() && !isModded() ){ return; }
            if( !isModded() ){ content = std::make_shared<std::string>( sourceFile ); return; }
            *content = "";
            for( size_t k = 0; k < 10 - pdg.length() ; ++k ){ *content += " "; }
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

    // ZW: struct for handling LHE format event block
        evHead event::getHead(){ return header; }
        std::vector<std::shared_ptr<lhePrt>> event::getPrts(){ return prts; }
        std::vector<std::shared_ptr<bodyWgt>> event::getWgts(){ return rwgt; }
        void event::setHead( evHead head ){ modded = true; header = head; }
        void event::addPrt( std::shared_ptr<lhePrt> prtcl ){ modded = true; prts.push_back( prtcl ); }
        void event::addPrt( lhePrt prtcl ){ modded = true; prts.push_back( std::make_shared<lhePrt>(prtcl) ); }
        void event::setPrts( std::vector<std::shared_ptr<lhePrt>> prtcls ){ modded = true; prts = prtcls; }
        void event::addWgt( bodyWgt nuWgt ){ addedWgt = true; rwgt.push_back( std::make_shared<bodyWgt>(nuWgt) ); }
        void event::addWgt( std::shared_ptr<bodyWgt> nuWgt ){ modded = true; rwgt.push_back( nuWgt ); }
        void event::addWgt( bodyWgt nuWgt, std::string& id ){ addedWgt = true; nuWgt.setId( id ); rwgt.push_back( std::make_shared<bodyWgt>(nuWgt) ); }
        void event::addWgt( std::shared_ptr<bodyWgt> nuWgt, std::string& id ){ modded = true; nuWgt->setId( id ); rwgt.push_back( nuWgt ); }
        bool event::newWeight(){ return addedWgt; }
        int event::getNprt(){ return prts.size(); }
        bool event::isModded() { return modded; }
        bool event::isModded( bool deep ) {
            if( !deep ){ return modded; }
            bool modStat = modded;
            for( auto child : children ){ if(modStat){ return modStat; }; modStat = (modStat || child->isModded( deep )); }
            modStat = (modStat || header.isModded());
            for( auto prt : prts ){ if(modStat){ return modStat; }; modStat = (modStat || prt->isModded()); }
            for( auto wgt : rwgt ){ if(modStat){ return modStat; }; modStat = (modStat || wgt->isModded()); }
            return modStat;
        }
        event::event(){ return; }
        event::event( std::vector<std::pair<int,int>>& prtInfo ){
            header.setNprt( std::to_string( prtInfo.size() ) );
            for( auto& prt : prtInfo ){
                prts.push_back( std::make_shared<lhePrt>( prt ) );
            }
        }
        event::event( std::vector<std::pair<std::string,std::string>>& prtInfo ){
            header.setNprt( prtInfo.size()  );
            for( auto& prt : prtInfo ){
                prts.push_back( std::make_shared<lhePrt>( prt ) );
            }
        }
        event::event( std::vector<std::shared_ptr<lhePrt>> prtInfo ){
            header.setNprt( std::to_string( prtInfo.size() ) );
            prts = prtInfo;
        }
        event::event( const std::string_view originFile, const size_t& begin, const std::vector<std::shared_ptr<xmlNode>>& childs ) 
        : xmlNode(originFile, begin, childs) {
            xmlFile = originFile; start = begin; children = childs; size_t trueStart = originFile.find_first_not_of(" \n\r\f\t\v", begin+1);
            if( trueStart == npos ){ return; }
            auto vals = lineFinder( originFile.substr( trueStart, originFile.find("<", trueStart +  3 ) - trueStart + 3 ));
            header = evHead(originFile, vals->at(0) + trueStart, vals->at(1) + trueStart + 1 );
            prts.reserve(vals->size());
            for( int k = 1 ; k < std::stoi(std::string(header.getNprt())) + 1; ++k)
            {
                prts.push_back( std::make_shared<lhePrt>(originFile, vals->at(k) + trueStart + 1, vals->at(k+1) + trueStart + 1) );
            }
        }
        event::event( const xmlNode& originFile )
        : xmlNode( originFile ) {
            size_t trueStart = xmlFile.find_first_not_of(" \n\r\f\t\v", start+1);
            auto vals = lineFinder( xmlFile.substr( trueStart, xmlFile.find("<", trueStart +  3 ) - trueStart + 3 ));
            header = evHead(xmlFile, vals->at(0) + trueStart, vals->at(1) + trueStart );
            prts.reserve(vals->size());
            for( int k = 1 ; k < std::stoi(std::string(header.getNprt())) + 1; ++k)
            {
                prts.push_back( std::make_shared<lhePrt>(xmlFile, vals->at(k) + trueStart + 1, vals->at(k+1) + trueStart) );
            }
        }
        event::event( const xmlNode* originFile )
        : xmlNode( *originFile ) {
            size_t trueStart = xmlFile.find_first_not_of(" \n\r\f\t\v", structure.getContStart() + 1);
            auto vals = lineFinder( xmlFile.substr( trueStart, xmlFile.find("<", trueStart +  3 ) - trueStart + 3 ));
            header = evHead(xmlFile, vals->at(0) + trueStart, vals->at(1) + trueStart );
            prts.reserve(vals->size());
            for( int k = 1 ; k < std::stoi(std::string(header.getNprt())) + 1; ++k)
            {
                prts.push_back( std::make_shared<lhePrt>(xmlFile, vals->at(k) + trueStart + 1, vals->at(k+1) + trueStart) );
            }
        }
        event::event( const std::shared_ptr<xmlNode>& originFile )
        : xmlNode( *originFile ) {
            size_t trueStart = xmlFile.find_first_not_of(" \n\r\f\t\v", structure.getContStart() + 1);
            auto vals = lineFinder( xmlFile.substr( trueStart, xmlFile.find("<", trueStart +  3 ) - trueStart + 3 ));
            header = evHead(xmlFile, vals->at(0) + trueStart, vals->at(1) + trueStart );
            prts.reserve(vals->size());
            for( int k = 1 ; k < std::stoi(std::string(header.getNprt())) + 1; ++k)
            {
                prts.push_back( std::make_shared<lhePrt>(xmlFile, vals->at(k) + trueStart + 1, vals->at(k+1) + trueStart) );
            }
        }
        event::event( xmlTree& originFile )
        : xmlNode( originFile ) {
            size_t trueStart = xmlFile.find_first_not_of(" \n\r\f\t\v", structure.getContStart() + 1);
            auto vals = lineFinder( xmlFile.substr( trueStart, xmlFile.find("<", trueStart +  3 ) - trueStart + 3 ));
            header = evHead(xmlFile, vals->at(0) + trueStart, vals->at(1) + trueStart );
            prts.reserve(vals->size());
            for( int k = 1 ; k < std::stoi(std::string(header.getNprt())) + 1; ++k)
            {
                prts.push_back( std::make_shared<lhePrt>(xmlFile, vals->at(k) + trueStart + 1, vals->at(k+1) + trueStart) );
            }
        }
        event::event( xmlTree* originFile )
        : xmlNode( *originFile ) {
            size_t trueStart = xmlFile.find_first_not_of(" \n\r\f\t\v", structure.getContStart() + 1);
            auto vals = lineFinder( xmlFile.substr( trueStart, xmlFile.find("<", trueStart +  3 ) - trueStart + 3 ));
            header = evHead(xmlFile, vals->at(0) + trueStart, vals->at(1) + trueStart );
            prts.reserve(vals->size());
            for( int k = 1 ; k < std::stoi(std::string(header.getNprt())) + 1; ++k)
            {
                prts.push_back( std::make_shared<lhePrt>(xmlFile, vals->at(k) + trueStart + 1, vals->at(k+1) + trueStart) );
            }
        }
        event::event( std::shared_ptr<xmlTree> originFile )
        : xmlNode( *originFile ) {
            size_t trueStart = xmlFile.find_first_not_of(" \n\r\f\t\v", structure.getContStart() + 1);
            auto vals = lineFinder( xmlFile.substr( trueStart, xmlFile.find("<", trueStart +  3 ) - trueStart + 3 ));
            header = evHead(xmlFile, vals->at(0) + trueStart, vals->at(1) + trueStart );
            prts.reserve(vals->size());
            for( int k = 1 ; k < std::stoi(std::string(header.getNprt())) + 1; ++k)
            {
                prts.push_back( std::make_shared<lhePrt>(xmlFile, vals->at(k) + trueStart + 1, vals->at(k+1) + trueStart) );
            }
        }
        bool event::prtsAreMod(){
            for( auto prt : prts ){ if( prt->isModded() ){ return true; } }
            return false;
        }
        bool event::headIsMod(){
            return header.isModded();
        }
        bool event::isSpecSort() const { return specSorted; }
        sortFcn event::getSortFcn() const { return eventSort; }
        statSort event::getStatSort() const { return specSort; }
        bool event::hasRwgt(){
            if( rwgt.size() > 0 ){ return true; }
            return false;
        }
        bool event::rwgtChild(){
            if( childRwgt != nullptr ){ return true; }
            for( auto child : children ){ if( clStringComp(child->getName(), std::string("rwgt") ) ){ childRwgt = child; return true; } }
            return false;
        }
        bool event::bothRwgt(){ return (hasRwgt() && rwgtChild() ); }
        bool event::eitherRwgt(){ return (hasRwgt() || rwgtChild() ); }
        bool event::initProcMap(bool hard)
        {
            if(!hard){ if( procMap.size() > 0 ){ return true; } }
            for( auto prt : prts ){
                procMap.insert({prt->getStatus(), std::vector<std::string_view>()});
                procOrder.insert({prt->getStatus(), std::vector<size_t>()});
            }
            for( auto prt : prts ){
                procMap[prt->getStatus()].push_back( prt->getPDG() );
            }
            for( auto stat = procMap.begin(); stat!= procMap.end(); ++stat ){
                procOrder[stat->first] = *stoiSort( stat->second );
            }
            hasBeenProc = true;
            return true;
        }
        bool event::initProcMap( sortFcn sorter, bool hard )
        {
            if(!hard){ if( procMap.size() > 0 ){ return true; } }
            specSorted = false;
            eventSort = sorter;
            for( auto prt : prts ){
                procMap.insert({prt->getStatus(), std::vector<std::string_view>()});
                procOrder.insert({prt->getStatus(), std::vector<size_t>()});
            }
            for( auto prt : prts ){
                procMap[prt->getStatus()].push_back( prt->getPDG() );
            }
            for( auto stat = procMap.begin(); stat!= procMap.end(); ++stat ){
                procOrder[stat->first] = *sorter( stat->second );
            }
            hasBeenProc = true;
            return true;
        }
        bool event::initProcMap( statSort sorter, bool hard )
        {
            if(!hard){ if( procMap.size() > 0 ){ return true; } }
            specSorted = true;
            specSort = sorter;
            for( auto prt : prts ){
                procMap.insert({prt->getStatus(), std::vector<std::string_view>()});
                procOrder.insert({prt->getStatus(), std::vector<size_t>()});
            }
            for( auto prt : prts ){
                procMap[prt->getStatus()].push_back( prt->getPDG() );
            }
            for( auto stat = procMap.begin(); stat!= procMap.end(); ++stat ){
                procOrder[stat->first] = *sorter(stat->first, stat->second );
            }
            hasBeenProc = true;
            return true;
        }
        bool event::inRwgtChild( std::string_view nameIn ){ 
            for( auto child : childRwgt->getChildren() ){ 
                for( auto tag : child->getTags() ){ if(clStringComp(tag->getVal(), nameIn)){ return true; } }
            }
            return false;
        }
        bool event::checkRwgtOverlap(){
            for( auto wgt : rwgt ){ 
                for( auto tag : wgt->getTags() ){ if( inRwgtChild( tag->getVal() ) ){ return true; } }
            }
            return false;
        }
        void event::childRwgtWriter(){
            if( rwgtChild() ){ nodeContent += *childRwgt->nodeWriter(); }
        }
        void event::vecRwgtWriter( bool midNode ){
            if( !midNode ){ nodeContent += "<rwgt>\n"; }
            for( auto wgt : rwgt ){ 
                nodeContent += *wgt->nodeWriter();
            }
            nodeContent += "</rwgt>\n";
        }
        void event::rwgtWriter(){
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
        void event::contWriter() {
            nodeContent = "\n" + *header.getContent();
            for( auto prt : prts ){
                nodeContent += *prt->getContent();
            }
        }
        void event::childWriter() {
            for( auto child : children ){
                if( clStringComp( child->getName(), std::string("wgt") ) ){ continue; }
                nodeContent += *child->nodeWriter();
            }
        }
        void event::fullWriter() {
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
        void event::fullWriter( bool deep ){
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
        void event::appendWgts(){
            if( !addedWgt ){ return; }
            writtenSelf->erase( writtenSelf->size() - 17, 17 );
            for( auto wgt : rwgt ){
                if( !wgt->isWritten() ){ wgt->appendWgt( writtenSelf ); }
            }
            *writtenSelf += "</rwgt>\n</event>\n";
        }
        std::shared_ptr<std::string> event::nodeWriter() {
            if( isModded(false) || !isWritten() ){ fullWriter(); return writtenSelf; }
            if( addedWgt ){ appendWgts(); }
            return writtenSelf;
        }
        std::shared_ptr<std::string> event::nodeWriter( bool recursive ){
            if( isModded( recursive ) || !isWritten() ){ fullWriter(); return writtenSelf; }
            if( addedWgt ){ appendWgts(); }
            return writtenSelf;
        }
        std::map<std::string_view, std::vector<std::string_view>> &event::getProc(){
            if( initProcMap() ){ return procMap; }
            else throw std::runtime_error("Error while parsing event node.");
        }
        std::map<std::string_view, std::vector<size_t>> &event::getProcOrder(){
            if( initProcMap() ){ return procOrder; }
            else throw std::runtime_error("Error while parsing event node.");
        }
        std::map<std::string_view, std::vector<std::string_view>> event::getProc() const {
            if ( hasBeenProc ){ return procMap; }
            else throw std::runtime_error("Const declaration of event node before it has been procesed.");
        }
        std::map<std::string_view, std::vector<size_t>> event::getProcOrder() const {
            if ( hasBeenProc ){ return procOrder; }
            else throw std::runtime_error("Const declaration of event node before it has been procesed.");
        }
        std::map<std::string_view, std::vector<std::string_view>> &event::getProc(sortFcn sorter){
            if( initProcMap(sorter) ){ return procMap; }
            else throw std::runtime_error("Error while parsing event node.");
        }
        std::map<std::string_view, std::vector<size_t>> &event::getProcOrder(sortFcn sorter){
            if( initProcMap(sorter) ){ return procOrder; }
            else throw std::runtime_error("Error while parsing event node.");
        }
        std::map<std::string_view, std::vector<std::string_view>> &event::getProc(statSort sorter){
            if( initProcMap(sorter) ){ return procMap; }
            else throw std::runtime_error("Error while parsing event node.");
        }
        std::map<std::string_view, std::vector<size_t>> &event::getProcOrder(statSort sorter){
            if( initProcMap(sorter) ){ return procOrder; }
            else throw std::runtime_error("Error while parsing event node.");
        }

    event& makeEv( std::vector<std::pair<int,int>>& particles ){
        static auto returnEvent = event( particles );
        return returnEvent;
    }

    std::vector<std::shared_ptr<lhePrt>> getParticles( event& ev ){
        return ev.getPrts();
    }

    // ZW: struct for handling the first line of
    // LHE format init tag
        bool lheInitHead::isWritten(){ return written; }
        bool lheInitHead::isModded(){ return modded; }
        std::shared_ptr<std::string> lheInitHead::getContent(){ 
            if( isModded() || !isWritten() ){ writer(); }
            return content; }
        lheInitHead::lheInitHead( std::string_view initHead ){
            auto vals = *nuBlankSplitter( initHead );
            if( vals.size() < 10 ){ return; }
            idbmup[0] = vals[0]; idbmup[1] = vals[1];
            ebmup[0] = vals[2]; ebmup[1] = vals[3];
            pdfgup[0] = vals[4]; pdfgup[1] = vals[5];
            pdfsup[0] = vals[6]; pdfsup[1] = vals[7];
            idwtup = vals[8]; nprup = vals[9];
        }
        lheInitHead::lheInitHead( xmlNode& initNode )
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
        void lheInitHead::writer(){
            *content = std::string(idbmup[0]) + " " + std::string(idbmup[1]) + " " + std::string(ebmup[0]) + " " + std::string(ebmup[1]) + " " + std::string(pdfgup[0]) 
    + " " + std::string(pdfgup[1]) + " " + std::string(pdfsup[0]) + " " + std::string(pdfsup[1]) + " " + std::string(idwtup) + " " + std::string(nprup) +"\n";
            written = true;
            modded = false;
        }

    // ZW: struct for handling process lines
    // in LHE format init tag
        bool lheInitLine::isWritten(){ return written; }
        bool lheInitLine::isModded(){ return modded; }
        std::shared_ptr<std::string> lheInitLine::getContent(){ 
            if( isModded() || !isWritten() ){ writer(); }
            return content; }
        lheInitLine::lheInitLine(){}
        lheInitLine::lheInitLine( std::string_view procLine )
        {
            auto vals = *nuBlankSplitter( procLine );
            if( vals.size() < 4 ){ return; }
            xsecup = vals[0];
            xerrup = vals[1];
            xmaxup = vals[2];
            lprup = vals[3];
        }
        void lheInitLine::writer(){
            *content = std::string(xsecup) + " " + std::string(xerrup) + " " + std::string(xmaxup) + " " + std::string(lprup) + "\n";
            written = true;
            modded = false;
        }

    // ZW: struct for handling single parameter line in
    // SLHA format parameter card
        void paramVal::parse(){
            id = std::stoi( std::string(idStr) );
            value = std::stod( std::string(valStr) );
        }
        paramVal::paramVal(){ realLine = ""; idStr = ""; valStr = ""; }
        paramVal::paramVal( std::string_view paramLine, bool parseOnline )
        {
            if( paramLine.find("\n") != npos ){
                auto startPos = paramLine.find_first_not_of(" \n", paramLine.find("\n"));
                if( startPos!= npos ){
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
        bool paramVal::isMod(){ return modded; }
        std::shared_ptr<std::string> paramVal::selfWrite(){
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

    // ZW: struct for handling single DECAY line
    // in SLHA format parameter card
        void decVal::parse() {
            auto vals = *nuBlankSplitter( realLine );
            id = std::stoi( std::string(vals[1]) );
            value = std::stod( std::string(vals[2]) );
            if( vals.size() > 3 )
            {
                auto comStart = realLine.find("#");
                comment = realLine.substr( comStart, realLine.find("\n", comStart) - comStart );
            }
        }
        decVal::decVal( std::string_view paramLine, bool parseOnline ) : paramVal( paramLine, false )
        {
            if( parseOnline ){ parse(); }
        }
        std::shared_ptr<std::string> decVal::selfWrite() {
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

    // ZW: struct for handling parameter block
    // in SLHA format parameter card
        void paramBlock::parse( bool parseOnline ){
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
        paramBlock::paramBlock(){ return; }
        paramBlock::paramBlock( std::string_view paramSet, bool parseOnline )
        {
            realBlock = paramSet;
            startPt = clStringFind( realBlock, std::string("\nB") );
            if( parseOnline ){ parse(parseOnline);  }
        }
        bool paramBlock::isMod(){ return modded; }
        std::shared_ptr<std::string> paramBlock::selfWrite(){
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
            else{ if( startPt == npos ){
                *writeBlock += realBlock;
            } else {
                *writeBlock = realBlock.substr( startPt );
            } }
            return writeBlock;
        }

    // ZW: struct for handling DECAY lines
    // in SLHA format parameter card
        void decBlock::parse( bool parseOnline ){
        if( realBlock.size() == 0 ){ return; }
            auto decLines = clFindEach( realBlock, std::string("\ndecay") );
            decays.reserve(decLines->size());
            if( realBlock.size() > 5 ){  if( clStringComp( realBlock.substr(0,5), std::string("decay")) )
            { decays.push_back( decVal(realBlock.substr( 0, realBlock.find("\n") ), parseOnline) ); } }
            for( auto pts : *decLines )
            {
                auto lineBr = realBlock.find( "\n", pts + 1 );
                if( lineBr == npos ){ decays.push_back( decVal( realBlock.substr( pts + 1), parseOnline ) ); continue; }
                decays.push_back( decVal( realBlock.substr( pts + 1, lineBr - pts - 1 ), parseOnline ) );
            }
        }
        void decBlock::parse( std::shared_ptr<std::vector<size_t>> decLines, bool parseOnline ) {
            decays.reserve(decLines->size());
            if( realBlock.size() > 5 ){ if( clStringComp( realBlock.substr(0,5), std::string("decay")) )
            { decays.push_back( decVal(realBlock.substr( 0, realBlock.find("\n") ), parseOnline) ); } }
            for( auto pts : *decLines )
            {
                auto lineBr = realBlock.find( "\n", pts + 1 );
                if( lineBr == npos ){ decays.push_back( decVal( realBlock.substr( pts + 1), parseOnline ) ); continue; }
                decays.push_back( decVal( realBlock.substr( pts + 1, lineBr - pts - 1 ), parseOnline ) );
            }
        }
        decBlock::decBlock( std::string_view paramSet, bool parseOnline ) : paramBlock( paramSet, parseOnline )
        {
            realBlock = paramSet;
            if( parseOnline ){ parse(parseOnline); }
        }
        std::shared_ptr<std::string> decBlock::selfWrite() {
            auto writeBlock = std::make_shared<std::string>("");
            *writeBlock += "\n";
            for ( auto val : decays )
            {
                *writeBlock += *val.selfWrite(); 
            }
            return writeBlock;
        }

    // ZW: struct for handling SLHA parameter cards
         void lesHouchesCard::parse( bool parseOnline )
        {
            if( parsed ){ return; }
            if( xmlFile.substr(start,1).find_first_of("BbDd#") == npos ){ start = clStringFindIf( xmlFile, std::string("\n"), lambdaNu ); }
            auto blockPts = clFindEach( xmlFile, std::string("\nblock") );
            auto decLines = clFindEach( xmlFile, std::string("\ndecay") );
            header = xmlFile.substr( start, std::min( blockPts->at(0), decLines->at(0) ) - start );
            for( size_t k  = 0 ; k < blockPts->size() - 1 ; ++k )
            {
                blocks.push_back( paramBlock( xmlFile.substr( blockPts->at(k), blockPts->at(k+1) - blockPts->at(k) ), parseOnline ) );
            }
            blocks.push_back(paramBlock(xmlFile.substr(blockPts->at(blockPts->size()-1), clStringFindIf( xmlFile, std::string("\n"),
             lambda, blockPts->at(blockPts->size()-1) + 1) - blockPts->at(blockPts->size()-1)), parseOnline));
            decays = decBlock( xmlFile );
            decays.parse( decLines, parseOnline );
            parsed = true;
        } 
        lesHouchesCard::lesHouchesCard( const std::string_view originFile, const size_t& begin, bool parseOnline ){ 
            xmlFile = originFile; start = begin;
            modded = false; blockStart = clStringFindIf( xmlFile, std::string("\n"), lambda, start + 1); end = xmlFile.find("</", blockStart);
            parsed = false;
            if( parseOnline ){ parse( parseOnline ); }
        }
        bool lesHouchesCard::isMod(){ return modded; }
        std::shared_ptr<std::string> lesHouchesCard::selfWrite(){
            auto writeCard = std::make_shared<std::string>(header);
            if( isMod() )
            { for( auto block : blocks )
                { *writeCard += *block.selfWrite(); }
                *writeCard += *decays.selfWrite(); }
            else{
                if( end != npos ){ *writeCard += std::string( xmlFile.substr( blockStart, end - blockStart ) );
                } else{ *writeCard += std::string( xmlFile.substr( blockStart ) ); }
            }
            return writeCard;
        }

        std::shared_ptr<lesHouchesCard> slhaNode::getParameters(){
            modded = true;
            return parameterCard;
        }
        slhaNode::slhaNode() : xmlNode(){}
        slhaNode::slhaNode( lesHouchesCard parameters ) : xmlNode(){
            parameterCard = std::make_shared<lesHouchesCard>( parameters );
            pCardInit = true;
        } 
        slhaNode::slhaNode( std::shared_ptr<lesHouchesCard> parameters ) : xmlNode(){
            parameterCard = parameters;
            pCardInit = true;
        }
        slhaNode::slhaNode( xmlNode& node, bool parseOnline ) : xmlNode( node ){ 
            parameterCard = std::make_shared<lesHouchesCard>( node.getFile(), node.getStart(), parseOnline );
        }
        slhaNode::slhaNode( xmlNode* node, bool parseOnline ) : xmlNode( *node ){ 
            parameterCard = std::make_shared<lesHouchesCard>( node->getFile(), node->getStart(), parseOnline );
        }
        slhaNode::slhaNode( std::shared_ptr<xmlNode> node, bool parseOnline ) : xmlNode( *node ){ 
            parameterCard = std::make_shared<lesHouchesCard>( node->getFile(), node->getStart(), parseOnline );
        }
        slhaNode::slhaNode( xmlTree tree, bool parseOnline ) : xmlNode( tree ){
            parameterCard = std::make_shared<lesHouchesCard>( tree.getOrigin(), tree.getStart(), parseOnline );
        }
        slhaNode::slhaNode( std::shared_ptr<xmlTree> tree, bool parseOnline ) : xmlNode( *tree ){
            parameterCard = std::make_shared<lesHouchesCard>( tree->getOrigin(), tree->getStart(), parseOnline );
        }
        slhaNode::slhaNode( xmlTree* tree, bool parseOnline ) : xmlNode( *tree ){
            parameterCard = std::make_shared<lesHouchesCard>( tree->getOrigin(), tree->getStart(), parseOnline );
        }
        slhaNode::slhaNode( const std::string_view originFile, const size_t& begin, bool parseOnline )
        : xmlNode( originFile, begin ){
            if( parse() ){ parameterCard = std::make_shared<lesHouchesCard>( content, begin, parseOnline ); pCardInit = true; }
        }
        void slhaNode::headWriter(){
            nodeHeader = "<slha";
            for( auto tag : tags ){
                nodeHeader += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            nodeHeader += ">";
        }
        void slhaNode::endWriter(){ nodeEnd += "</slha>\n"; }
        void slhaNode::contWriter(){
            if( pCardInit ){
                nodeContent = *parameterCard->selfWrite();
            } else {
                nodeContent = content;
            }
        }

    // ZW: struct for handling LHE init nodes
        std::shared_ptr<lheInitHead> initNode::getHead(){ return initHead; }
        std::vector<std::shared_ptr<lheInitLine>> initNode::getLines(){ return initLines; }
        void initNode::setHead( std::shared_ptr<lheInitHead> head ){ modded = true; initHead = head; }
        void initNode::setLines( std::vector<std::shared_ptr<lheInitLine>> lines ){ modded = true; initLines = lines; initHead->nprup = std::to_string( initLines.size() ); }
        void initNode::addLine( std::shared_ptr<lheInitLine> line ){ modded = true; initLines.push_back( line ); initHead->nprup = std::to_string( initLines.size() ); }
        initNode::initNode() : xmlNode(){ name = "init"; }
        initNode::initNode( const std::string_view originFile, const size_t& begin, bool parseOnline )
        : xmlNode( originFile, begin ){
            content = originFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            if( parseOnline ){ parse( parseOnline ); }
        }
        initNode::initNode( xmlNode& node, bool parseOnline ) : xmlNode( node ){
            content = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            if( parseOnline ){ parse( parseOnline ); }
        }
        initNode::initNode( xmlNode* node, bool parseOnline ) : xmlNode( *node ){
            content = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            if( parseOnline ){ parse( parseOnline ); }
        }
        initNode::initNode( std::shared_ptr<xmlNode> node, bool parseOnline ) : xmlNode( *node ){
            content = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            if( parseOnline ){ parse( parseOnline ); }
        }
        initNode::initNode( xmlTree tree, bool parseOnline ) : xmlNode( tree ){
            content = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            if( parseOnline ){ parse( parseOnline ); }
        }
        initNode::initNode( std::shared_ptr<xmlTree> tree, bool parseOnline ) : xmlNode( *tree ){
            content = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            if( parseOnline ){ parse( parseOnline ); }
        }
        initNode::initNode( xmlTree* tree, bool parseOnline ) : xmlNode( *tree ){
            content = xmlFile.substr( structure.getContStart(), structure.getContEnd() - structure.getContStart() );
            if( parseOnline ){ parse( parseOnline ); }
        }
        bool initNode::parseContent(){
            if( content.size() == 0 ){ return false; }
            auto linebreaks = lineFinder( content );
            if( linebreaks->size() == 0 ){ return false; }
            initHead = std::make_shared<lheInitHead>(content.substr( 0, linebreaks->at(0) ) );
            for( size_t k = 0 ; k < linebreaks->size() - 1 ; ++k ){
                initLines.push_back( std::make_shared<lheInitLine>( content.substr( linebreaks->at(k), linebreaks->at(k+1) - linebreaks->at(k) ) ) );
            }
            return true;
        }
        void initNode::contWriter(){
            if( isModded() ){nodeContent = std::string( content ); return; }
            nodeContent = *initHead->getContent();
            for( auto line : initLines ){
                nodeContent += *line->getContent();
            }
        }
    
    // ZW: struct for explicitly handling LHE header nodes
        size_t lheHead::addWgtGroup( std::shared_ptr<weightGroup>& wgtGroup ){
            hasRwgt = true; 
            modded = true; 
            if( wgtGrpInit( wgtGroup ) ){
                rwgtNodes->addGroup( wgtGroup );
            }
            return (rwgtNodes->noGrps() - 1);
        }
        size_t lheHead::addWgtGroup( weightGroup wgtGroup ){
            hasRwgt = true;
            modded = true;
            auto wgtGrpPtr = std::make_shared<weightGroup>( wgtGroup );
            if( wgtGrpInit( wgtGrpPtr ) ){
                rwgtNodes->addGroup( std::make_shared<weightGroup>( wgtGroup ) );
            }
            return (rwgtNodes->noGrps() - 1);
        }
        void lheHead::addWgt( size_t index, std::shared_ptr<headWeight> nuWgt ){
            if( index >= (size_t)rwgtNodes->getGroups().size() )
                throw std::range_error( "Appending weight to uninitialised weightgroup." );
            hasRwgt = true;
            modded = true;
            rwgtNodes->addWgt( index, nuWgt );
        }
        void lheHead::addWgt( size_t index, headWeight nuWgt ){
            if( index >= (size_t)rwgtNodes->getGroups().size() )
                throw std::range_error( "Appending weight to uninitialised weightgroup." );
            hasRwgt = true;
            modded = true;
            rwgtNodes->addWgt( index, nuWgt );
        }
        void lheHead::addWgt( size_t index, std::shared_ptr<headWeight> nuWgt, std::string idTagg ){
            if( index >= (size_t)rwgtNodes->getGroups().size() )
                throw std::range_error( "Appending weight to uninitialised weightgroup." );
            hasRwgt = true;
            modded = true;
            nuWgt->setId( idTagg );
            rwgtNodes->addWgt( index, nuWgt );
        }
        void lheHead::addWgt( size_t index, headWeight nuWgt, std::string idTagg ){
            if( index >= (size_t)rwgtNodes->getGroups().size() )
                throw std::range_error( "Appending weight to uninitialised weightgroup." );
            hasRwgt = true;
            modded = true;
            nuWgt.setId( idTagg );
            rwgtNodes->addWgt( index, nuWgt );
        }
        void lheHead::setInitRwgt( initRwgt initWgt ){  hasRwgt = true; modded = true; rwgtNodes = std::make_shared<initRwgt>(initWgt); }
        void lheHead::setInitRwgt( std::shared_ptr<initRwgt> initWgt ){ hasRwgt = true; modded = true; rwgtNodes = initWgt; }
        std::vector<std::shared_ptr<weightGroup>> lheHead::getWgtGroups(){ return rwgtNodes->getGroups(); }
        std::shared_ptr<initRwgt> lheHead::getInitRwgt(){ return rwgtNodes; }
        std::shared_ptr<slhaNode> lheHead::getParameters(){ return parameters; }
        void lheHead::setParameters( std::shared_ptr<slhaNode> params ){ parameters = params; }
        bool lheHead::rwgtInc(){ return hasRwgt; }
        lheHead::lheHead(){ return; }
        lheHead::lheHead( const std::string_view originFile, const size_t& begin, const std::vector<std::shared_ptr<xmlNode>>& childs )
        : xmlNode(originFile, begin, childs){
            xmlFile = originFile; start = begin; children = childs; size_t trueStart = originFile.find_first_not_of(" ", begin+1);
            if( trueStart != npos ){name = originFile.substr( trueStart, originFile.find_first_of(">/ ", trueStart) - trueStart );}
            for( auto child : children ){
                if (child->getName() == "slha" ){ parameters = std::make_shared<slhaNode>( *child ); continue; }
                if (child->getName() == "initrwgt" ){ rwgtNodes = std::make_shared<initRwgt>( *child ); continue; }
            }
        }
        lheHead::lheHead( xmlNode& node ) : xmlNode(node){
            for( auto child : node.getChildren() ){
                if ( child->getName() == "slha" ){ parameters = std::make_shared<slhaNode>( *child ); continue; }
                if ( child->getName() == "initrwgt" ){ rwgtNodes = std::make_shared<initRwgt>( *child ); continue; }
            }
        }
        lheHead::lheHead( xmlNode* node ) : xmlNode(*node){
            for( auto child : node->getChildren() ){
                if ( child->getName() == "slha" ){ parameters = std::make_shared<slhaNode>( *child ); continue; }
                if ( child->getName() == "initrwgt" ){ rwgtNodes = std::make_shared<initRwgt>( *child ); continue; }
            }
        }
        lheHead::lheHead( std::shared_ptr<xmlNode> node ) : xmlNode( *node ){
            for( auto child : node->getChildren() ){
                if ( child->getName() == "slha" ){ parameters = std::make_shared<slhaNode>( *child ); continue; }
                if ( child->getName() == "initrwgt" ){ rwgtNodes = std::make_shared<initRwgt>( *child ); continue; }
            }
        }
        lheHead::lheHead( xmlTree tree ) : xmlNode( tree ){
            for( auto child : children ){
                if ( child->getName() == "slha" ){ parameters = std::make_shared<slhaNode>( *child ); continue; }
                if ( child->getName() == "initrwgt" ){ rwgtNodes = std::make_shared<initRwgt>( *child ); continue; }
            }
        }
        lheHead::lheHead( std::shared_ptr<xmlTree> tree ) : xmlNode( *tree ){
            for( auto child : children ){
                if ( child->getName() == "slha" ){ parameters = std::make_shared<slhaNode>( *child ); continue; }
                if ( child->getName() == "initrwgt" ){ rwgtNodes = std::make_shared<initRwgt>( *child ); continue; }
            }
        }
        lheHead::lheHead( xmlTree* tree ) : xmlNode( *tree ){
            for( auto child : children ){
                if ( child->getName() == "slha" ){ parameters = std::make_shared<slhaNode>( *child ); continue; }
                if ( child->getName() == "initrwgt" ){ rwgtNodes = std::make_shared<initRwgt>( *child ); continue; }
            }
        }
        bool lheHead::wgtGrpInit( std::shared_ptr<weightGroup>& wgtGrp ){
            if( wgtGrpIsInit ){ return true; }
            if( rwgtNodes == nullptr ){
                rwgtNodes = std::make_shared<initRwgt>();
                wgtGrpIsInit = true;
                rwgtNodes->addGroup( wgtGrp );
                return false;
            } else throw std::runtime_error( "Error while initiating return LHE file header (initrwgt node is defined in an unrecognised manner)." );
        }
        void lheHead::setRelChild(){
            if( relChildSet ){ return; }
            relChild.reserve( children.size() );
            for( size_t k = 0 ; k < children.size() ; ++k ){
                auto child = &children[k];
                if( (*child)->getName() == "slha" ){ continue; }
                if( (*child)->getName() == "initrwgt" ){ continue; }
                relChild.push_back( k );
            }
            relChildSet = true;
        }
        bool lheHead::parseChildren( bool recursive ){
            bool status = true;
            for( auto child : children ){
                if( child->getName() == "slha" || child->getName() == "initrwgt" ){ continue; }
                child->parser( recursive );
                status = (status && child->isParsed() );
                deepParsed = true;
            }
            return status;
        }
        void lheHead::headWriter(){
            nodeHeader =  "<header";
            for( auto tag : tags ){
                nodeHeader += " " + std::string(tag->getId()) + "=\"" + std::string(tag->getVal()) + "\"";
            }
            nodeHeader += ">\n";
        }
        void lheHead::childWriter(){
            setRelChild();
            for( auto relKid : relChild ){
                nodeContent += *(children[relKid]->nodeWriter());
            }
            if( parameters != nullptr ){ nodeContent += *parameters->nodeWriter(); }
            if( hasRwgt ){ 
                nodeContent += *rwgtNodes->nodeWriter();
            }
        }
        void lheHead::fullWriter(){
            if( isModded() ){
            headWriter();
            contWriter();
            childWriter();
            endWriter();
            writtenSelf = std::make_shared<std::string>( nodeHeader + nodeContent + nodeEnd );
            written = true;
            }
        }

    // ZW: struct for keeping track of appended weights in LHE node,
    // since weight information is stored both in the header 
    // and in the individual events
        newWgt::newWgt( std::shared_ptr<headWeight> heaWgt, std::vector<std::shared_ptr<bodyWgt>> bodWgts ){
            headWgt = heaWgt; bodyWgts = bodWgts;
        }
        newWgt::newWgt( std::shared_ptr<headWeight> heaWgt, std::shared_ptr<std::vector<double>> wgts ){
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
        newWgt::newWgt( std::string_view parameters, std::shared_ptr<std::vector<double>> wgts, std::string idTag ){
            headWgt = std::make_shared<headWeight>(parameters, idTag);
            bodyWgts = std::vector<std::shared_ptr<bodyWgt>>(wgts->size());
            for( size_t i = 0 ; i < wgts->size() ; ++i ){
                bodyWgts[i] = std::make_shared<bodyWgt>(wgts->at(i), idTag);
            }
        }
        newWgt::newWgt( std::string_view parameters, int idNum, std::shared_ptr<std::vector<double>> wgts, std::string idTag  ){
            std::string newTag = std::string( idTag ) + "_" + std::to_string( idNum );
            headWgt = std::make_shared<headWeight>(parameters, newTag);
            bodyWgts = std::vector<std::shared_ptr<bodyWgt>>(wgts->size());
            for( size_t i = 0 ; i < wgts->size() ; ++i ){
                bodyWgts[i] = std::make_shared<bodyWgt>(wgts->at(i), newTag);
            }
        }
        newWgt::newWgt( std::string& parameters ){
            headWgt = std::make_shared<headWeight>(parameters);
        }
        newWgt::newWgt( std::string& parameters, std::string& idTag ){
            headWgt = std::make_shared<headWeight>(parameters, idTag);
        }
        std::shared_ptr<headWeight> newWgt::getHeadWgt(){ return headWgt; }
        std::vector<std::shared_ptr<bodyWgt>> newWgt::getBodyWgts(){ return bodyWgts; }
        void newWgt::addBdyWgts( std::shared_ptr<std::vector<double>> wgts ){
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

    // ZW: general struct for handling LHE files explicitly
        lheNode::lheNode() : xmlNode(){}
        lheNode::lheNode( const std::string_view originFile, const size_t& begin, const std::vector<std::shared_ptr<xmlNode>>& childs )
        : xmlNode(originFile, begin, childs){
            //xmlFile = originFile; start = begin; children = childs; size_t trueStart = originFile.find_first_not_of(" ", begin+1);
            //if( trueStart != npos ){name = originFile.substr( trueStart, originFile.find_first_of(">/ ", trueStart) - trueStart );}
            for( auto child : children ){
                if( child->getName() == "header" ){ header = std::make_shared<lheHead>( *child ); continue; }
                if( child->getName() == "init" ){ init = std::make_shared<initNode>( *child ); continue; }
                if( child->getName() == "event" ){ events.push_back( std::make_shared<event>( *child ) ); continue; }
            }
        }
        std::shared_ptr<lheHead> lheNode::getHeader(){ return header; }
        std::shared_ptr<initNode> lheNode::getInit(){ return init; }
        std::vector<std::shared_ptr<event>> lheNode::getEvents(){ return events; }
        bool lheNode::isModded(){ return modded; }
        bool lheNode::isModded( bool deep ){
            if( !deep ){ return isModded(); }
            bool modStat = isModded();
            for( auto child : children ){ modStat = ( modStat || child->isModded( deep ) ); }
            for( auto event : events ){ modStat = ( modStat || event->isModded( deep ) ); }
            return modStat;
        }
        void lheNode::setInit( std::shared_ptr<initNode> initNod ){ init = initNod; }
        void lheNode::setHeader( std::shared_ptr<lheHead> headNod ){ header = headNod; }
        void lheNode::addWgt( size_t index, newWgt& addedWgt ){
            header->addWgt( index, addedWgt.getHeadWgt() );
            auto wgtsVec = addedWgt.getBodyWgts();
            for( size_t k = 0 ; k < wgtsVec.size() ; ++k ){
                events[k]->addWgt( wgtsVec[k] );
            }
        }
        void lheNode::addWgt( size_t index, newWgt& addedWgt, std::string& idTag ){
            header->addWgt( index, addedWgt.getHeadWgt(), idTag );
            auto wgtsVec = addedWgt.getBodyWgts();
            for( size_t k = 0 ; k < wgtsVec.size() ; ++k ){
                events[k]->addWgt( wgtsVec[k] );
            }
        }
        void lheNode::setRelStats( std::vector<std::string_view>& particles ){
            relStat = particles;
        }
        std::vector<std::string_view>& lheNode::getRelStats(){
            return relStat;
        }
        void lheNode::setSameSort( sortFcn& sortF ){
            particleSort = sortF;
        }
        sortFcn& lheNode::getSameSort(){
            return particleSort;
        }
        void lheNode::setStatSort( statSort& statS ){
            statParticleSort = statS;
        }
        statSort& lheNode::getStatSort(){
            return statParticleSort;
        }
        void lheNode::headerWriter(){
            nodeContent += "\n" + *header->nodeWriter();
        }
        void lheNode::initWriter(){
            nodeContent += *init->nodeWriter();
        }
        void lheNode::eventWriter(){
            for( auto event : events ){
                nodeContent += *event->nodeWriter();
            }
        }
        void lheNode::contWriter(){
            nodeContent = "";
            headerWriter();
            initWriter();
            eventWriter();
        }
        void lheNode::fullWriter(){
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
        std::shared_ptr<std::string> lheNode::nodeWriter() {
            if( isModded( true ) || !isWritten() ){ fullWriter(); }
            return writtenSelf;
        }

    // ZW: function for extracting event information from
    // LHE files
    std::vector<std::shared_ptr<std::vector<double>>> valExtraction( lheNode& lheFile )
    {
        bool getGs = true;
        auto momVec = std::make_shared<std::vector<double>>();
        auto wgtVec = std::make_shared<std::vector<double>>();
        auto gVec = std::make_shared<std::vector<double>>();
        auto events = lheFile.getEvents();
        momVec->reserve( events.size() * 4 * std::stoi(std::string(events[0]->getHead().getNprt())) );
        wgtVec->reserve( events.size() );
        gVec->reserve( events.size() );
        if( getGs ){
        for( auto event : events )
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
        for( auto event : events )
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
    // and return a REX format event object
    std::shared_ptr<event> evPtrParsor( std::string_view parseFile, size_t& initPos, size_t& endPos )
    {
        auto currNode = std::make_shared<event>(parseFile, initPos);
        initPos = nodeStartFind( parseFile, initPos + 1 );
        while( initPos < endPos )
        {
            currNode->addChild(xmlPtrParser( parseFile, initPos, endPos ));
        }
        size_t equalSign = parseFile.find_first_of("=>", initPos);
        size_t nodeInitEnd = parseFile.find(">", initPos);
        while( equalSign < nodeInitEnd ){
            currNode->addTag( xmlTagParser(parseFile, equalSign) );
        }
        initPos = nodeStartFind( parseFile, endPos );
        endPos = nodeEndFind( parseFile, endPos + 1 );
        return currNode;
    }

    // ZW: fcn for parsing an LHE format header
    // and return a REX format lheHead object
    std::shared_ptr<lheHead> lheHeadParser( std::string_view parseFile, size_t& initPos, size_t& endPos )
    {
        auto currNode = std::make_shared<lheHead>(parseFile, initPos);
        initPos = nodeStartFind( parseFile, initPos + 1 );
        while( initPos < endPos )
        {
            currNode->addChild(xmlPtrParser( parseFile, initPos, endPos ));
            if( currNode->getChildren()[ currNode->getChildren().size() - 1 ]->getName() == "init" ){ continue; }
            if( currNode->getChildren()[ currNode->getChildren().size() - 1 ]->getName() == "slha" ){
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
        initPos = nodeStartFind( parseFile, endPos );
        endPos = nodeEndFind( parseFile, endPos + 1 );
        return currNode;
    }

    // ZW: fcn for parsing an LHE format file
    // and return a REX format LHE node object
    std::shared_ptr<lheNode> lheParser( std::string_view parseFile, size_t& initPos, size_t& endPos )
    {
        auto currNode = std::make_shared<lheNode>(parseFile, initPos);
        initPos = nodeStartFind( parseFile, initPos + 1 );
        while( initPos < endPos )
        {
            if( parseFile.substr( initPos, 6 ) == "<event" ){
                currNode->getEvents().push_back( evPtrParsor( parseFile, initPos, endPos ) );
                continue;
            } else if( parseFile.substr( initPos, 7 ) == "<header"  ){
                currNode->setHeader(lheHeadParser( parseFile, initPos, endPos ));
                continue;
            } else if( parseFile.substr( initPos, 5 ) == "<init"  ){
                currNode->setInit( std::make_shared<initNode>( parseFile, initPos ) );
                initPos = nodeStartFind( parseFile, endPos );
                endPos = nodeEndFind( parseFile, nodeEndFind( parseFile, endPos + 1 ) + 1);
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
        initPos = nodeStartFind( parseFile, endPos );
        endPos = nodeEndFind( parseFile, endPos + 1 );
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
        std::vector<size_t> orderMOne;
        std::vector<size_t> orderOne;
        std::vector<size_t> orderMTwo;
        std::vector<size_t> orderTwo;
        std::vector<size_t> orderThree;
        std::vector<size_t> orderNine;
        std::map<std::string_view,std::vector<std::string_view>> valVecs{{"-1", minusOne}, {"1", plusOne}, {"-2", minusTwo}, {"2", plusTwo}, {"3", plusThree}, {"-9", minusNine}};
        std::map<std::string_view,std::vector<size_t>> orderVecs{{"-1", orderMOne}, {"1", orderOne}, {"-2", orderMTwo}, {"2", orderTwo}, {"3", orderThree}, {"9",orderNine}};
        lheProc( event& eventNode )
        {
            for( auto prt : eventNode.getPrts() )
            {
                valVecs[prt->getStatus()].push_back(prt->getPDG());
            }
            for( auto valVec = valVecs.begin() ; valVec!= valVecs.end() ; ++valVec ){
                if( valVec->second.size() == 0 ){ continue; }
                orderVecs[valVec->first] = *stoiSort( valVec->second );
            }
        }
        std::shared_ptr<std::string> writer(){
            auto written = std::make_shared<std::string>();
            for( auto inits : valVecs["-1"] ){
                written->append(inits);
                written->append(" ");
            }
            if( valVecs["2"].size() > 0 ){
                written->append("> ");
                for( auto inits : valVecs["2"] ){
                    written->append(inits);
                    written->append(" ");
                }
            }
            written->append("> ");
            for( auto inits : valVecs["1"] ){
                written->append(inits);
                written->append(" ");
            }
            return written;
        }
    };

    // ZW: fcn for uploading text files to the program
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

    // ZW: fcn for extracting the full
    // process information from an LHE event
    std::shared_ptr<std::map<std::string_view, std::vector<std::string_view>>> pdgXtract( event& currEv )
    {
        auto currProc = std::make_shared<std::map<std::string_view, std::vector<std::string_view>>>();
        auto &useProc = *currProc;
        for( auto prt : currEv.getPrts() )
        {
            useProc[ prt->getStatus() ].push_back(prt->getPDG());
        }
        return currProc;
    }
    
    template <typename T>
    bool chaoticVecComp( const std::vector<T>& vec1, const std::vector<size_t> order1, const std::vector<T>& vec2, const std::vector<size_t> order2 )
    {
        if( vec1.size()!= vec2.size() ){ return false; }
        for( size_t i = 0; i < vec1.size(); i++ ){
            if( vec1[order1[i]]!= vec2[order2[i]] ){ return false; }
        }
        return true;
    }

    // ZW: fcn for comparing two processes in the
    // format output by pdgXtract
    bool sameProcString( std::map<std::string_view, std::vector<std::string_view>>& firstVec, std::map<std::string_view,
     std::vector<std::string_view>>& secVec, const std::vector<std::string_view>& statVec )
    {
        if( firstVec.size() != secVec.size() ){return false;}
        for(auto code : statVec )
        {
            if( firstVec[code] != secVec[code] ){ return false; }
        }
        return true;
    }

    bool sameProcString( std::map<std::string_view, std::vector<std::string_view>>& firstVec, std::map<std::string_view, std::vector<size_t>>& firstOrder,
     std::map<std::string_view, std::vector<std::string_view>>& secVec, std::map<std::string_view, std::vector<size_t>>& secondOrder,
     std::vector<std::string_view>& statVec )
    {
        if( firstVec.size() != secVec.size() ){return false;}
        for(auto code : statVec )
        {
            if( !chaoticVecComp(firstVec[code], firstOrder[code], secVec[code], secondOrder[code]) ){ return false; }
        }
        return true;
    }

    // ZW: fcn for processes in the lheProc struct format
    bool procComp( lheProc& firstProc, lheProc& secProc,  std::vector<std::string_view> statVec )
    {
        for( auto stat : statVec )
        {
            if( firstProc.valVecs.at(stat).size() != secProc.valVecs.at(stat).size() ){ return false; }
            if( !chaoticVecComp( firstProc.valVecs[stat], firstProc.orderVecs[stat], secProc.valVecs[stat], secProc.orderVecs[stat] ) ){ return false; }
        }
        return true;
    }

    bool evProcComp( event& firstEv, event& secEv, std::vector<std::string_view> statVec = {"-1", "1"} )
    {
        for( auto stat : statVec )
        {
            if( firstEv.getProc()[stat].size()!= secEv.getProc()[stat].size() ){ return false; }
            if(!chaoticVecComp( firstEv.getProc()[stat], firstEv.getProcOrder()[stat], 
                secEv.getProc()[stat], secEv.getProcOrder()[stat] ) ){ return false; }
        }
        return true;
    }

    bool evProcComp( event& firstEv, event& secEv, std::vector<std::string_view> statVec, 
    sortFcn sorter )
    {
        for( auto stat : statVec )
        {
            if( firstEv.getProc(sorter)[stat].size()!= secEv.getProc(sorter)[stat].size() ){ return false; }
            if(!chaoticVecComp( firstEv.getProc(sorter)[stat], firstEv.getProcOrder(sorter)[stat], 
                secEv.getProc(sorter)[stat], secEv.getProcOrder(sorter)[stat] ) ){ return false; }
        }
        return true;
    }

    bool evProcComp( event& firstEv, event& secEv, std::vector<std::string_view> statVec, 
    statSort sorter )
    {
        for( auto stat : statVec )
        {
            if( firstEv.getProc(sorter)[stat].size()!= secEv.getProc(sorter)[stat].size() ){ return false; }
            if(!chaoticVecComp( firstEv.getProc(sorter)[stat], firstEv.getProcOrder(sorter)[stat], 
                secEv.getProc(sorter)[stat], secEv.getProcOrder(sorter)[stat] ) ){ return false; }
        }
        return true;
    }

    bool evProcComp( const event& firstEv, const event& secEv, std::vector<std::string_view> statVec = {"-1", "1"} )
    {
        for( auto stat : statVec )
        {
            if( firstEv.getProc().at(stat).size()!= secEv.getProc().at(stat).size() ){ return false; }
            if(!chaoticVecComp( firstEv.getProc().at(stat), firstEv.getProcOrder().at(stat), 
                secEv.getProc().at(stat), secEv.getProcOrder().at(stat) ) ){ return false; }
        }
        return true;
    }

    bool evProcComp( const event& firstEv, const event& secEv, std::vector<std::string_view> statVec, 
    sortFcn sorter )
    {
        for( auto stat : statVec )
        {
            if( firstEv.getProc().at(stat).size()!= secEv.getProc().at(stat).size() ){ return false; }
            if(!chaoticVecComp( firstEv.getProc().at(stat), firstEv.getProcOrder().at(stat), 
                secEv.getProc().at(stat), secEv.getProcOrder().at(stat) ) ){ return false; }
        }
        return true;
    }

    bool evProcComp( const event& firstEv, const event& secEv, std::vector<std::string_view> statVec, 
    statSort sorter )
    {
        for( auto stat : statVec )
        {
            if( firstEv.getProc().at(stat).size()!= secEv.getProc().at(stat).size() ){ return false; }
            if(!chaoticVecComp( firstEv.getProc().at(stat), firstEv.getProcOrder().at(stat), 
                secEv.getProc().at(stat), secEv.getProcOrder().at(stat) ) ){ return false; }
        }
        return true;
    }

        bool eventComp::operator()( event& firstEv, event& secEv){
            if( firstEv.isSpecSort() ) {return evProcComp( firstEv, secEv, {"-1", "1"}, firstEv.getStatSort());}
            else {return evProcComp( firstEv, secEv, {"-1", "1"}, firstEv.getSortFcn() );}
        }
        bool eventComp::operator()( const event& firstEv, const event& secEv) const {
            if( firstEv.isSpecSort() ) {return evProcComp( firstEv, secEv, {"-1", "1"}, firstEv.getStatSort());}
            else {return evProcComp( firstEv, secEv, {"-1", "1"}, firstEv.getSortFcn() );}
        }
        bool eventComp::operator()(event& firstEv, event& secEv, std::vector<std::string_view> statVec){
            if( firstEv.isSpecSort() ) {return evProcComp( firstEv, secEv, statVec, firstEv.getStatSort());}
            else {return evProcComp( firstEv, secEv, statVec, firstEv.getSortFcn() );}
        }

    // ZW: fcn for checking whether a list of pdgXtract format
    // processes sourceProcList contains a given process newProc
    bool procVecContains( std::vector<std::shared_ptr<std::map<std::string_view, std::vector<std::string_view>>>>& sourceProcList, 
    std::map<std::string_view, std::vector<std::string_view>>& newProc, const std::vector<std::string_view>& statVec  )
    {\
        for( auto proc : sourceProcList )
        {
            if( sameProcString( *proc, newProc, statVec ) ){ return true; }
        }
        return false;
    }

    // ZW: fcn for checking whether a vector of lheProc structs
    // procList contains a given lheProc nuProc
    bool procListComp( const std::vector<std::shared_ptr<lheProc>>& procList, lheProc& nuProc, std::vector<std::string_view> statVec )
    {
        if( procList.size() != 0 ){
            for(auto proc : procList )
            {
                if( procComp( *proc, nuProc, statVec ) ){ return true; }
            }
        }
        return false;
    }

    bool evProcListComp( std::vector<std::shared_ptr<event>>& procList, event& nuEv, std::vector<std::string_view> statVec )
    {
        if( procList.size()!= 0 ){
            for( auto ev : procList )
            {
                if( evProcComp( *ev, nuEv, statVec ) ){ return true; }
            }
        }
        return false;
    }

    bool evProcListComp( std::vector<std::shared_ptr<event>>& procList, event& nuEv, std::vector<std::string_view> statVec, 
    sortFcn sorter )
    {
        if( procList.size()!= 0 ){
            for( auto ev : procList )
            {
                if( evProcComp( *ev, nuEv, statVec, sorter ) ){ return true; }
            }
        }
        return false;
    }

    bool evProcListComp( std::vector<std::shared_ptr<event>>& procList, event& nuEv, std::vector<std::string_view> statVec, 
    statSort sorter )
    {
        if( procList.size()!= 0 ){
            for( auto ev : procList )
            {
                if( evProcComp( *ev, nuEv, statVec, sorter ) ){ return true; }
            }
        }
        return false;
    }

    // ZW: fcn for extracting the different processes
    // in a given REX format LHE file in the pdgXtract format
    std::vector<std::shared_ptr<std::map<std::string_view, std::vector<std::string_view>>>> procExtractor( lheNode& lheFile )
    {
        std::vector<std::shared_ptr<std::map<std::string_view, std::vector<std::string_view>>>> procList;
        const static std::vector<std::string_view> statVec = { "-1", "1", "-2", "2", "3", "-9" };
        for( auto event : lheFile.getEvents() )
        {
            auto currProc = pdgXtract( *event );
            if( procVecContains( procList, *currProc, statVec ) ){ continue; }
            procList.push_back(currProc);
        }
        return procList;
    }

    // ZW: fcn for extracting the different processes
    // in a given REX format LHE file in the lheProc format
    std::vector<std::shared_ptr<lheProc>> processPull( lheNode& lheFile, 
    std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        //const static std::vector<std::string_view> statVec = { "-1", "1", "-2", "2", "3", "-9" };
        std::vector<std::shared_ptr<lheProc>> procsList{};
        for( auto event : lheFile.getEvents() )
        {
            auto currProc =  std::make_shared<lheProc>( *event );
            if( procListComp( procsList, *currProc, statVec ) ){ continue; }
            procsList.push_back( currProc );
        }
        return procsList;
    }

    std::vector<std::shared_ptr<event>> evProcessPull( lheNode& lheFile, std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        //const static std::vector<std::string_view> statVec = { "-1", "1", "-2", "2", "3", "-9" };
        std::vector<std::shared_ptr<event>> procsList{};
        for( auto currEv : lheFile.getEvents() )
        {
            if( evProcListComp( procsList, *currEv, statVec ) ){ continue; }
            procsList.push_back( currEv );
        }
        return procsList;
    }

    std::vector<std::shared_ptr<event>> evProcessPull( lheNode& lheFile, 
    sortFcn sorter,
    std::vector<std::string_view> statVec = { "-1", "1" })
    {
        //const static std::vector<std::string_view> statVec = { "-1", "1", "-2", "2", "3", "-9" };
        std::vector<std::shared_ptr<event>> procsList{};
        lheFile.setSameSort(sorter);
        for( auto currEv : lheFile.getEvents() )
        {
            if( evProcListComp( procsList, *currEv, statVec, sorter ) ){ continue; }
            procsList.push_back( currEv );
        }
        return procsList;
    }

    std::vector<std::shared_ptr<event>> evProcessPull( lheNode& lheFile, 
    statSort sorter,
    std::vector<std::string_view> statVec = { "-1", "1" })
    {
        //const static std::vector<std::string_view> statVec = { "-1", "1", "-2", "2", "3", "-9" };
        std::vector<std::shared_ptr<event>> procsList{};
        lheFile.setStatSort(sorter);
        for( auto currEv : lheFile.getEvents() )
        {
            if( evProcListComp( procsList, *currEv, statVec, sorter ) ){ continue; }
            procsList.push_back( currEv );
        }
        return procsList;
    }

    // ZW: fcn for keeping track of subprocess ordering
    // in LHE file
    size_t procPos( const std::vector<std::shared_ptr<lheProc>>& evtSet, lheProc& currProc, 
        std::vector<std::string_view>& statVec )
    {
        for( size_t k = 0 ; k < evtSet.size() ; ++k )
        {
            for( auto stat : statVec )
            {
                if( evtSet[k]->valVecs[stat] != currProc.valVecs[stat] ){ break; }   
            }
            return k;
        }
        return evtSet.size();
    }

    size_t evProcPos( const std::vector<std::shared_ptr<event>>& evtSet, event& currEv, 
        std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        for( size_t k = 0 ; k < evtSet.size() ; ++k )
        { 
            if( evProcComp(*evtSet[k], currEv, statVec) ){ return k; }
        }
        return evtSet.size();
    }

    size_t evProcPos( const std::vector<std::shared_ptr<event>>& evtSet, event& currEv,  
    sortFcn sorter, std::vector<std::string_view> statVec = {"-1", "1"} )
    {
        for( size_t k = 0 ; k < evtSet.size() ; ++k )
        { 
            if( evProcComp(*evtSet[k], currEv, statVec, sorter) ){ return k; }
        }
        return evtSet.size();
    }

    size_t evProcPos( const std::vector<std::shared_ptr<event>>& evtSet, event& currEv, 
    statSort sorter, std::vector<std::string_view> statVec = {"-1", "1"} )
    {
        for( size_t k = 0 ; k < evtSet.size() ; ++k )
        { 
            if( evProcComp(*evtSet[k], currEv, statVec, sorter) ){ return k; }
        }
        return evtSet.size();
    }

    // ZW: fcn for extracting the subprocess ordering
    // of LHE file
    std::vector<std::shared_ptr<std::vector<bool>>> procOrder( lheNode& lheFile, const std::vector<std::shared_ptr<lheProc>>& evtSet,
        std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        //const static std::vector<std::string> statVec = { "-1", "1", "-2", "2", "3", "-9" };
        std::vector<std::shared_ptr<std::vector<bool>>> eventBools( evtSet.size(), std::make_shared<std::vector<bool>> ( lheFile.getEvents().size() ));
        //std::vector<std::vector<bool>> pracBools( evtSet.size(), std::vector<bool> ( lheFile.getEvents().size() ));
        for( auto boolSets : eventBools ){
            std::fill( boolSets->begin(), boolSets->end(), false );
        }
        for( size_t k = 0 ; k < lheFile.getEvents().size() ; ++k )
        {
            auto currProc = lheProc(*lheFile.getEvents()[k]);
            eventBools[ procPos(evtSet, currProc, statVec) ]->at( k ) = true;
        }
        //for( size_t k = 0 ; k < eventBools.size() ; ++k )
        //{
        //    eventBools[k] = std::make_shared<std::vector<bool>>( pracBools[k] );
        //}
        return eventBools; 
    }

    std::vector<std::shared_ptr<std::vector<bool>>> evProcOrder( lheNode& lheFile, const std::vector<std::shared_ptr<event>>& evtSet,
        std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        std::vector<std::shared_ptr<std::vector<bool>>> eventBools;
        eventBools.reserve(evtSet.size());
        for (size_t i = 0; i < evtSet.size(); ++i) {
            eventBools.push_back(std::make_shared<std::vector<bool>>(lheFile.getEvents().size(), false));
        }
        for( size_t k = 0 ; k < lheFile.getEvents().size() ; ++k )
        {
            eventBools[ evProcPos(evtSet, *lheFile.getEvents()[k], statVec) ]->at( k ) = true;
        }
        return eventBools;
    }

    std::vector<std::shared_ptr<std::vector<bool>>> evProcOrder( lheNode& lheFile, const std::vector<std::shared_ptr<event>>& evtSet, 
    sortFcn sorter,
        std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        std::vector<std::shared_ptr<std::vector<bool>>> eventBools;
        eventBools.reserve(evtSet.size());
        for (size_t i = 0; i < evtSet.size(); ++i) {
            eventBools.push_back(std::make_shared<std::vector<bool>>(lheFile.getEvents().size(), false));
        }
        for( size_t k = 0 ; k < lheFile.getEvents().size() ; ++k )
        {
            eventBools[ evProcPos(evtSet, *lheFile.getEvents()[k], sorter, statVec) ]->at( k ) = true;
        }
        return eventBools;
    }

    std::vector<std::shared_ptr<std::vector<bool>>> evProcOrder( lheNode& lheFile, const std::vector<std::shared_ptr<event>>& evtSet, 
    statSort sorter,
        std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        std::vector<std::shared_ptr<std::vector<bool>>> eventBools;
        eventBools.reserve(evtSet.size());
        for (size_t i = 0; i < evtSet.size(); ++i) {
            eventBools.push_back(std::make_shared<std::vector<bool>>(lheFile.getEvents().size(), false));
        }
        for( size_t k = 0 ; k < lheFile.getEvents().size() ; ++k )
        {
            eventBools[ evProcPos(evtSet, *lheFile.getEvents()[k], sorter, statVec) ]->at( k ) = true;
        }
        return eventBools;
    }

    // ZW: fcn for reordering LHE file based on subprocess
    std::shared_ptr<std::vector<std::shared_ptr<event>>> eventReOrder( lheNode& lheFile, std::vector<bool> relProc )
    {
        auto reOrdered = std::make_shared<std::vector<std::shared_ptr<event>>>();
        reOrdered->reserve( std::count( relProc.begin(), relProc.end(), true ) );
        for( size_t k = 0 ; k < relProc.size() ; ++k )
        {
            if(!relProc[k]){continue;}
            reOrdered->push_back( lheFile.getEvents()[k] );
        }
        return reOrdered;
    }

    // ZW: wrapper for eventReOrder
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> lheReOrder( lheNode& lheFile,
        std::vector<std::string_view> statVec = { "-1", "1" }  )
    {
        auto procSets = processPull( lheFile, statVec ); 
        auto relProcs = procOrder( lheFile, procSets, statVec ); 
        std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> ordProcs(procSets.size());
        for( size_t k = 0 ; k < relProcs.size() ; ++k )
        { 
            ordProcs[k] = eventReOrder( lheFile, *relProcs[k] );
        }
        return ordProcs;
    }

    std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> lheEvReOrder( lheNode& lheFile,
        std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        auto procSets = evProcessPull( lheFile, statVec ); 
        auto relProcs = evProcOrder( lheFile, procSets, statVec ); 
        std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> ordProcs(procSets.size());
        for( size_t k = 0 ; k < relProcs.size() ; ++k )
        { 
            ordProcs[k] = eventReOrder( lheFile, *relProcs[k] );
        }
        return ordProcs;
    }

    std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> lheEvReOrder( lheNode& lheFile,
        std::vector<std::shared_ptr<event>> procSets, std::vector<std::shared_ptr<std::vector<bool>>> relProcs,
        std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        //auto procSets = evProcessPull( lheFile, statVec ); 
        //auto relProcs = evProcOrder( lheFile, procSets, statVec ); 
        std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> ordProcs(procSets.size());
        for( size_t k = 0 ; k < relProcs.size() ; ++k )
        { 
            ordProcs[k] = eventReOrder( lheFile, *relProcs[k] );
        }
        return ordProcs;
    }

    std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> lheEvReOrder( lheNode& lheFile, 
    sortFcn sorter,
    std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        auto procSets = evProcessPull( lheFile, sorter, statVec ); 
        auto relProcs = evProcOrder( lheFile, procSets, sorter, statVec ); 
        std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> ordProcs(procSets.size());
        for( size_t k = 0 ; k < relProcs.size() ; ++k )
        { 
            ordProcs[k] = eventReOrder( lheFile, *relProcs[k] );
        }
        return ordProcs;
    }

    std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> lheEvReOrder( lheNode& lheFile, 
    std::vector<std::shared_ptr<event>> procSets, std::vector<std::shared_ptr<std::vector<bool>>> relProcs,
    sortFcn sorter, std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        //auto procSets = evProcessPull( lheFile, sorter, statVec ); 
        //auto relProcs = evProcOrder( lheFile, procSets, sorter, statVec ); 
        std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> ordProcs(procSets.size());
        for( size_t k = 0 ; k < relProcs.size() ; ++k )
        { 
            ordProcs[k] = eventReOrder( lheFile, *relProcs[k] );
        }
        return ordProcs;
    }

    std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> lheEvReOrder( lheNode& lheFile, 
    statSort sorter,
    std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        auto procSets = evProcessPull( lheFile, sorter, statVec ); 
        auto relProcs = evProcOrder( lheFile, procSets, sorter, statVec ); 
        std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> ordProcs(procSets.size());
        for( size_t k = 0 ; k < relProcs.size() ; ++k )
        { 
            ordProcs[k] = eventReOrder( lheFile, *relProcs[k] );
        }
        return ordProcs;
    }

    std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> lheEvReOrder( lheNode& lheFile, 
    std::vector<std::shared_ptr<event>> procSets, std::vector<std::shared_ptr<std::vector<bool>>> relProcs,
    statSort sorter, std::vector<std::string_view> statVec = { "-1", "1" } )
    {
        //auto procSets = evProcessPull( lheFile, sorter, statVec ); 
        //auto relProcs = evProcOrder( lheFile, procSets, sorter, statVec ); 
        std::vector<std::shared_ptr<std::vector<std::shared_ptr<event>>>> ordProcs(procSets.size());
        for( size_t k = 0 ; k < relProcs.size() ; ++k )
        { 
            ordProcs[k] = eventReOrder( lheFile, *relProcs[k] );
        }
        return ordProcs;
    }

    // ZW: transposed event information struct
        evtInfo::evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile ){
            int nEvt = lheFile.size();
            wgts.reserve(nEvt); scales.reserve(nEvt); aQEDs.reserve(nEvt); aQCDs.reserve(nEvt); nprts.reserve(nEvt); procIDs.reserve(nEvt);
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
        evtInfo::evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const std::vector<std::string_view>& statVec ){
            int nEvt = lheFile.size();
            wgts.reserve(nEvt); scales.reserve(nEvt); aQEDs.reserve(nEvt); aQCDs.reserve(nEvt); relNPrts.reserve(nEvt); procIDs.reserve(nEvt);
            for( auto evt : lheFile )
            {
                wgts.push_back(evt->getHead().getWeight());
                scales.push_back(evt->getHead().getScale());
                aQEDs.push_back(evt->getHead().getAQED());
                aQCDs.push_back(evt->getHead().getAQCD());
                size_t nPrt = 0;
                for( auto stat : statVec ){ nPrt += evt->getProc()[stat].size(); }
                relNPrts.push_back(nPrt);
                procIDs.push_back(evt->getHead().getProcID());
            }
        }
        evtInfo::evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const std::vector<std::string_view>& statVec, 
        sortFcn sorter ){
            int nEvt = lheFile.size();
            wgts.reserve(nEvt); scales.reserve(nEvt); aQEDs.reserve(nEvt); aQCDs.reserve(nEvt); relNPrts.reserve(nEvt); procIDs.reserve(nEvt);
            for( auto evt : lheFile )
            {
                wgts.push_back(evt->getHead().getWeight());
                scales.push_back(evt->getHead().getScale());
                aQEDs.push_back(evt->getHead().getAQED());
                aQCDs.push_back(evt->getHead().getAQCD());
                size_t nPrt = 0;
                for( auto stat : statVec ){ nPrt += evt->getProc(sorter)[stat].size(); }
                relNPrts.push_back(nPrt);
                procIDs.push_back(evt->getHead().getProcID());
            }
        }
        evtInfo::evtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const std::vector<std::string_view>& statVec, 
        statSort sorter ){
            int nEvt = lheFile.size();
            wgts.reserve(nEvt); scales.reserve(nEvt); aQEDs.reserve(nEvt); aQCDs.reserve(nEvt); relNPrts.reserve(nEvt); procIDs.reserve(nEvt);
            for( auto evt : lheFile )
            {
                wgts.push_back(evt->getHead().getWeight());
                scales.push_back(evt->getHead().getScale());
                aQEDs.push_back(evt->getHead().getAQED());
                aQCDs.push_back(evt->getHead().getAQCD());
                size_t nPrt = 0;
                for( auto stat : statVec ){ nPrt += evt->getProc(sorter)[stat].size(); }
                relNPrts.push_back(nPrt);
                procIDs.push_back(evt->getHead().getProcID());
            }
        }

    // ZW: transposed particle information struct
        prtInfo::prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt ){
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
                    for( size_t k = 0 ; k < 2 ; ++k )
                    {
                        moms.push_back( prt->getMom()[k] );
                        mothers.push_back( prt->getMothers()[k] );
                        icols.push_back( prt->getColor()[k] );
                    }
                    moms.push_back( prt->getMom()[2] );
                }
            }
        }
        prtInfo::prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, const std::vector<std::string_view>& statVec ){
            int nEvt = lheFile.size(); 
            moms.reserve(4*nPrt*nEvt); vtims.reserve(nPrt*nEvt); masses.reserve(nPrt*nEvt); pdgs.reserve(nPrt*nEvt); 
            spins.reserve(nPrt*nEvt); statuses.reserve(nPrt*nEvt); mothers.reserve(2*nPrt*nEvt); icols.reserve(2*nPrt*nEvt);
            for( auto evt : lheFile )
            {
                for( auto stat : statVec )
                {
                    for( auto i : evt->getProcOrder()[stat] )
                    {
                        auto prt = evt->getPrts()[i];
                        moms.push_back( prt->getE() );
                        masses.push_back( prt->getMass() );
                        vtims.push_back( prt->getVTim() );
                        spins.push_back( prt->getSpin() );
                        statuses.push_back( prt->getStatus() );
                        pdgs.push_back( prt->getPDG() );
                        for( size_t k = 0 ; k < 2 ; ++k )
                        {
                            moms.push_back( prt->getMom()[k] );
                            mothers.push_back( prt->getMothers()[k] );
                            icols.push_back( prt->getColor()[k] );
                        }
                        moms.push_back( prt->getMom()[2] );
                    }
                }
            }
        }
        prtInfo::prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, const std::vector<std::string_view>& statVec, 
        sortFcn sorter ){
            int nEvt = lheFile.size(); 
            moms.reserve(4*nPrt*nEvt); vtims.reserve(nPrt*nEvt); masses.reserve(nPrt*nEvt); pdgs.reserve(nPrt*nEvt); 
            spins.reserve(nPrt*nEvt); statuses.reserve(nPrt*nEvt); mothers.reserve(2*nPrt*nEvt); icols.reserve(2*nPrt*nEvt);
            for( auto evt : lheFile )
            {
                for( auto stat : statVec )
                {
                    for( auto i : evt->getProcOrder(sorter)[stat] )
                    {
                        auto prt = evt->getPrts()[i];
                        moms.push_back( prt->getE() );
                        masses.push_back( prt->getMass() );
                        vtims.push_back( prt->getVTim() );
                        spins.push_back( prt->getSpin() );
                        statuses.push_back( prt->getStatus() );
                        pdgs.push_back( prt->getPDG() );
                        for( size_t k = 0 ; k < 2 ; ++k )
                        {
                            moms.push_back( prt->getMom()[k] );
                            mothers.push_back( prt->getMothers()[k] );
                            icols.push_back( prt->getColor()[k] );
                        }
                        moms.push_back( prt->getMom()[2] );
                    }
                }
            }
        }
        prtInfo::prtInfo( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, const std::vector<std::string_view>& statVec, 
        statSort sorter ){
            int nEvt = lheFile.size(); 
            moms.reserve(4*nPrt*nEvt); vtims.reserve(nPrt*nEvt); masses.reserve(nPrt*nEvt); pdgs.reserve(nPrt*nEvt); 
            spins.reserve(nPrt*nEvt); statuses.reserve(nPrt*nEvt); mothers.reserve(2*nPrt*nEvt); icols.reserve(2*nPrt*nEvt);
            for( auto evt : lheFile )
            {
                for( auto stat : statVec )
                {
                    for( auto i : evt->getProcOrder(sorter)[stat] )
                    {
                        auto prt = evt->getPrts()[i];
                        moms.push_back( prt->getE() );
                        masses.push_back( prt->getMass() );
                        vtims.push_back( prt->getVTim() );
                        spins.push_back( prt->getSpin() );
                        statuses.push_back( prt->getStatus() );
                        pdgs.push_back( prt->getPDG() );
                        for( size_t k = 0 ; k < 2 ; ++k )
                        {
                            moms.push_back( prt->getMom()[k] );
                            mothers.push_back( prt->getMothers()[k] );
                            icols.push_back( prt->getColor()[k] );
                        }
                        moms.push_back( prt->getMom()[2] );
                    }
                }
            }
        }

    // ZW: transposed LHE file with a single process type
        transMonoLHE::transMonoLHE( const std::vector<std::shared_ptr<REX::event>>& lheFile , const int nPrt ){
            evtsHead = evtInfo(lheFile);
            evtsData = prtInfo(lheFile, nPrt);
            process = lheFile[0];
        }
        transMonoLHE::transMonoLHE( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, const std::vector<std::string_view>& statVec ){
            evtsHead = evtInfo(lheFile, statVec);
            evtsData = prtInfo(lheFile, nPrt, statVec);
            process = lheFile[0];
        }
        transMonoLHE::transMonoLHE( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, 
        sortFcn sorter,
        std::vector<std::string_view> statVec ){
            evtsHead = evtInfo(lheFile, statVec);
            evtsData = prtInfo(lheFile, nPrt, statVec, sorter);
            process = lheFile[0];
        }
        transMonoLHE::transMonoLHE( const std::vector<std::shared_ptr<REX::event>>& lheFile, const int nPrt, 
        statSort sorter,
        std::vector<std::string_view> statVec){
            evtsHead = evtInfo(lheFile, statVec);
            evtsData = prtInfo(lheFile, nPrt, statVec, sorter);
            process = lheFile[0];
        }

    // ZW: transposed LHE file ordered by subprocess
        transLHE::transLHE(){ return; }
        transLHE::transLHE( lheNode& lheFile )
        {
            procSets = evProcessPull( lheFile ); 
            relProcs = evProcOrder( lheFile, procSets );
            xmlFile = lheFile.getFile();
            auto procsOrdered = lheEvReOrder( lheFile, procSets, relProcs ); 
            subProcs = std::vector<std::shared_ptr<transMonoLHE>>( procsOrdered.size() ); 
            for( size_t k = 0 ; k < procsOrdered.size() ; ++k )
            { 
                subProcs[k] = std::make_shared<transMonoLHE>( *procsOrdered[k], procsOrdered[k]->at(0)->getNprt() );
            }
        }
        transLHE::transLHE( lheNode& lheFile, 
        sortFcn sorter, 
        const std::vector<std::string_view>& statVec  )
        {
            procSets = evProcessPull( lheFile, sorter, statVec ); 
            relProcs = evProcOrder( lheFile, procSets, sorter, statVec );
            xmlFile = lheFile.getFile();
            auto procsOrdered = lheEvReOrder( lheFile, procSets, relProcs, sorter, statVec ); 
            subProcs = std::vector<std::shared_ptr<transMonoLHE>>( procsOrdered.size() ); 
            for( size_t k = 0 ; k < procsOrdered.size() ; ++k )
            { 
                subProcs[k] = std::make_shared<transMonoLHE>( *procsOrdered[k], procsOrdered[k]->at(0)->getNprt(), sorter, statVec );
            }
        }
        transLHE::transLHE( lheNode& lheFile, 
        statSort sorter, 
        const std::vector<std::string_view>& statVec)
        {
            procSets = evProcessPull( lheFile, sorter, statVec ); 
            relProcs = evProcOrder( lheFile, procSets, sorter, statVec );
            xmlFile = lheFile.getFile();
            auto procsOrdered = lheEvReOrder( lheFile, procSets, relProcs, sorter, statVec ); 
            subProcs = std::vector<std::shared_ptr<transMonoLHE>>( procsOrdered.size() ); 
            for( size_t k = 0 ; k < procsOrdered.size() ; ++k )
            { 
                subProcs[k] = std::make_shared<transMonoLHE>( *procsOrdered[k], procsOrdered[k]->at(0)->getNprt(), sorter, statVec );
            }
        }
        transLHE::transLHE( lheNode& lheFile, const std::vector<std::string_view>& statVec )
        {
            procSets = evProcessPull( lheFile, statVec ); 
            relProcs = evProcOrder( lheFile, procSets, statVec );
            xmlFile = lheFile.getFile();
            auto procsOrdered = lheEvReOrder( lheFile, procSets, relProcs, statVec ); 
            subProcs = std::vector<std::shared_ptr<transMonoLHE>>( procsOrdered.size() ); 
            for( size_t k = 0 ; k < procsOrdered.size() ; ++k )
            { 
                subProcs[k] = std::make_shared<transMonoLHE>( *procsOrdered[k], procsOrdered[k]->at(0)->getNprt(), statVec );
            }
        }
//        template <typename T>
        std::shared_ptr<std::vector<double>> transLHE::vectorFlat( std::vector<std::shared_ptr<std::vector<double>>> vecVec )
        {
            if( vecVec.size() != relProcs.size() ) throw std::range_error("vectorFlat: input vector size does not match number of subprocesses");
            for( size_t k = 0 ; k < vecVec.size() ; ++k){
                if( vecVec[k]->size() == relProcs[k]->size() ) continue;
                else throw std::range_error("vectorFlat: input vector size does not match number of events for subprocess");
            }
            auto flatVec = std::make_shared<std::vector<double>>(relProcs[0]->size());
            for( size_t k = 0 ; k < relProcs.size() ; ++k ){
                size_t currInd = 0;
                for( size_t j = 0 ; j < relProcs[k]->size() ; ++j ){
                    if( relProcs[k]->at(j) ){
                        flatVec->at(currInd) = vecVec[k]->at(currInd);
                        ++currInd;
                    }
                }
            }
            return flatVec;
        }

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
        std::vector<bool> lheRetDs::getBools(){
            return { ebmup, xsecup, xerrup, xmaxup, xwgtup, scalup, aqedup, aqcdup, 
        pup, mass, vtimup, spinup };
        }

    // ZW: bool struct to define which int values
    // to extract transposed from LHE file
        std::vector<bool> lheRetInts::getBools(){
            return { idbmup, pdfgup, pdfsup, idwtup, nprup, lprup,
            nup, idprup, idup, istup, mothup, icolup };
        }

    // ZW: function for extracting transposed double values
    // from LHE file
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> lheValDoubles( lheNode& lheFile, lheRetDs vals )
    {
        // ZW: hard-setting returning g_S instead of a_S for now
        bool aStogS = true;
        auto boolVec = vals.getBools();
        const int noVals = std::count(boolVec.begin(), boolVec.end(), true);
        auto lheAOS = transLHE( lheFile );
        auto lheDos = std::make_shared<std::vector<std::shared_ptr<std::vector<double>>>>(noVals * lheAOS.subProcs.size() ); 
        std::vector<std::shared_ptr<std::vector<double>>> &lheDs = *lheDos;
        int currInd = 0;
        if( boolVec[0] ){ lheDs[currInd] = vecStoD( { lheFile.getInit()->getHead()->ebmup[0], lheFile.getInit()->getHead()->ebmup[1] } ); ++currInd; }
        if( boolVec[1] ){ 
            std::vector<std::string_view> xsecVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xsecVec.push_back(line->xsecup);
            }
            lheDs[currInd] = vecStoD( xsecVec );
             ++currInd; }
        if( boolVec[2] ){ 
            std::vector<std::string_view> xerrVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xerrVec.push_back(line->xerrup);
            }
            lheDs[currInd] = vecStoD( xerrVec );
             ++currInd; }
        if( boolVec[3] ){ 
            std::vector<std::string_view> xmaxVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xmaxVec.push_back(line->xmaxup);
            }
            lheDs[currInd] = vecStoD( xmaxVec );
             ++currInd; } 
        for( size_t k = 0 ; k < lheAOS.subProcs.size() ; ++k )
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

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> lheValDoubles(transLHE& lheAOS, lheRetDs vals )
    {
        // ZW: hard-setting returning g_S instead of a_S for now
        bool aStogS = true;
        auto boolVec = vals.getBools();
        const int noVals = std::count(boolVec.begin(), boolVec.end(), true);
        //auto lheAOS = transLHE( lheFile );
        auto lheDos = std::make_shared<std::vector<std::shared_ptr<std::vector<double>>>>(noVals * lheAOS.subProcs.size() ); 
        std::vector<std::shared_ptr<std::vector<double>>> &lheDs = *lheDos;
        int currInd = 0;
        for( size_t k = 0 ; k < lheAOS.subProcs.size() ; ++k )
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

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> lheValDoubles( lheNode& lheFile, 
    const std::vector<std::string_view>& statVec, lheRetDs vals = lheRetDs() )
    {
        // ZW: hard-setting returning g_S instead of a_S for now
        bool aStogS = true;
        auto boolVec = vals.getBools();
        const int noVals = std::count(boolVec.begin(), boolVec.end(), true);
        auto lheAOS = transLHE( lheFile, statVec );
        auto lheDos = std::make_shared<std::vector<std::shared_ptr<std::vector<double>>>>(noVals * lheAOS.subProcs.size() ); 
        std::vector<std::shared_ptr<std::vector<double>>> &lheDs = *lheDos;
        int currInd = 0;
        if( boolVec[0] ){ lheDs[currInd] = vecStoD( { lheFile.getInit()->getHead()->ebmup[0], lheFile.getInit()->getHead()->ebmup[1] } ); ++currInd; }
        if( boolVec[1] ){ 
            std::vector<std::string_view> xsecVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xsecVec.push_back(line->xsecup);
            }
            lheDs[currInd] = vecStoD( xsecVec );
             ++currInd; }
        if( boolVec[2] ){ 
            std::vector<std::string_view> xerrVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xerrVec.push_back(line->xerrup);
            }
            lheDs[currInd] = vecStoD( xerrVec );
             ++currInd; }
        if( boolVec[3] ){ 
            std::vector<std::string_view> xmaxVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xmaxVec.push_back(line->xmaxup);
            }
            lheDs[currInd] = vecStoD( xmaxVec );
             ++currInd; } 
        for( size_t k = 0 ; k < lheAOS.subProcs.size() ; ++k )
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

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> lheValDoubles( lheNode& lheFile, 
    sortFcn sorter, 
    const std::vector<std::string_view>& statVec = {"-1", "1"}, lheRetDs vals = lheRetDs() )
    {
        // ZW: hard-setting returning g_S instead of a_S for now
        bool aStogS = true;
        auto boolVec = vals.getBools();
        const int noVals = std::count(boolVec.begin(), boolVec.end(), true);
        auto lheAOS = transLHE( lheFile, sorter, statVec );
        auto lheDos = std::make_shared<std::vector<std::shared_ptr<std::vector<double>>>>(noVals * lheAOS.subProcs.size() ); 
        std::vector<std::shared_ptr<std::vector<double>>> &lheDs = *lheDos;
        int currInd = 0;
        if( boolVec[0] ){ lheDs[currInd] = vecStoD( { lheFile.getInit()->getHead()->ebmup[0], lheFile.getInit()->getHead()->ebmup[1] } ); ++currInd; }
        if( boolVec[1] ){ 
            std::vector<std::string_view> xsecVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xsecVec.push_back(line->xsecup);
            }
            lheDs[currInd] = vecStoD( xsecVec );
             ++currInd; }
        if( boolVec[2] ){ 
            std::vector<std::string_view> xerrVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xerrVec.push_back(line->xerrup);
            }
            lheDs[currInd] = vecStoD( xerrVec );
             ++currInd; }
        if( boolVec[3] ){ 
            std::vector<std::string_view> xmaxVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xmaxVec.push_back(line->xmaxup);
            }
            lheDs[currInd] = vecStoD( xmaxVec );
             ++currInd; } 
        for( size_t k = 0 ; k < lheAOS.subProcs.size() ; ++k )
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

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<double>>>> lheValDoubles( lheNode& lheFile, 
    statSort sorter, 
    const std::vector<std::string_view>& statVec = {"-1", "1"}, lheRetDs vals = lheRetDs() )
    {
        // ZW: hard-setting returning g_S instead of a_S for now
        bool aStogS = true;
        auto boolVec = vals.getBools();
        const int noVals = std::count(boolVec.begin(), boolVec.end(), true);
        auto lheAOS = transLHE( lheFile, sorter, statVec );
        auto lheDos = std::make_shared<std::vector<std::shared_ptr<std::vector<double>>>>(noVals * lheAOS.subProcs.size() ); 
        std::vector<std::shared_ptr<std::vector<double>>> &lheDs = *lheDos;
        int currInd = 0;
        if( boolVec[0] ){ lheDs[currInd] = vecStoD( { lheFile.getInit()->getHead()->ebmup[0], lheFile.getInit()->getHead()->ebmup[1] } ); ++currInd; }
        if( boolVec[1] ){ 
            std::vector<std::string_view> xsecVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xsecVec.push_back(line->xsecup);
            }
            lheDs[currInd] = vecStoD( xsecVec );
             ++currInd; }
        if( boolVec[2] ){ 
            std::vector<std::string_view> xerrVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xerrVec.push_back(line->xerrup);
            }
            lheDs[currInd] = vecStoD( xerrVec );
             ++currInd; }
        if( boolVec[3] ){ 
            std::vector<std::string_view> xmaxVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                xmaxVec.push_back(line->xmaxup);
            }
            lheDs[currInd] = vecStoD( xmaxVec );
             ++currInd; } 
        for( size_t k = 0 ; k < lheAOS.subProcs.size() ; ++k )
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
        if( boolVec[0] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->idbmup[0], lheFile.getInit()->getHead()->idbmup[1] } ); ++currInd; }
        if( boolVec[1] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->pdfgup[0], lheFile.getInit()->getHead()->pdfgup[1] } ); ++currInd; }
        if( boolVec[2] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->pdfsup[0], lheFile.getInit()->getHead()->pdfsup[1] } ); ++currInd; }
        if( boolVec[3] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->idwtup } ); ++currInd; }
        if( boolVec[4] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->nprup } ); ++currInd; }
        if( boolVec[5] ){ 
            std::vector<std::string_view> lprVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                lprVec.push_back(line->lprup);
            }
            lheDs[currInd] = vecStoI( lprVec );
             ++currInd; }
        for( size_t k = 0 ; k < lheAOS.subProcs.size() ; ++k )
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

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> lheValInts( lheNode& lheFile, std::vector<std::string_view> statVec,
    lheRetInts vals = lheRetInts() )
    {
        auto boolVec = vals.getBools();
        const int noVals = std::count(boolVec.begin(), boolVec.end(), true);
        auto lheAOS = transLHE( lheFile, statVec );
        auto lheIs = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(noVals * lheAOS.subProcs.size() );
        std::vector<std::shared_ptr<std::vector<int>>> &lheDs = *lheIs;
        int currInd = 0;
        if( boolVec[0] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->idbmup[0], lheFile.getInit()->getHead()->idbmup[1] } ); ++currInd; }
        if( boolVec[1] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->pdfgup[0], lheFile.getInit()->getHead()->pdfgup[1] } ); ++currInd; }
        if( boolVec[2] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->pdfsup[0], lheFile.getInit()->getHead()->pdfsup[1] } ); ++currInd; }
        if( boolVec[3] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->idwtup } ); ++currInd; }
        if( boolVec[4] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->nprup } ); ++currInd; }
        if( boolVec[5] ){ 
            std::vector<std::string_view> lprVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                lprVec.push_back(line->lprup);
            }
            lheDs[currInd] = vecStoI( lprVec );
             ++currInd; }
        for( size_t k = 0 ; k < lheAOS.subProcs.size() ; ++k )
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

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> lheValInts( lheNode& lheFile, 
    sortFcn sorter, 
    std::vector<std::string_view> statVec = {"-1", "1"}, lheRetInts vals = lheRetInts() )
    {
        auto boolVec = vals.getBools();
        const int noVals = std::count(boolVec.begin(), boolVec.end(), true);
        auto lheAOS = transLHE( lheFile, sorter, statVec );
        auto lheIs = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(noVals * lheAOS.subProcs.size() );
        std::vector<std::shared_ptr<std::vector<int>>> &lheDs = *lheIs;
        int currInd = 0;
        if( boolVec[0] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->idbmup[0], lheFile.getInit()->getHead()->idbmup[1] } ); ++currInd; }
        if( boolVec[1] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->pdfgup[0], lheFile.getInit()->getHead()->pdfgup[1] } ); ++currInd; }
        if( boolVec[2] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->pdfsup[0], lheFile.getInit()->getHead()->pdfsup[1] } ); ++currInd; }
        if( boolVec[3] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->idwtup } ); ++currInd; }
        if( boolVec[4] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->nprup } ); ++currInd; }
        if( boolVec[5] ){ 
            std::vector<std::string_view> lprVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                lprVec.push_back(line->lprup);
            }
            lheDs[currInd] = vecStoI( lprVec );
             ++currInd; }
        for( size_t k = 0 ; k < lheAOS.subProcs.size() ; ++k )
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

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> lheValInts( lheNode& lheFile, 
    statSort sorter, 
    std::vector<std::string_view> statVec = {"-1", "1"}, lheRetInts vals = lheRetInts() )
    {
        auto boolVec = vals.getBools();
        const int noVals = std::count(boolVec.begin(), boolVec.end(), true);
        auto lheAOS = transLHE( lheFile, sorter, statVec );
        auto lheIs = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(noVals * lheAOS.subProcs.size() );
        std::vector<std::shared_ptr<std::vector<int>>> &lheDs = *lheIs;
        int currInd = 0;
        if( boolVec[0] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->idbmup[0], lheFile.getInit()->getHead()->idbmup[1] } ); ++currInd; }
        if( boolVec[1] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->pdfgup[0], lheFile.getInit()->getHead()->pdfgup[1] } ); ++currInd; }
        if( boolVec[2] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->pdfsup[0], lheFile.getInit()->getHead()->pdfsup[1] } ); ++currInd; }
        if( boolVec[3] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->idwtup } ); ++currInd; }
        if( boolVec[4] ){ lheDs[currInd] = vecStoI( { lheFile.getInit()->getHead()->nprup } ); ++currInd; }
        if( boolVec[5] ){ 
            std::vector<std::string_view> lprVec( lheFile.getInit()->getLines().size() );
            for( auto line : lheFile.getInit()->getLines() )
            {
                lprVec.push_back(line->lprup);
            }
            lheDs[currInd] = vecStoI( lprVec );
             ++currInd; }
        for( size_t k = 0 ; k < lheAOS.subProcs.size() ; ++k )
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

#endif
