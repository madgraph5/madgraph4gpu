#include "teawREX.hpp"
#include <cstdlib>
#include <typeinfo>

std::shared_ptr<std::vector<double>> meEval( std::vector<double>& x, std::vector<double>& y){
    int random = rand() % 10;
    if( random == 0 ){ random = 11; }
    auto thisIsIt = std::make_shared<std::vector<double>>( y.size(), random );
    return thisIsIt;
}

std::shared_ptr<std::vector<size_t>> sortFunc(std::vector<std::string_view> arguments){
    return REX::stoiSort(arguments);
}

std::shared_ptr<std::vector<size_t>> sorterFunc(std::string_view dummy, std::vector<std::string_view> arguments){
    return REX::stoiSort(arguments);
}

int main( int argc, char* argv[] ){
    std::string lheFilePath;
    
    // READ COMMAND LINE ARGUMENTS
    for( int arg = 0; arg < argc; arg++ )
    {
        auto currArg = std::string( argv[arg] );
        if( currArg.substr(0,9) == "--lhefile" || currArg.substr(0,4) == "-lhe" )
        {
            lheFilePath = currArg.substr( currArg.find( "=" ) + 1 ); 
        }
    }


    std::string currPath = argv[0];
    auto sembler = std::function<std::shared_ptr<std::vector<size_t>>(std::vector<std::string_view>)>(sortFunc);
    auto sembler2 = std::function<std::shared_ptr<std::vector<size_t>>(std::string_view, std::vector<std::string_view>)>(sorterFunc);
    auto lheFile = REX::filePuller(lheFilePath);
    //std::cout << lheFile->substr(0, 1) << "\n";
    //std::cout << bool(lheFile->compare(0, 1, "<")) << "\n";
    //std::cout << lheFile->substr(1968, 1999 - 1968) << "\n";
    auto parseLhe = REX::lheNode(*lheFile);
    //std::cout << *parseLhe.nodeWriter() << "\n";
    auto treeMan = parseLhe.getTree();
    //std::cout << parseLhe.getChildren().size() << " & "  << parseLhe.getEvents().size() << " & "  << treeMan.getChildren()->size() << "\n";
    auto proceses = REX::lheReOrder(parseLhe, {"-1", "1", "2"} );
    auto processes2 = REX::lheEvReOrder(parseLhe, {"-1", "1", "2"} );
    //std::cout << proceses.size() << " & " << processes2.size() << "\n";
    bool comp = REX::evProcComp( *parseLhe.getEvents()[0], *parseLhe.getEvents()[1], {"-1", "1"} );
    if( comp ){ std::cout << "true\n"; }
    else{ std::cout << "false\n"; }
    auto evlist = REX::evProcessPull( parseLhe, {"-1", "1"} );
    //auto evsVals = lheValDoubles(parseLhe);
    auto evsVals = lheValDoubles(parseLhe, sembler2);
    int siz = 0;
    for( auto& ev : *evsVals ){
        siz += ev->size();
    }
    std::cout << evsVals->size() << "\n";
    std::cout << siz << "\n";
    return 0;

}