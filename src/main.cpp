#include "snapshot.hpp"

int main(int argc, char* argv[]){
    Snapshot ss;
    std::cout<<CONSOLE_YELLOW<<"Initialized paramters ..."<<CONSOLE_WHITE<<std::endl;
    if (argv[4])
        ss.initializeParameters(std::stod(argv[3]), std::stod(argv[4]));
    else
        ss.initializeParameters(std::stoi(argv[3]));
    std::cout<<CONSOLE_YELLOW<<"Build alpha shape ..."<<CONSOLE_WHITE<<std::endl;
    ss.buildAlphaShape(argv[1]);
    std::cout<<CONSOLE_YELLOW<<"Read filtration ..."<<CONSOLE_WHITE<<std::endl;
    ss.readFiltration(argv[2]);    
    std::cout<<CONSOLE_YELLOW<<"Preprocess ..."<<CONSOLE_WHITE<<std::endl;
    ss.preprocess();

    std::cout<<CONSOLE_YELLOW<<"Take snapshots ..."<<CONSOLE_WHITE<<std::endl;
    ss.takeSnapshots();
    std::cout<<CONSOLE_YELLOW<<"Write ..."<<CONSOLE_WHITE<<std::endl;
    ss.write();
    
    std::cout<<CONSOLE_YELLOW<<"Debug ..."<<CONSOLE_WHITE<<std::endl;
    ss.debug();

	return 0;
}
