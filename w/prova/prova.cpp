#include <iostream>

#include "Scenario.hpp"
#include "ObservableProcessor.hpp"


int main(void)
{
    std::cout << "Prova." << std::endl;

    std::cout << "Creating Scenario." << std::endl;
    Scenario* Scen;
    Scen = new Scenario();

    std::cout << "Scenario Name: " << Scen->GetName() << std::endl;

    ObservableProcessor* Proc;
    Proc = new ObservableProcessorA();

    Scen->Run(Proc);

    std::cout << "Now the name of the scenario has changed to:" << Scen->GetName() << std::endl;
    std::cout << "Fi. " << std::endl;
    return 0;
}
