#ifndef PLOTS_HPP_
#define PLOTS_HPP_

/**
    \file   Plots.hpp
    \brief  Header file for th Plots class file.
*/

/**
    \brief  The Plots class file defines the stadard types of plots, so that they can be reused by the 
            different classes. This class relies 100% on the matplotlibcpp plotting class.
*/

#include "matplotlibcpp.hpp"

// plot size
#define PLT_SIZE_X 1200
#define PLT_SIZE_Y  780

// colors
#define PLT_COL0 "r"
#define PLT_COL1 "g"
#define PLT_COL2 "b"
#define PLT_COL3 "m"
#define PLT_COL4 "c"
#define PLT_COL5 "y"
#define PLT_COL6 "k"

// lines
#define PLT_LIN0 "o-"
#define PLT_LIN1 "x--"
#define PLT_LIN2 "*:"
#define PLT_LIN3 "+-."

//
#define PLT_MAX_N 500 /* maximum number of samples in plot */

class Plots
{
  public:
    static const char *col[]; 
    static const char *lin[];

    static void PlotSingle( const std::string fname,
                            const std::vector<std::vector<double>> x,
                            const std::vector<std::vector<double>> y,
                            const std::vector<std::string> xlabels,
                            const std::vector<std::string> ylabels,
                            const std::vector<std::string> titles);
    static void PlotMultiple( const std::string fname, 
                              const std::vector<std::vector<double>> x,
                              const std::vector<std::vector<double>> y,
                              const std::vector<std::string> xlabels,
                              const std::vector<std::string> ylabels,
                              const std::vector<std::string> titles);
    static void PlotImage(const std::string fname,
                          const std::vector<std::vector<double>> x,
                          const std::vector<std::vector<double>> y,
                          const std::vector<std::vector<double>> z,
                          const std::string xlabel,
                          const std::string ylabel,
                          const std::string zlabel,
                          const std::string title);
};

#endif // PLOTS_HPP_
