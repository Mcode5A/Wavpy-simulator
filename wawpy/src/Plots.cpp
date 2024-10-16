#ifndef PLOTS_CPP_
#define PLOTS_CPP_

#include "Plots.hpp"

#include <iostream>

const char *mychar[] = {"a", "b"};
// Colors
const char *Plots::col[] = {PLT_COL0, PLT_COL1, PLT_COL2, PLT_COL3, PLT_COL4, PLT_COL5, PLT_COL6};

// Ĺine styles
const char *Plots::lin[] = {PLT_LIN0, PLT_LIN1, PLT_LIN2, PLT_LIN3};

void Plots::PlotSingle( const std::string fname,
                        const std::vector<std::vector<double>> x,
                        const std::vector<std::vector<double>> y,
                        const std::vector<std::string> xlabels,
                        const std::vector<std::string> ylabels,
                        const std::vector<std::string> titles)
{
  std::cout << "PlotSingle is TBD." << std::endl;
  int N = x.size();
  matplotlibcpp::figure_size(PLT_SIZE_X, PLT_SIZE_Y);

  for (int k=0; k<N; k++)
  {
    matplotlibcpp::named_plot(ylabels[k], x[k], y[k], std::string(Plots::col[k%N]) + std::string(Plots::lin[k%N]));
  }
  matplotlibcpp::xlabel(xlabels[0]);
  matplotlibcpp::title(titles[0]);
  matplotlibcpp::grid(true);
  matplotlibcpp::legend();
 
  // Visualize plot or save to file
  if (fname.empty())
  {
    matplotlibcpp::show();
  }
  else
  {
    matplotlibcpp::save(fname);
    matplotlibcpp::close();
  }

  return;
}


void Plots::PlotMultiple( const std::string fname, 
                          const std::vector<std::vector<double>> x,
                          const std::vector<std::vector<double>> y,
                          const std::vector<std::string> xlabels,
                          const std::vector<std::string> ylabels,
                          const std::vector<std::string> titles)
{
  //std::cout << "PlotMultiple is TBD." << std::endl;
  
  int N = x.size();
  matplotlibcpp::figure_size(PLT_SIZE_X, PLT_SIZE_Y);
  
  for (int k=0; k<N; k++)
  {
    matplotlibcpp::subplot(N, 1, k+1);
    matplotlibcpp::plot(x[k], y[k], std::string(Plots::col[k%N]) + std::string(Plots::lin[k%N]));
    matplotlibcpp::ylabel(ylabels[k]);
    matplotlibcpp::xlabel(xlabels[k]);
    matplotlibcpp::title(titles[k]);
    matplotlibcpp::grid(true);
  }

  // Visualize plot or save to file
  if (fname.empty())
  {
    matplotlibcpp::show();
  }
  else
  {
    matplotlibcpp::save(fname);
    matplotlibcpp::close();
  }

  return;
}


void Plots::PlotImage(const std::string fname,
                      const std::vector<std::vector<double>> x,
                      const std::vector<std::vector<double>> y,
                      const std::vector<std::vector<double>> z,
                      const std::string xlabel,
                      const std::string ylabel,
                      const std::string zlabel,
                      const std::string title)
{
  matplotlibcpp::plot_surface(x,y,z);
  std::cout << "fname= " << fname << std::endl << std::endl << std::endl;
  
  if (fname.empty())
  {
    matplotlibcpp::show();
  }
  else
  {
    matplotlibcpp::save(fname);
  }

  return;
}


#endif // PLOTS_CPP_
