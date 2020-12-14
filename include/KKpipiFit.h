// Martin Duy Tat 12th December 2020
/**
 * KKpipiFit is a namespace with helper functions for doing the fitting and pull studies
 */

#ifndef KKPIPIFIT
#define KKPIPIFIT

#include<string>
#include"KKpipiFit.h"
#include"PhaseSpaceParameterisation.h"
#include"NaivePhaseSpace.h"
#include"RectangularPhaseSpace.h"
#include"SophisticatedPhaseSpace.h"
#include"AmplitudePhaseSpace.h"
#include"BinList.h"
#include"CPParameters.h"
#include"Gamma.h"
#include"FitGamma.h"
#include"TTree.h"
#include"TFile.h"

namespace KKpipiFit {
  /**
   * Function for allocating the correct binning scheme
   * @param binning_choice Name of binning scheme chosen
   * @return Pointer to a PhaseSpaceParameterisation object, where the chosen binning scheme is allocated
   */
  PhaseSpaceParameterisation* PickBinningScheme(const std::string &binning_choice);
  /**
   * Function for loading events into bins from trees
   * @param TreeBplus TTree with \f$B^+\f$ events
   * @param TreeBminus TTree with \f$B^-\f$ events
   * @param binlist BinList object to load input data into
   * @param StartEvent Index of first event
   * @param SampleSize Size of datasample
   */
  void LoadTreesIntoBins(TTree *TreeBplus, TTree *TreeBminus, BinList &binlist, const int &StartEvent = 0, const int &TotalEvents = -1);
  /**
   * Function for loading input data into the BinList object
   * @param Bplusfile TFile with input data file with \f$B^+\f$ events
   * @param Bminusfile TFile with input data file with \f$B^-\f$ events
   * @param binlist BinList object to load input data into
   * @param StartEvent Index of first event
   * @param SampleSize Size of datasample
   */
  void LoadInputDataIntoBins(const std::string &BplusFilename, const std::string &BminusFilename, BinList &binlist, const int &StartEvent = 0, const int &TotalEvents = -1);
  /**
   * Function for printing fit results for \f$x_\pm\f$ and \f$y_\pm\f$
   * @param cpparameters CPParameters object with fit results
   */
  void PrintXY (const CPParameters &cpparameters);
  /**
   * Function for printing fit results for \f$r_B\f$, \f$\delta_B\f$ and \f$\gamma\f$
   * @param gammaparams Gamma object with fit results
   */
  void PrintGamma (const Gamma &gammaparams);
  /**
   * Function that asks if user wants to draw contours of \f$r_B\f$, \f$\delta_B\f$ and \f$\gamma\f$ after fitting
   * If user inputs "yes", the contours are saved as PNG files
   * @param fitgamma FitGamma object with fit results
   */
  void DrawContours(const FitGamma &fitgamma);
}

#endif
