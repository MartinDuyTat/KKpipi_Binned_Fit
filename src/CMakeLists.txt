add_library(KKpipiBinnedFit Amplitude.cpp AmplitudeBinningOptimizer.cpp AmplitudePhaseSpace.cpp Bin.cpp BinList.cpp CPParameters.cpp DDecayParameters.cpp Event.cpp EventList.cpp FitGamma.cpp Fitter.cpp Gamma.cpp Generator.cpp KKpipiFit.cpp KKpipiMath.cpp Likelihood.cpp NaivePhaseSpace.cpp PhaseSpaceParameterisation.cpp RectangularPhaseSpace.cpp SophisticatedPhaseSpace.cpp XYLikelihood.cpp)

target_include_directories(KKpipiBinnedFit PUBLIC ../include)

target_link_libraries(KKpipiBinnedFit PUBLIC ROOT::Physics ROOT::Minuit2 ROOT::Gpad)

target_link_libraries(KKpipiBinnedFit PUBLIC -ldl)

if(OpenMP_CXX_FOUND)
    target_link_libraries(KKpipiBinnedFit PUBLIC OpenMP::OpenMP_CXX)
endif()
