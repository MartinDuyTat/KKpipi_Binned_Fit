add_library(KKpipi Amplitude.cpp Bin.cpp BinList.cpp CPParameters.cpp DDecayParameters.cpp Event.cpp EventList.cpp Fitter.cpp Generator.cpp Likelihood.cpp PhaseSpaceParameterisation.cpp)

target_include_directories(KKpipi PUBLIC ../include)

target_link_libraries(KKpipi PUBLIC ROOT::Physics ROOT::Minuit2)

target_link_libraries(KKpipi PUBLIC -ldl)
