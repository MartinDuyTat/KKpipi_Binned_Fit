add_executable(FitNaive FitNaive.cpp)
add_executable(FitRectangular FitRectangular.cpp)
add_executable(FitSophisticated FitSophisticated.cpp)
add_executable(ComputeDParameters ComputeDParameters.cpp)
add_executable(PullStudy PullStudy.cpp)
add_executable(PlotDParameters PlotDParameters.cpp)
add_executable(BinNumberStudy BinNumberStudy.cpp)
add_executable(ComputeMeanPhases ComputeMeanPhases.cpp)
add_executable(MapStrongPhases MapStrongPhases.cpp)

target_link_libraries(FitNaive PUBLIC KKpipi)
target_link_libraries(FitNaive PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitRectangular PUBLIC KKpipi)
target_link_libraries(FitRectangular PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitSophisticated PUBLIC KKpipi)
target_link_libraries(FitSophisticated PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(BinNumberStudy PUBLIC KKpipi)
target_link_libraries(BinNumberStudy PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(ComputeDParameters PUBLIC KKpipi)
target_link_libraries(ComputeDParameters PUBLIC ROOT::Tree)

target_link_libraries(PlotDParameters PUBLIC KKpipi)

target_link_libraries(PullStudy PUBLIC KKpipi)
target_link_libraries(PullStudy PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(ComputeMeanPhases PUBLIC KKpipi)
target_link_libraries(ComputeMeanPhases PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)
target_link_libraries(MapStrongPhases PUBLIC KKpipi)

set_target_properties(FitNaive PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
set_target_properties(FitRectangular PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
set_target_properties(FitSophisticated PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
set_target_properties(ComputeDParameters PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
set_target_properties(PullStudy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
set_target_properties(PlotDParameters PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
set_target_properties(BinNumberStudy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
set_target_properties(ComputeMeanPhases PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
set_target_properties(MapStrongPhases PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
