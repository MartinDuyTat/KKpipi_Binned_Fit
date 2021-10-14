set(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(BinningOptimizer BinningOptimizer.cpp)
add_executable(BinningQ BinningQ.cpp)
add_executable(BinNumberStudy BinNumberStudy.cpp)
add_executable(ComputeDParameters ComputeDParameters.cpp)
add_executable(ComputeMeanPhases ComputeMeanPhases.cpp)
add_executable(ConvertFlatCoordinates ConvertFlatCoordinates.cpp)
add_executable(CopyRedundantCS CopyRedundantCS.cpp)
add_executable(FastBinningTest FastBinningTest.cpp)
add_executable(FitKKpipi FitKKpipi.cpp)
add_executable(GenerateFlat GenerateFlat.cpp)
add_executable(MapStrongPhases MapStrongPhases.cpp)
add_executable(PullStudy PullStudy.cpp)
add_executable(PlotDParameters PlotDParameters.cpp)

target_link_libraries(BinningOptimizer PUBLIC KKpipiBinnedFit)
target_link_libraries(BinningOptimizer PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(BinningQ PUBLIC KKpipiBinnedFit)
target_link_libraries(BinningQ PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(BinNumberStudy PUBLIC KKpipiBinnedFit)
target_link_libraries(BinNumberStudy PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(ComputeDParameters PUBLIC KKpipiBinnedFit)
target_link_libraries(ComputeDParameters PUBLIC ROOT::Tree)

target_link_libraries(ComputeMeanPhases PUBLIC KKpipiBinnedFit)
target_link_libraries(ComputeMeanPhases PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(ConvertFlatCoordinates PUBLIC KKpipiBinnedFit)
target_link_libraries(ConvertFlatCoordinates PUBLIC ROOT::Tree)

target_link_libraries(CopyRedundantCS PUBLIC KKpipiBinnedFit)
target_link_libraries(CopyRedundantCS PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FastBinningTest PUBLIC KKpipiBinnedFit)
target_link_libraries(FastBinningTest PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitKKpipi PUBLIC KKpipiBinnedFit)
target_link_libraries(FitKKpipi PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(GenerateFlat PUBLIC KKpipiBinnedFit)
target_link_libraries(GenerateFlat PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(MapStrongPhases PUBLIC KKpipiBinnedFit ROOT::Tree)

target_link_libraries(PlotDParameters PUBLIC KKpipiBinnedFit)
target_link_libraries(PlotDParameters PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(PullStudy PUBLIC KKpipiBinnedFit)
target_link_libraries(PullStudy PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

install(TARGETS BinningOptimizer BinningQ BinNumberStudy ComputeDParameters ComputeMeanPhases ConvertFlatCoordinates CopyRedundantCS FastBinningTest FitKKpipi GenerateFlat MapStrongPhases PlotDParameters PullStudy DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
