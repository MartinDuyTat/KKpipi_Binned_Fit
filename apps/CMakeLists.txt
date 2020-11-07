add_executable(Fitting Fitting.cpp)
add_executable(ComputeDParameters ComputeDParameters.cpp)
add_executable(PoolStudy PoolStudy.cpp)

target_link_libraries(Fitting PUBLIC KKpipi)
target_link_libraries(Fitting PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(ComputeDParameters PUBLIC KKpipi)

target_link_libraries(PoolStudy PUBLIC KKpipi)
target_link_libraries(PoolStudy PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

set_target_properties(Fitting PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
set_target_properties(ComputeDParameters PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
set_target_properties(PoolStudy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
