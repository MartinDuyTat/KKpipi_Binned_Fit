add_executable(Fitting Fitting.cpp)

add_executable(ComputeDParameters ComputeDParameters.cpp)

target_link_libraries(Fitting PUBLIC KKpipi)

target_link_libraries(Fitting PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(ComputeDParameters PUBLIC KKpipi)

set_target_properties(Fitting PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)

set_target_properties(ComputeDParameters PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
