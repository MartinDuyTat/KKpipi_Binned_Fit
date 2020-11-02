add_executable(Fitting Fitting.cpp)

target_link_libraries(Fitting PUBLIC KKpipi)

set_target_properties(Fitting PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)