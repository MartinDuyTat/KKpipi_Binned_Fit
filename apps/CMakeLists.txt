add_executable(Fitting Fitting.cpp)

target_link_libraries(Fitting PUBLIC KKpipi)

find_package(ROOT 6.22 CONFIG REQUIRED)

target_link_libraries(Fitting PUBLIC ROOT::Physics)

set_target_properties(Fitting PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../bin)
