cmake_minimum_required(VERSION 3.17)
project(KKpipi_Binned_Fit
	VERSION 1.0
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options (-Wall)
add_compile_options (-Wunused)
add_compile_options (-Wextra)
add_compile_options (-Werror)

add_subdirectory(apps)
add_subdirectory(include)
add_subdirectory(src)