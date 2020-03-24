cmake_minimum_required(VERSION 3.2)

if(WIN32 OR WIN64)
    message(STATUS "Windows detected")
    project(tf_gd_lib VERSION 0.0.1000)
    #set(BOOST_ROOT "C:/Boost")
else()
    project(tf_gd_lib VERSION 0.0.$ENV{TRAVIS_BUILD_NUMBER})
endif()

find_package(Boost COMPONENTS unit_test_framework REQUIRED)

if (WIN32 OR WIN64)
    # disable autolinking in boost
    add_definitions( -DBOOST_ALL_NO_LIB )

    set(Boost_USE_STATIC_LIBS OFF CACHE BOOL "use static libraries from Boost")

    add_definitions(-DBOOST_ALL_DYN_LINK )
endif()

#add_executable(tf_gd_lib_cli main.cpp UnitSpline.h UnitTableFunctions.h UnitGradDescent.h)

add_library(tf_gd_lib SHARED UnitSpline.h UnitSpline.cpp 
                             UnitTableFunctions.h UnitTableFunctions.cpp 
                             UnitGradDescent.h UnitGradDescent.cpp)

add_executable(tf_gd_lib_tests tests.cpp UnitSpline.h UnitTableFunctions.h UnitGradDescent.h)

# add tf_gd_lib_cli if it will be used
if(WIN32 OR WIN64)
    set_target_properties(tf_gd_lib tf_gd_lib_tests PROPERTIES
            CXX_STANDARD 17
            CXX_STANDARD_REQUIRED ON
            COMPILE_OPTIONS "-Wpedantic;-Wall;-Wextra")
else()
    set_target_properties(tf_gd_lib tf_gd_lib_tests PROPERTIES
            CXX_STANDARD 17
            CXX_STANDARD_REQUIRED ON
            COMPILE_OPTIONS "-Wpedantic;-Wall;-Wextra")

    set(CMAKE_BUILD_TYPE Release)
endif()

#target_link_libraries(tf_gd_lib_cli
#    tf_gd_lib
#)

if(WIN32 OR WIN64)
	set_target_properties(tf_gd_lib_tests PROPERTIES
		COMPILE_DEFINITIONS BOOST_TEST_DYN_LINK
		INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR}
	)
endif()


target_link_libraries(tf_gd_lib_tests
    ${Boost_LIBRARIES}
	tf_gd_lib
)

install(TARGETS tf_gd_lib LIBRARY DESTINATION bin)
#install(TARGETS tf_gd_lib RUNTIME DESTINATION bin)

set(CPACK_GENERATOR DEB)

set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")

set(CPACK_PACKAGE_CONTACT "znseday@yandex.ru")

set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

include(CPack)

enable_testing()
add_test(tf_gd_lib_tests_ tf_gd_lib_tests)
