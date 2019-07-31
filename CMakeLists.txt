cmake_minimum_required(VERSION 2.8.9)
  project(generic_pixel-level_editing)

  # Tell CMake to run moc when necessary:
  set(CMAKE_AUTOMOC ON)
  # As moc files are generated in the binary dir, tell CMake
  # to always look for includes there:
  set(CMAKE_INCLUDE_CURRENT_DIR ON)

  ############################# #####################################################
  #     PATHS
  ##################################################################################

  # Path for output executable files
  IF(NOT (IS_DIRECTORY "${generic_pixel-level_editing_SOURCE_DIR}/bin"))
     file(MAKE_DIRECTORY "${generic_pixel-level_editing_SOURCE_DIR}/bin")
  endif()
  set(EXECUTABLE_OUTPUT_PATH "${generic_pixel-level_editing_SOURCE_DIR}/bin")

  ################################
  # PATH Qt: http://www.qt.io/download/
  ################################


  ################################
  # PATH OpenCV
  ################################


  ##################################################################################
  #     Qt
  ##################################################################################

  # Widgets finds its own dependencies (QtGui and QtCore).

  set( Qt5Core_DIR ${Qt5_DIR}/../Qt5Core)
  set( Qt5Widgets_DIR  ${Qt5_DIR}/../Qt5Widgets)
  set( Qt5Gui_DIR  ${Qt5_DIR}/../Qt5Gui)

  # Find packages Qt
  # find_package(Qt5)
  find_package(Qt5Core)
  find_package(Qt5Widgets)
  find_package(Qt5Gui)

  if(${Qt5Core_FOUND})
     # The Qt5Widgets_INCLUDES also includes the include directories for
     # dependencies QtCore and QtGui
     include_directories(${Qt5Core_INCLUDES})
     include_directories(${Qt5Widgets_INCLUDES})

     # We need add -DQT_WIDGETS_LIB when using QtWidgets in Qt 5.
     add_definitions(${Qt5Core_DEFINITIONS})
     add_definitions(${Qt5Widgets_DEFINITIONS})

     # Executables fail to build with Qt 5 in the default configuration
     # without -fPIE. We add that here.
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Core_EXECUTABLE_COMPILE_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
  endif()

  ##################################################################################
  #     OpenCV
  ##################################################################################

  find_package(OpenCV)
  if(${OpenCV_FOUND})
     include_directories(${OpenCV_INCLUDES})
     link_directories(${OpenCV})
  endif()


  ##################################################################################
  #     Boost
  ##################################################################################

  find_package(Boost COMPONENTS system filesystem program_options)

  ##################################################################################
  #     Compiler FLAGS
  ##################################################################################
 # set(CMAKE_CXX_COMPILER_ID "GNU")
 if ( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") )
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mtune=native")
    if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
       set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Dsrandom=srand -Drandom=rand -D_USE_MATH_DEFINES -D_WIN")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message("Using Visual Studio, are you sure?")
endif()

message(STATUS "Compiler  = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "System    = ${CMAKE_SYSTEM_NAME}")
message(STATUS "Prefix    = ${CMAKE_PREFIX_PATH}")
message(STATUS "Flags     = ${CMAKE_CXX_FLAGS}")
message(STATUS "Build     = ${CMAKE_BUILD_TYPE}")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
   message(STATUS "R.Flags   = ${CMAKE_CXX_FLAGS_RELEASE}")
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
   message(STATUS "D.Flags   = ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

  ##################################################################################
  # Eigen: Download from github
  ##################################################################################
  include(ExternalProject)

  set(EXTERNAL_INSTALL_LOCATION ${generic_pixel-level_editing_SOURCE_DIR}/libs)

	 ExternalProject_Add(Eigen
	   GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
	   PREFIX ${EXTERNAL_INSTALL_LOCATION}/Eigen
	   BUILD_COMMAND ""
	   CONFIGURE_COMMAND ""
	   INSTALL_COMMAND ""
	 )

  ##################################################################################
  #     paths: include and lib
  ##################################################################################

  include_directories(${generic_pixel-level_editing_SOURCE_DIR}/src)
  include_directories(${generic_pixel-level_editing_SOURCE_DIR}/libs)
  include_directories(${EXTERNAL_INSTALL_LOCATION}/Eigen/src/Eigen)


  ##################################################################################
  #     setting variables for source code and libs
  ##################################################################################
  set(GUI_SRCS
        ${generic_pixel-level_editing_SOURCE_DIR}/src/gui/labelimage.h
        ${generic_pixel-level_editing_SOURCE_DIR}/src/gui/imagemousebrush.h
        ${generic_pixel-level_editing_SOURCE_DIR}/src/gui/multiimageviewer.h
        ${generic_pixel-level_editing_SOURCE_DIR}/src/gui/widgetfilter.h
        ${generic_pixel-level_editing_SOURCE_DIR}/src/gui/windowfilter.h
        ${generic_pixel-level_editing_SOURCE_DIR}/libs/imageHDR.cpp
        ${generic_pixel-level_editing_SOURCE_DIR}/libs/imageHDR.h
  )

  set(GUI_LIBS Qt5::Gui Qt5::Core Qt5::Widgets)

  set(SUPERPIXELS_SRCS
        ${generic_pixel-level_editing_SOURCE_DIR}/src/superpixels/superpixel.cpp
        ${generic_pixel-level_editing_SOURCE_DIR}/src/superpixels/superpixels.cpp
        ${generic_pixel-level_editing_SOURCE_DIR}/src/denseLabeling/denseLabeling.cpp
        ${generic_pixel-level_editing_SOURCE_DIR}/libs/lib_lsc/LSC.cpp
  )

  ##################################################################################
  #     src
  ##################################################################################
  add_subdirectory(src/apps/dehazing)
  add_subdirectory(src/apps/depthoffield)
  add_subdirectory(src/apps/tonemapping)

