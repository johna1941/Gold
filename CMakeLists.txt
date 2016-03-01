# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(Gold1)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include /usr/local/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(Gold1 Gold1.cc ${sources} ${headers})
target_link_libraries(Gold1 ${Geant4_LIBRARIES} /usr/local/lib/libgsl.dylib)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Gold1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(Gold1_SCRIPTS
  Gold1.in
  Gold1.out
  init_vis.mac
  init_run.mac
  vis.mac
  gun.mac
  gui.mac
  icons.mac
  harrypotter5.png
  harrypotter1.png
  monkey.png
  illuminati.png
  )
#
# protongun.mac, photongun.mac, squaregun.mac gui.mac, icons.mac, and 
# harrypotter5.png have been added to CMakeLists.txt to ensure that they 
# are copied across to the build directory when the "make" command is 
# entered. harrypotter5.png is a 128x128 pixel icon used as the icon 
# for /run/beamOn 100. 
# John: this is just humour on our part, you're welcome to ignore it. 
#
#
foreach(_script ${Gold1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(Gold1 DEPENDS Gold1)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS Gold1 DESTINATION bin)


