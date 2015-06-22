
option(BUILD_SICONOS_PLUGIN "Building SiconosPlugin" OFF)

if(NOT BUILD_SICONOS_PLUGIN)
  return()
endif()

set(target CnoidSiconosPlugin)

set(sources
  SiconosPlugin.cpp
  SPSimulatorItem.cpp
  SPConstraintForceSolver.cpp
  SPCore.cpp
  )

set(headers
  )

include_directories(
  /usr/local/include/Siconos/Numerics
)

link_directories(/usr/local/lib/x86_64-linux-gnu)

 
make_gettext_mofiles(${target} mofiles)
add_cnoid_plugin(${target} SHARED ${sources} ${headers} ${mofiles})
target_link_libraries(${target} CnoidBodyPlugin)
apply_common_setting_for_plugin(${target} "${headers}")


target_link_libraries(${target} SiconosNumerics)


if(ENABLE_PYTHON)
#  add_subdirectory(python)
endif()
