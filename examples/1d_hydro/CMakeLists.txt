# Copyright (c) 2007-2012 Hartmut Kaiser
# Copyright (c) 2011-2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

set(example_programs
   )

if(HPX_WITH_CXX11_LAMBDAS)
  set(example_programs
      ${example_programs}
      1d_hydro_upwind)
  set(1d_hydro_upwind_FLAGS DEPENDENCIES iostreams_component)
endif()

foreach(example_program ${example_programs})
  set(sources
      ${example_program}.cpp)

  source_group("Source Files" FILES ${sources})

  # add example executable
  add_hpx_executable(${example_program}
                     SOURCES ${sources}
                     ${${example_program}_FLAGS}
                     FOLDER "Examples/1D Hydro/${example_program}")

  # add a custom target for this example
  add_hpx_pseudo_target(examples.1d_hydro.${example_program})

  # make pseudo-targets depend on master pseudo-target
  add_hpx_pseudo_dependencies(examples.1d_hydro
                              examples.1d_hydro.${example_program})

  # add dependencies to pseudo-target
  add_hpx_pseudo_dependencies(examples.1d_hydro.${example_program}
                              ${example_program}_exe)
endforeach()

