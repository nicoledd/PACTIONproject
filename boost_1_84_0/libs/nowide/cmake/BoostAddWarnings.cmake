# Enable warnings for the given target
# Arguments:
#   target: Non-interface library or executable
#   level: Which warnings to add. Valid values: off, on, all, extra
#   warningsAsErrors: Optional, ON to treat warnings as errors
# Optionally can pass a bool as warnings-as-errors
function(boost_add_warnings target level)
  if(ARGC GREATER 2)
    set(warningsAsErrors ${ARGV2})
  else()
    set(warningsAsErrors OFF)
  endif()
  set(allowed_levels off on all extra pedantic)
  if(NOT level IN_LIST allowed_levels)
    message(FATAL_ERROR "${level} is not a valid warning level (${allowed_levels})")
  endif()
  if(MSVC)
    set(warn_off /W0)
    set(warn_on /W3)
    foreach(_lvl IN ITEMS all extra pedantic)
      set(warn_${_lvl} /W4)
    endforeach()
    set(werror /WX)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(warn_off -w)
    set(warn_on  -Wall)
    set(warn_all -Wall)
    set(warn_extra -Wall -Wextra)
    set(warn_pedantic -Wall -Wextra -pedantic)
    set(werror -Werror)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(warn_off -w0)
    foreach(_lvl IN ITEMS on all extra pedantic)
      set(warn_${_lvl} -w1)
    endforeach()
    set(werror "")
  endif()
  target_compile_options(${target} PRIVATE ${warn_${level}})
  if(warningsAsErrors)
    target_compile_options(${target} PRIVATE ${werror})
  endif()
endfunction()
