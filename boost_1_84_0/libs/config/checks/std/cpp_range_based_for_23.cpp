//  This file was automatically generated on Fri Oct 13 19:09:38 2023
//  by libs/config/tools/generate.cpp
//  Copyright John Maddock 2002-21.
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/config for the most recent version.//
//  Revision $Id$
//

#ifdef __has_include
#if __has_include(<version>)
#include <version>
#endif
#endif

#ifndef __cpp_range_based_for
#error "Macro << __cpp_range_based_for is not set"
#endif

#if __cpp_range_based_for < 202211
#error "Macro __cpp_range_based_for had too low a value"
#endif

int main( int, char *[] )
{
   return 0;
}

