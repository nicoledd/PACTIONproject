// Copyright Louis Dionne 2013-2022
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)

#include <boost/hana/assert.hpp>
#include <boost/hana/optional.hpp>
namespace hana = boost::hana;


static_assert(hana::just(1).value_or('x') == 1, "");
static_assert(hana::nothing.value_or('x') == 'x', "");

int main() { }
