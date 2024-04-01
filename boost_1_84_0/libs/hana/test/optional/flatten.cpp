// Copyright Louis Dionne 2013-2022
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)

#include <boost/hana/assert.hpp>
#include <boost/hana/equal.hpp>
#include <boost/hana/flatten.hpp>
#include <boost/hana/optional.hpp>

#include <laws/base.hpp>
namespace hana = boost::hana;
using hana::test::ct_eq;


int main() {
    BOOST_HANA_CONSTANT_CHECK(hana::equal(
        hana::flatten(hana::nothing),
        hana::nothing
    ));
    BOOST_HANA_CONSTANT_CHECK(hana::equal(
        hana::flatten(hana::just(hana::nothing)),
        hana::nothing
    ));
    BOOST_HANA_CONSTANT_CHECK(hana::equal(
        hana::flatten(hana::just(hana::just(ct_eq<4>{}))),
        hana::just(ct_eq<4>{})
    ));
}
