// Copyright (c) 2008-2009 Ben Hanson
// Copyright (c) 2008-2011 Hartmut Kaiser
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file licence_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Auto-generated by boost::lexer, do not edit

#if !defined(BOOST_SPIRIT_LEXER_NEXT_TOKEN_CONJURE_STATIC_JUL_25_2011_07_03_08)
#define BOOST_SPIRIT_LEXER_NEXT_TOKEN_CONJURE_STATIC_JUL_25_2011_07_03_08

#include <boost/spirit/home/support/detail/lexer/char_traits.hpp>

////////////////////////////////////////////////////////////////////////////////
// the generated table of state names and the tokenizer have to be
// defined in the boost::spirit::lex::lexertl::static_ namespace
namespace boost { namespace spirit { namespace lex { namespace lexertl { namespace static_ {

////////////////////////////////////////////////////////////////////////////////
// this table defines the names of the lexer states
char const* const lexer_state_names_conjure_static[1] = 
{
    "INITIAL"
};

////////////////////////////////////////////////////////////////////////////////
// this variable defines the number of lexer states
std::size_t const lexer_state_count_conjure_static = 1;

////////////////////////////////////////////////////////////////////////////////
// this function returns the next matched token
template<typename Iterator>
std::size_t next_token_conjure_static (std::size_t& /*start_state_*/, bool& /*bol_*/, 
    Iterator &start_token_, Iterator const& end_, std::size_t& unique_id_)
{
    enum {end_state_index, id_index, unique_id_index, state_index, bol_index,
        eol_index, dead_state_index, dfa_offset};

    static std::size_t const npos = static_cast<std::size_t>(~0);
    static std::size_t const lookup_[256] = {
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 7, 7, 41, 41, 7, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        7, 8, 41, 41, 41, 41, 9, 41,
        10, 11, 12, 13, 14, 15, 41, 16,
        17, 17, 17, 17, 17, 17, 17, 17,
        17, 17, 41, 19, 20, 21, 22, 41,
        41, 18, 18, 18, 18, 18, 18, 18,
        18, 18, 18, 18, 18, 18, 18, 18,
        18, 18, 18, 18, 18, 18, 18, 18,
        18, 18, 18, 41, 41, 41, 41, 18,
        41, 23, 18, 18, 24, 25, 26, 18,
        27, 28, 18, 18, 29, 18, 30, 31,
        18, 18, 32, 33, 34, 35, 36, 37,
        18, 18, 18, 38, 39, 40, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41,
        41, 41, 41, 41, 41, 41, 41, 41 };
    static std::size_t const dfa_alphabet_ = 42;
    static std::size_t const dfa_[2604] = {
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 27, 13, 11, 20, 21, 18, 16,
        24, 17, 19, 2, 26, 25, 14, 12,
        15, 26, 26, 7, 4, 26, 6, 26,
        26, 26, 9, 26, 3, 26, 5, 8,
        22, 10, 23, 0, 1, 35, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 2, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 32,
        28, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 26,
        26, 0, 0, 0, 0, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 28, 26,
        26, 26, 26, 26, 0, 0, 0, 0,
        1, 32, 28, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 26, 26, 0, 0, 0, 0, 29,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 0, 0,
        0, 0, 1, 32, 28, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 26, 26, 0, 0, 0,
        0, 26, 26, 26, 26, 26, 26, 26,
        26, 30, 26, 26, 26, 26, 26, 26,
        0, 0, 0, 0, 1, 32, 28, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 26, 26, 0,
        0, 0, 0, 26, 26, 26, 32, 26,
        26, 26, 31, 26, 26, 26, 26, 26,
        26, 26, 0, 0, 0, 0, 1, 32,
        28, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 26,
        26, 0, 0, 0, 0, 26, 26, 26,
        26, 26, 26, 33, 26, 26, 26, 26,
        26, 26, 26, 26, 0, 0, 0, 0,
        1, 32, 28, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 26, 26, 0, 0, 0, 0, 26,
        26, 26, 26, 34, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 0, 0,
        0, 0, 1, 32, 28, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 26, 26, 0, 0, 0,
        0, 26, 26, 35, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 36, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 37,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        1, 61, 26, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 38, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 262177, 20, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 39,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 131091, 12, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 40, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 131093,
        14, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 41, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        1, 393241, 16, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 393242, 17, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 131099, 18, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 131100,
        19, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 42, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        1, 40, 21, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 41, 22, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 123, 23, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 125,
        24, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        1, 44, 25, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 59, 27, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 32, 28, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 26, 26, 0,
        0, 0, 0, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 0, 0, 0, 0, 1, 34,
        30, 0, 0, 0, 0, 27, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        1, 32, 28, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 26, 26, 0, 0, 0, 0, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 43, 26, 26, 0, 0,
        0, 0, 1, 32, 28, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 26, 26, 0, 0, 0,
        0, 26, 26, 26, 26, 26, 26, 44,
        26, 26, 26, 26, 26, 26, 26, 26,
        0, 0, 0, 0, 1, 32, 28, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 26, 26, 0,
        0, 0, 0, 26, 26, 26, 26, 26,
        45, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 0, 0, 0, 0, 1, 32,
        28, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 26,
        26, 0, 0, 0, 0, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        46, 26, 26, 26, 0, 0, 0, 0,
        1, 65538, 4, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 26, 26, 0, 0, 0, 0, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 0, 0,
        0, 0, 1, 32, 28, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 26, 26, 0, 0, 0,
        0, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 47, 26, 26, 26, 26,
        0, 0, 0, 0, 1, 32, 28, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 26, 26, 0,
        0, 0, 0, 26, 26, 26, 26, 26,
        48, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 0, 0, 0, 0, 1, 32,
        28, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 26,
        26, 0, 0, 0, 0, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        49, 26, 26, 26, 0, 0, 0, 0,
        1, 131084, 8, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 131085, 9, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 131089, 10, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 131090,
        11, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        1, 131092, 13, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 131094, 15, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 42, 42, 42, 42, 42,
        50, 42, 42, 42, 42, 42, 42, 42,
        42, 42, 42, 42, 42, 42, 42, 42,
        42, 42, 42, 42, 42, 42, 42, 42,
        42, 42, 42, 42, 42, 42, 1, 32,
        28, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 26,
        26, 0, 0, 0, 0, 26, 26, 51,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 0, 0, 0, 0,
        1, 32, 28, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 26, 26, 0, 0, 0, 0, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 43, 26, 26, 26, 26, 0, 0,
        0, 0, 1, 32, 28, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 26, 26, 0, 0, 0,
        0, 26, 52, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        0, 0, 0, 0, 1, 65537, 3, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 26, 26, 0,
        0, 0, 0, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 0, 0, 0, 0, 1, 32,
        28, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 26,
        26, 0, 0, 0, 0, 26, 26, 53,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 0, 0, 0, 0,
        1, 32, 28, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 26, 26, 0, 0, 0, 0, 26,
        26, 26, 26, 26, 26, 54, 26, 26,
        26, 26, 26, 26, 26, 26, 0, 0,
        0, 0, 1, 32, 28, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 26, 26, 0, 0, 0,
        0, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 55, 26, 26,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 56, 56, 56, 56, 56,
        50, 56, 56, 56, 57, 56, 56, 56,
        56, 56, 56, 56, 56, 56, 56, 56,
        56, 56, 56, 56, 56, 56, 56, 56,
        56, 56, 56, 56, 56, 56, 1, 36,
        1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 26,
        26, 0, 0, 0, 0, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 0, 0, 0, 0,
        1, 65536, 2, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 26, 26, 0, 0, 0, 0, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 0, 0,
        0, 0, 1, 65539, 5, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 26, 26, 0, 0, 0,
        0, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        0, 0, 0, 0, 1, 32, 28, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 26, 26, 0,
        0, 0, 0, 26, 26, 58, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 0, 0, 0, 0, 1, 32,
        28, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 26,
        26, 0, 0, 0, 0, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 59, 26,
        26, 26, 26, 26, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 56,
        56, 56, 56, 56, 60, 56, 56, 56,
        56, 56, 56, 56, 56, 56, 56, 56,
        56, 56, 56, 56, 56, 56, 56, 56,
        56, 56, 56, 56, 56, 56, 56, 56,
        56, 56, 1, 33, 29, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 65540, 6, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 26, 26, 0,
        0, 0, 0, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 0, 0, 0, 0, 1, 32,
        28, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 26,
        26, 0, 0, 0, 0, 26, 26, 26,
        26, 26, 26, 26, 61, 26, 26, 26,
        26, 26, 26, 26, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 56,
        56, 56, 56, 56, 60, 56, 56, 56,
        57, 56, 56, 56, 56, 56, 56, 56,
        56, 56, 56, 56, 56, 56, 56, 56,
        56, 56, 56, 56, 56, 56, 56, 56,
        56, 56, 1, 65541, 7, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 26, 26, 0, 0, 0,
        0, 26, 26, 26, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26,
        0, 0, 0, 0 };

    if (start_token_ == end_)
    {
        unique_id_ = npos;
        return 0;
    }

    std::size_t const* ptr_ = dfa_ + dfa_alphabet_;
    Iterator curr_ = start_token_;
    bool end_state_ = *ptr_ != 0;
    std::size_t id_ = *(ptr_ + id_index);
    std::size_t uid_ = *(ptr_ + unique_id_index);
    Iterator end_token_ = start_token_;

    while (curr_ != end_)
    {
        std::size_t const state_ =
            ptr_[lookup_[static_cast<unsigned char>(*curr_++)]];

        if (state_ == 0) break;

        ptr_ = &dfa_[state_ * dfa_alphabet_];

        if (*ptr_)
        {
            end_state_ = true;
            id_ = *(ptr_ + id_index);
            uid_ = *(ptr_ + unique_id_index);
            end_token_ = curr_;
        }
    }

    if (end_state_)
    {
        // return longest match
        start_token_ = end_token_;
    }
    else
    {
        id_ = npos;
        uid_ = npos;
    }

    unique_id_ = uid_;
    return id_;
}

////////////////////////////////////////////////////////////////////////////////
// this defines a generic accessors for the information above
struct lexer_conjure_static
{
    // version number and feature-set of compatible static lexer engine
    enum
    {
        static_version = 65536,
        supports_bol = false,
        supports_eol = false
    };

    // return the number of lexer states
    static std::size_t state_count()
    {
        return lexer_state_count_conjure_static; 
    }

    // return the name of the lexer state as given by 'idx'
    static char const* state_name(std::size_t idx)
    {
        return lexer_state_names_conjure_static[idx]; 
    }

    // return the next matched token
    template<typename Iterator>
    static std::size_t next(std::size_t &start_state_, bool& bol_
      , Iterator &start_token_, Iterator const& end_, std::size_t& unique_id_)
    {
        return next_token_conjure_static(start_state_, bol_, start_token_, end_, unique_id_);
    }
};

}}}}}  // namespace boost::spirit::lex::lexertl::static_

#endif
