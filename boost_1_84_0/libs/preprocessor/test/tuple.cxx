# /* **************************************************************************
#  *                                                                          *
#  *     (C) Copyright Paul Mensonides 2002.
#  *     Distributed under the Boost Software License, Version 1.0. (See
#  *     accompanying file LICENSE_1_0.txt or copy at
#  *     http://www.boost.org/LICENSE_1_0.txt)
#  *                                                                          *
#  ************************************************************************** */
#
# /* Revised by Edward Diener (2011,2014,2020) */
#
# /* See http://www.boost.org for most recent version. */
#
# include <boost/preprocessor/config/limits.hpp>
# include <boost/preprocessor/cat.hpp>
# include <boost/preprocessor/control/if.hpp>
# include <boost/preprocessor/tuple.hpp>
# include <boost/preprocessor/array/elem.hpp>
# include <boost/preprocessor/array/size.hpp>
# include <boost/preprocessor/list/at.hpp>
# include <boost/preprocessor/list/size.hpp>
# include <boost/preprocessor/seq/elem.hpp>
# include <boost/preprocessor/seq/size.hpp>
# include <boost/preprocessor/facilities/is_empty.hpp>
# include <boost/preprocessor/variadic/size.hpp>
# include <boost/preprocessor/variadic/elem.hpp>
# include <boost/preprocessor/variadic/has_opt.hpp>
# include <libs/preprocessor/test/test.h>
# include <libs/preprocessor/test/tuple_elem_bug_test.cxx>

# define TUPLE (0, 1, 2, 3, 4, 5)
# define TUPLE_NONE ()
# define TUPLE_LARGE (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32)
# define TUPLE_VERY_LARGE (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63)
# define T2 (+3, /2, +6)

#if BOOST_PP_LIMIT_TUPLE > 64

# define TUPLE_LARGE_128 \
    ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, \
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103 )
    
# define TUPLE_VERY_LARGE_128 \
    ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, \
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127 )
    
#endif

#if BOOST_PP_LIMIT_TUPLE > 128

# define TUPLE_LARGE_256 \
    ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, \
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, \
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141 )

# define TUPLE_VERY_LARGE_256 \
    ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, \
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, \
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, \
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255 )

# define TUPLE_VERY_LARGE_255 \
    ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, \
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, \
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, \
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254 )

#endif

# define CALC(x) BOOST_PP_TUPLE_ELEM(0, x) BOOST_PP_TUPLE_ELEM(1, x) BOOST_PP_TUPLE_ELEM(2, x)
# define TEST_EAT BOOST_PP_TUPLE_EAT()(1, 2) 4
# define TEST_EAT_NONE BOOST_PP_TUPLE_EAT()() 17
# define TEST_EAT_LARGE BOOST_PP_TUPLE_EAT()(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32) 6
# define TEST_EAT_VERY_LARGE BOOST_PP_TUPLE_EAT()(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63) 8

#if BOOST_PP_LIMIT_TUPLE > 64

# define TEST_EAT_LARGE_128 \
    BOOST_PP_TUPLE_EAT() \
    ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, \
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99 ) \
    113
    
# define TEST_EAT_VERY_LARGE_128 \
    BOOST_PP_TUPLE_EAT() \
    ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, \
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127 ) \
    123

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

# define TEST_EAT_LARGE_256 \
    BOOST_PP_TUPLE_EAT() \
    ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, \
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, \
    128, 129, 130, 131, 132, 133, 134 ) \
    75
    
# define TEST_EAT_VERY_LARGE_256 \
    BOOST_PP_TUPLE_EAT() \
    ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, \
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, \
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, \
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255 ) \
    224

#endif

// elem

BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_TUPLE_ELEM(1, 0, TUPLE_NONE)) == 1 END
BEGIN BOOST_PP_TUPLE_ELEM(6, 3, TUPLE) == 3 END
BEGIN BOOST_PP_TUPLE_ELEM(6, 5, TUPLE) == 5 END
BEGIN BOOST_PP_TUPLE_ELEM(33, 15, TUPLE_LARGE) == 15 END
BEGIN BOOST_PP_TUPLE_ELEM(33, 27, TUPLE_LARGE) == 27 END
BEGIN BOOST_PP_TUPLE_ELEM(33, 32, TUPLE_LARGE) == 32 END
BEGIN BOOST_PP_TUPLE_ELEM(64, 22, TUPLE_VERY_LARGE) == 22 END
BEGIN BOOST_PP_TUPLE_ELEM(64, 47, TUPLE_VERY_LARGE) == 47 END
BEGIN BOOST_PP_TUPLE_ELEM(64, 63, TUPLE_VERY_LARGE) == 63 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_ELEM(104, 73, TUPLE_LARGE_128) == 73 END
BEGIN BOOST_PP_TUPLE_ELEM(104, 89, TUPLE_LARGE_128) == 89 END
BEGIN BOOST_PP_TUPLE_ELEM(104, 101, TUPLE_LARGE_128) == 101 END
BEGIN BOOST_PP_TUPLE_ELEM(128, 95, TUPLE_VERY_LARGE_128) == 95 END
BEGIN BOOST_PP_TUPLE_ELEM(128, 110, TUPLE_VERY_LARGE_128) == 110 END
BEGIN BOOST_PP_TUPLE_ELEM(128, 126, TUPLE_VERY_LARGE_128) == 126 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_ELEM(142, 83, TUPLE_LARGE_256) == 83 END
BEGIN BOOST_PP_TUPLE_ELEM(142, 131, TUPLE_LARGE_256) == 131 END
BEGIN BOOST_PP_TUPLE_ELEM(142, 140, TUPLE_LARGE_256) == 140 END
BEGIN BOOST_PP_TUPLE_ELEM(256, 174, TUPLE_VERY_LARGE_256) == 174 END
BEGIN BOOST_PP_TUPLE_ELEM(256, 226, TUPLE_VERY_LARGE_256) == 226 END
BEGIN BOOST_PP_TUPLE_ELEM(256, 253, TUPLE_VERY_LARGE_256) == 253 END

#endif

// reverse

BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_TUPLE_ELEM(1, 0, BOOST_PP_TUPLE_REVERSE(1,TUPLE_NONE))) == 1 END
BEGIN BOOST_PP_TUPLE_ELEM(6, 2, BOOST_PP_TUPLE_REVERSE(6,TUPLE)) == 3 END
BEGIN BOOST_PP_TUPLE_ELEM(33, 27, BOOST_PP_TUPLE_REVERSE(33,TUPLE_LARGE)) == 5 END
BEGIN BOOST_PP_TUPLE_ELEM(64, 43, BOOST_PP_TUPLE_REVERSE(64,TUPLE_VERY_LARGE)) == 20 END
BEGIN CALC(T2) == 7 END
BEGIN CALC(BOOST_PP_TUPLE_REVERSE(3, T2)) == 6 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_ELEM(104, 83, BOOST_PP_TUPLE_REVERSE(104,TUPLE_LARGE_128)) == 20 END
BEGIN BOOST_PP_TUPLE_ELEM(128, 119, BOOST_PP_TUPLE_REVERSE(128,TUPLE_VERY_LARGE_128)) == 8 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_ELEM(142, 56, BOOST_PP_TUPLE_REVERSE(142,TUPLE_LARGE_256)) == 85 END
BEGIN BOOST_PP_TUPLE_ELEM(256, 212, BOOST_PP_TUPLE_REVERSE(256,TUPLE_VERY_LARGE_256)) == 43 END

#endif

// eat

BEGIN TEST_EAT == 4 END
BEGIN TEST_EAT_NONE == 17 END
BEGIN TEST_EAT_LARGE == 6 END
BEGIN TEST_EAT_VERY_LARGE == 8 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN TEST_EAT_LARGE_128 == 113 END
BEGIN TEST_EAT_VERY_LARGE_128 == 123 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN TEST_EAT_LARGE_256 == 75 END
BEGIN TEST_EAT_VERY_LARGE_256 == 224 END

#endif

// to_array

BEGIN BOOST_PP_ARRAY_ELEM(3,BOOST_PP_TUPLE_TO_ARRAY(6,TUPLE)) == 3 END
BEGIN BOOST_PP_ARRAY_ELEM(29,BOOST_PP_TUPLE_TO_ARRAY(33,TUPLE_LARGE)) == 29 END
BEGIN BOOST_PP_ARRAY_ELEM(61,BOOST_PP_TUPLE_TO_ARRAY(64,TUPLE_VERY_LARGE)) == 61 END
BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_ARRAY_ELEM(0,BOOST_PP_TUPLE_TO_ARRAY(1,TUPLE_NONE))) == 1 END
BEGIN BOOST_PP_ARRAY_SIZE(BOOST_PP_TUPLE_TO_ARRAY(1,TUPLE_NONE)) == 1 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_ARRAY_ELEM(65,BOOST_PP_TUPLE_TO_ARRAY(104,TUPLE_LARGE_128)) == 65 END
BEGIN BOOST_PP_ARRAY_ELEM(117,BOOST_PP_TUPLE_TO_ARRAY(128,TUPLE_VERY_LARGE_128)) == 117 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_ARRAY_ELEM(131,BOOST_PP_TUPLE_TO_ARRAY(142,TUPLE_LARGE_256)) == 131 END
BEGIN BOOST_PP_ARRAY_ELEM(197,BOOST_PP_TUPLE_TO_ARRAY(256,TUPLE_VERY_LARGE_256)) == 197 END

#endif

// to_list

BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(6,TUPLE), 2) == 2 END
BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(33,TUPLE_LARGE), 19) == 19 END
BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(64,TUPLE_VERY_LARGE), 62) == 62 END
BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(1,TUPLE_NONE), 0)) == 1 END
BEGIN BOOST_PP_LIST_SIZE(BOOST_PP_TUPLE_TO_LIST(1,TUPLE_NONE)) == 1 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(104,TUPLE_LARGE_128), 88) == 88 END
BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(128,TUPLE_VERY_LARGE_128), 113) == 113 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(142,TUPLE_LARGE_256), 137) == 137 END
BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(256,TUPLE_VERY_LARGE_256), 235) == 235 END

#endif

// to_seq

BEGIN BOOST_PP_SEQ_ELEM(4,BOOST_PP_TUPLE_TO_SEQ(6,TUPLE)) == 4 END
BEGIN BOOST_PP_SEQ_ELEM(31,BOOST_PP_TUPLE_TO_SEQ(33,TUPLE_LARGE)) == 31 END
BEGIN BOOST_PP_SEQ_ELEM(55,BOOST_PP_TUPLE_TO_SEQ(64,TUPLE_VERY_LARGE)) == 55 END
BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_SEQ_ELEM(0,BOOST_PP_TUPLE_TO_SEQ(1,TUPLE_NONE))) == 1 END
BEGIN BOOST_PP_SEQ_SIZE(BOOST_PP_TUPLE_TO_SEQ(1,TUPLE_NONE)) == 1 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_SEQ_ELEM(97,BOOST_PP_TUPLE_TO_SEQ(104,TUPLE_LARGE_128)) == 97 END
BEGIN BOOST_PP_SEQ_ELEM(123,BOOST_PP_TUPLE_TO_SEQ(128,TUPLE_VERY_LARGE_128)) == 123 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_SEQ_ELEM(53,BOOST_PP_TUPLE_TO_SEQ(142,TUPLE_LARGE_256)) == 53 END
BEGIN BOOST_PP_SEQ_ELEM(181,BOOST_PP_TUPLE_TO_SEQ(256,TUPLE_VERY_LARGE_256)) == 181 END

#endif

// elem

BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_TUPLE_ELEM(0, TUPLE_NONE)) == 1 END
BEGIN BOOST_PP_TUPLE_ELEM(3, TUPLE) == 3 END
BEGIN BOOST_PP_TUPLE_ELEM(5, TUPLE) == 5 END
BEGIN BOOST_PP_TUPLE_ELEM(15, TUPLE_LARGE) == 15 END
BEGIN BOOST_PP_TUPLE_ELEM(27, TUPLE_LARGE) == 27 END
BEGIN BOOST_PP_TUPLE_ELEM(32, TUPLE_LARGE) == 32 END
BEGIN BOOST_PP_TUPLE_ELEM(22, TUPLE_VERY_LARGE) == 22 END
BEGIN BOOST_PP_TUPLE_ELEM(47, TUPLE_VERY_LARGE) == 47 END
BEGIN BOOST_PP_TUPLE_ELEM(63, TUPLE_VERY_LARGE) == 63 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_ELEM(73, TUPLE_LARGE_128) == 73 END
BEGIN BOOST_PP_TUPLE_ELEM(89, TUPLE_LARGE_128) == 89 END
BEGIN BOOST_PP_TUPLE_ELEM(101, TUPLE_LARGE_128) == 101 END
BEGIN BOOST_PP_TUPLE_ELEM(95, TUPLE_VERY_LARGE_128) == 95 END
BEGIN BOOST_PP_TUPLE_ELEM(110, TUPLE_VERY_LARGE_128) == 110 END
BEGIN BOOST_PP_TUPLE_ELEM(126, TUPLE_VERY_LARGE_128) == 126 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_ELEM(83, TUPLE_LARGE_256) == 83 END
BEGIN BOOST_PP_TUPLE_ELEM(131, TUPLE_LARGE_256) == 131 END
BEGIN BOOST_PP_TUPLE_ELEM(140, TUPLE_LARGE_256) == 140 END
BEGIN BOOST_PP_TUPLE_ELEM(174, TUPLE_VERY_LARGE_256) == 174 END
BEGIN BOOST_PP_TUPLE_ELEM(226, TUPLE_VERY_LARGE_256) == 226 END
BEGIN BOOST_PP_TUPLE_ELEM(253, TUPLE_VERY_LARGE_256) == 253 END

#endif

// reverse

BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_REVERSE(TUPLE_NONE))) == 1 END
BEGIN BOOST_PP_TUPLE_ELEM(2, BOOST_PP_TUPLE_REVERSE(TUPLE)) == 3 END
BEGIN BOOST_PP_TUPLE_ELEM(27, BOOST_PP_TUPLE_REVERSE(TUPLE_LARGE)) == 5 END
BEGIN BOOST_PP_TUPLE_ELEM(43, BOOST_PP_TUPLE_REVERSE(TUPLE_VERY_LARGE)) == 20 END
BEGIN CALC(BOOST_PP_TUPLE_REVERSE(T2)) == 6 END
BEGIN BOOST_PP_VARIADIC_ELEM(2,BOOST_PP_TUPLE_ENUM(BOOST_PP_TUPLE_REVERSE(TUPLE))) == 3 END
BEGIN BOOST_PP_VARIADIC_ELEM(27,BOOST_PP_TUPLE_ENUM(BOOST_PP_TUPLE_REVERSE(TUPLE_LARGE))) == 5 END
BEGIN BOOST_PP_VARIADIC_ELEM(45,BOOST_PP_TUPLE_ENUM(BOOST_PP_TUPLE_REVERSE(TUPLE_VERY_LARGE))) == 18 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_ELEM(83, BOOST_PP_TUPLE_REVERSE(TUPLE_LARGE_128)) == 20 END
BEGIN BOOST_PP_TUPLE_ELEM(119, BOOST_PP_TUPLE_REVERSE(TUPLE_VERY_LARGE_128)) == 8 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_ELEM(56, BOOST_PP_TUPLE_REVERSE(TUPLE_LARGE_256)) == 85 END
BEGIN BOOST_PP_TUPLE_ELEM(212, BOOST_PP_TUPLE_REVERSE(TUPLE_VERY_LARGE_256)) == 43 END

#endif

// enum

# if BOOST_PP_VARIADIC_HAS_OPT()

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_ENUM(TUPLE_NONE)) == 0 END

# else

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_ENUM(TUPLE_NONE)) == 1 END

# endif

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_ENUM(TUPLE)) == 6 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_ENUM(TUPLE_LARGE)) == 33 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_ENUM(TUPLE_VERY_LARGE)) == 64 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_ENUM(TUPLE_LARGE_128)) == 104 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_ENUM(TUPLE_VERY_LARGE_128)) == 128 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_ENUM(TUPLE_LARGE_256)) == 142 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_ENUM(TUPLE_VERY_LARGE_256)) == 256 END

#endif

// insert

BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_INSERT(TUPLE_NONE,0,40)) == 40 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_INSERT(TUPLE_NONE,0,40)) == 2 END

BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_INSERT(TUPLE,2,40)) == 0 END
BEGIN BOOST_PP_TUPLE_ELEM(1, BOOST_PP_TUPLE_INSERT(TUPLE,1,40)) == 40 END
BEGIN BOOST_PP_TUPLE_ELEM(2, BOOST_PP_TUPLE_INSERT(TUPLE,1,40)) == 1 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_INSERT(TUPLE,1,40)) == 7 END

BEGIN BOOST_PP_TUPLE_ELEM(8, BOOST_PP_TUPLE_INSERT(TUPLE_LARGE,22,1000)) == 8 END
BEGIN BOOST_PP_TUPLE_ELEM(22, BOOST_PP_TUPLE_INSERT(TUPLE_LARGE,22,1000)) == 1000 END
BEGIN BOOST_PP_TUPLE_ELEM(26, BOOST_PP_TUPLE_INSERT(TUPLE_LARGE,22,1000)) == 25 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_INSERT(TUPLE_LARGE,22,1000)) == 34 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_ELEM(69, BOOST_PP_TUPLE_INSERT(TUPLE_LARGE_128,76,1000)) == 69 END
BEGIN BOOST_PP_TUPLE_ELEM(101, BOOST_PP_TUPLE_INSERT(TUPLE_LARGE_128,101,1000)) == 1000 END
BEGIN BOOST_PP_TUPLE_ELEM(98, BOOST_PP_TUPLE_INSERT(TUPLE_LARGE_128,96,1000)) == 97 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_INSERT(TUPLE_LARGE_128,84,1000)) == 105 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_ELEM(73, BOOST_PP_TUPLE_INSERT(TUPLE_LARGE_256,134,1000)) == 73 END
BEGIN BOOST_PP_TUPLE_ELEM(141, BOOST_PP_TUPLE_INSERT(TUPLE_LARGE_256,141,1000)) == 1000 END
BEGIN BOOST_PP_TUPLE_ELEM(133, BOOST_PP_TUPLE_INSERT(TUPLE_LARGE_256,39,1000)) == 132 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_INSERT(TUPLE_LARGE_256,121,1000)) == 143 END
BEGIN BOOST_PP_TUPLE_ELEM(227,BOOST_PP_TUPLE_INSERT(TUPLE_VERY_LARGE_255,212,1000)) == 226 END
BEGIN BOOST_PP_TUPLE_ELEM(212,BOOST_PP_TUPLE_INSERT(TUPLE_VERY_LARGE_255,212,1000)) == 1000 END

#endif

// pop_back

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_BACK(TUPLE)) == 5 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_BACK(TUPLE_LARGE)) == 32 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_BACK(TUPLE_VERY_LARGE)) == 63 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_BACK(TUPLE_LARGE_128)) == 103 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_BACK(TUPLE_VERY_LARGE_128)) == 127 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_BACK(TUPLE_LARGE_256)) == 141 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_BACK(TUPLE_VERY_LARGE_256)) == 255 END

#endif

// pop_front

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_FRONT(TUPLE)) == 5 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_FRONT(TUPLE_LARGE)) == 32 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_FRONT(TUPLE_VERY_LARGE)) == 63 END

BEGIN BOOST_PP_TUPLE_ELEM(1, BOOST_PP_TUPLE_POP_FRONT(TUPLE)) == 2 END
BEGIN BOOST_PP_TUPLE_ELEM(31, BOOST_PP_TUPLE_POP_FRONT(TUPLE_LARGE)) == 32 END
BEGIN BOOST_PP_TUPLE_ELEM(55, BOOST_PP_TUPLE_POP_FRONT(TUPLE_VERY_LARGE)) == 56 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_FRONT(TUPLE_LARGE_128)) == 103 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_FRONT(TUPLE_VERY_LARGE_128)) == 127 END
BEGIN BOOST_PP_TUPLE_ELEM(84, BOOST_PP_TUPLE_POP_FRONT(TUPLE_LARGE_128)) == 85 END
BEGIN BOOST_PP_TUPLE_ELEM(117, BOOST_PP_TUPLE_POP_FRONT(TUPLE_VERY_LARGE_128)) == 118 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_FRONT(TUPLE_LARGE_256)) == 141 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_POP_FRONT(TUPLE_VERY_LARGE_256)) == 255 END
BEGIN BOOST_PP_TUPLE_ELEM(129, BOOST_PP_TUPLE_POP_FRONT(TUPLE_LARGE_256)) == 130 END
BEGIN BOOST_PP_TUPLE_ELEM(248, BOOST_PP_TUPLE_POP_FRONT(TUPLE_VERY_LARGE_256)) == 249 END

#endif

// push_back

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_BACK(TUPLE_NONE, 1)) == 2 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_BACK(TUPLE, 6)) == 7 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_BACK(TUPLE_LARGE, 33)) == 34 END
BEGIN BOOST_PP_TUPLE_ELEM(1, BOOST_PP_TUPLE_PUSH_BACK(TUPLE_NONE, 1)) == 1 END
BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_PUSH_BACK(TUPLE, 6)) == 0 END
BEGIN BOOST_PP_TUPLE_ELEM(33, BOOST_PP_TUPLE_PUSH_BACK(TUPLE_LARGE, 33)) == 33 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_BACK(TUPLE_LARGE_128, 66)) == 105 END
BEGIN BOOST_PP_TUPLE_ELEM(104, BOOST_PP_TUPLE_PUSH_BACK(TUPLE_LARGE_128, 101)) == 101 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_BACK(TUPLE_LARGE_256, 192)) == 143 END
BEGIN BOOST_PP_TUPLE_ELEM(142, BOOST_PP_TUPLE_PUSH_BACK(TUPLE_LARGE_256, 77)) == 77 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_BACK(TUPLE_VERY_LARGE_255, 255)) == 256 END

#endif

// push_front

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_FRONT(TUPLE_NONE, 55)) == 2 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_FRONT(TUPLE, 555)) == 7 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_FRONT(TUPLE_LARGE, 666)) == 34 END
BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_PUSH_FRONT(TUPLE_NONE, 55)) == 55 END
BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_PUSH_FRONT(TUPLE, 555)) == 555 END
BEGIN BOOST_PP_TUPLE_ELEM(33, BOOST_PP_TUPLE_PUSH_FRONT(TUPLE_LARGE, 33)) == 32 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_FRONT(TUPLE_LARGE_128, 666)) == 105 END
BEGIN BOOST_PP_TUPLE_ELEM(103, BOOST_PP_TUPLE_PUSH_FRONT(TUPLE_LARGE_128, 29)) == 102 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_FRONT(TUPLE_LARGE_256, 333)) == 143 END
BEGIN BOOST_PP_TUPLE_ELEM(136, BOOST_PP_TUPLE_PUSH_FRONT(TUPLE_LARGE_256, 47)) == 135 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_PUSH_FRONT(TUPLE_VERY_LARGE_255, 4)) == 256 END

#endif

// rem

#if BOOST_PP_VARIADICS_MSVC && _MSC_VER <= 1400

# if BOOST_PP_VARIADIC_HAS_OPT()

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM(1)()) == 0 END

# else

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM(1)()) == 1 END

# endif

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM(33)(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32)) == 33 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM(64)(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63)) == 64 END

#else

#if BOOST_PP_VARIADICS_MSVC
BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_TUPLE_REM_CAT() TUPLE_NONE) == 1 END
#else
BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_TUPLE_REM() TUPLE_NONE) == 1 END
#endif

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM()(0, 1, 2, 3, 4, 5)) == 6 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM()(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32)) == 33 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM()(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63)) == 64 END

#endif

// rem_ctor

# if BOOST_PP_VARIADIC_HAS_OPT()

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM_CTOR(TUPLE_NONE)) == 0 END

# else

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM_CTOR(TUPLE_NONE)) == 1 END

# endif

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM_CTOR(TUPLE)) == 6 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM_CTOR(TUPLE_LARGE)) == 33 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM_CTOR(TUPLE_VERY_LARGE)) == 64 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM_CTOR(TUPLE_LARGE_128)) == 104 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM_CTOR(TUPLE_VERY_LARGE_128)) == 128 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM_CTOR(TUPLE_LARGE_256)) == 142 END
BEGIN BOOST_PP_VARIADIC_SIZE(BOOST_PP_TUPLE_REM_CTOR(TUPLE_VERY_LARGE_256)) == 256 END

#endif

// remove

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REMOVE(TUPLE, 1)) == 5 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REMOVE(TUPLE_LARGE, 17)) == 32 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REMOVE(TUPLE_VERY_LARGE, 27)) == 63 END
BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_REMOVE(TUPLE, 2)) == 0 END
BEGIN BOOST_PP_TUPLE_ELEM(29, BOOST_PP_TUPLE_REMOVE(TUPLE_LARGE, 25)) == 30 END
BEGIN BOOST_PP_TUPLE_ELEM(62, BOOST_PP_TUPLE_REMOVE(TUPLE_VERY_LARGE, 48)) == 63 END

#if BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_STRICT()

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REMOVE(TUPLE_LARGE, 0)) == 32 END
BEGIN BOOST_PP_TUPLE_ELEM(25, BOOST_PP_TUPLE_REMOVE(TUPLE_VERY_LARGE, 0)) == 26 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REMOVE(TUPLE_LARGE_128, 100)) == 103 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REMOVE(TUPLE_VERY_LARGE_128, 123)) == 127 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REMOVE(TUPLE_VERY_LARGE_128, 0)) == 127 END
BEGIN BOOST_PP_TUPLE_ELEM(102, BOOST_PP_TUPLE_REMOVE(TUPLE_LARGE_128, 97)) == 103 END
BEGIN BOOST_PP_TUPLE_ELEM(76, BOOST_PP_TUPLE_REMOVE(TUPLE_LARGE_128, 0)) == 77 END
BEGIN BOOST_PP_TUPLE_ELEM(119, BOOST_PP_TUPLE_REMOVE(TUPLE_VERY_LARGE_128, 115)) == 120 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REMOVE(TUPLE_LARGE_256, 133)) == 141 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REMOVE(TUPLE_LARGE_256, 0)) == 141 END
BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REMOVE(TUPLE_VERY_LARGE_256, 241)) == 255 END
BEGIN BOOST_PP_TUPLE_ELEM(140, BOOST_PP_TUPLE_REMOVE(TUPLE_LARGE_256, 138)) == 141 END
BEGIN BOOST_PP_TUPLE_ELEM(181, BOOST_PP_TUPLE_REMOVE(TUPLE_VERY_LARGE_256, 166)) == 182 END
BEGIN BOOST_PP_TUPLE_ELEM(236, BOOST_PP_TUPLE_REMOVE(TUPLE_VERY_LARGE_256, 0)) == 237 END

#endif

// replace

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE, 27, 1000)) == 64 END
BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_REPLACE(TUPLE_NONE, 0, 71)) == 71 END
BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_REPLACE(TUPLE, 1, 44)) == 0 END
BEGIN BOOST_PP_TUPLE_ELEM(29, BOOST_PP_TUPLE_REPLACE(TUPLE_LARGE, 29, 999)) == 999 END
BEGIN BOOST_PP_TUPLE_ELEM(38, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE, 37, 1)) == 38 END
BEGIN BOOST_PP_TUPLE_ELEM(28, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE, 28, 1)) == 1 END

#if BOOST_PP_CONFIG_FLAGS() & BOOST_PP_CONFIG_STRICT()

BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE, 0, 64)) == 64 END
BEGIN BOOST_PP_TUPLE_ELEM(17, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE, 0, 256)) == 17 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE_128, 93, 1000)) == 128 END
BEGIN BOOST_PP_TUPLE_ELEM(89, BOOST_PP_TUPLE_REPLACE(TUPLE_LARGE_128, 89, 111)) == 111 END
BEGIN BOOST_PP_TUPLE_ELEM(73, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE_128, 66, 1)) == 73 END
BEGIN BOOST_PP_TUPLE_ELEM(122, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE_128, 122, 1)) == 1 END
BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE_128, 0, 128)) == 128 END
BEGIN BOOST_PP_TUPLE_ELEM(95, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE_128, 0, 128)) == 95 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_SIZE(BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE_256, 217, 1000)) == 256 END
BEGIN BOOST_PP_TUPLE_ELEM(136, BOOST_PP_TUPLE_REPLACE(TUPLE_LARGE_256, 136, 999)) == 999 END
BEGIN BOOST_PP_TUPLE_ELEM(192, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE_256, 185, 1)) == 192 END
BEGIN BOOST_PP_TUPLE_ELEM(237, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE_256, 237, 1)) == 1 END
BEGIN BOOST_PP_TUPLE_ELEM(0, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE_256, 0, 256)) == 256 END
BEGIN BOOST_PP_TUPLE_ELEM(167, BOOST_PP_TUPLE_REPLACE(TUPLE_VERY_LARGE_256, 0, 256)) == 167 END

#endif

// size

BEGIN BOOST_PP_TUPLE_SIZE(TUPLE_NONE) == 1 END
BEGIN BOOST_PP_TUPLE_SIZE(TUPLE) == 6 END
BEGIN BOOST_PP_TUPLE_SIZE(TUPLE_LARGE) == 33 END
BEGIN BOOST_PP_TUPLE_SIZE(TUPLE_VERY_LARGE) == 64 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_TUPLE_SIZE(TUPLE_LARGE_128) == 104 END
BEGIN BOOST_PP_TUPLE_SIZE(TUPLE_VERY_LARGE_128) == 128 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_TUPLE_SIZE(TUPLE_LARGE_256) == 142 END
BEGIN BOOST_PP_TUPLE_SIZE(TUPLE_VERY_LARGE_256) == 256 END

#endif

// to_array

BEGIN BOOST_PP_ARRAY_ELEM(3,BOOST_PP_TUPLE_TO_ARRAY(TUPLE)) == 3 END
BEGIN BOOST_PP_ARRAY_ELEM(29,BOOST_PP_TUPLE_TO_ARRAY(TUPLE_LARGE)) == 29 END
BEGIN BOOST_PP_ARRAY_ELEM(61,BOOST_PP_TUPLE_TO_ARRAY(TUPLE_VERY_LARGE)) == 61 END
BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_ARRAY_ELEM(0,BOOST_PP_TUPLE_TO_ARRAY(TUPLE_NONE))) == 1 END
BEGIN BOOST_PP_ARRAY_SIZE(BOOST_PP_TUPLE_TO_ARRAY(TUPLE_NONE)) == 1 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_ARRAY_ELEM(65,BOOST_PP_TUPLE_TO_ARRAY(TUPLE_LARGE_128)) == 65 END
BEGIN BOOST_PP_ARRAY_ELEM(117,BOOST_PP_TUPLE_TO_ARRAY(TUPLE_VERY_LARGE_128)) == 117 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_ARRAY_ELEM(131,BOOST_PP_TUPLE_TO_ARRAY(TUPLE_LARGE_256)) == 131 END
BEGIN BOOST_PP_ARRAY_ELEM(197,BOOST_PP_TUPLE_TO_ARRAY(TUPLE_VERY_LARGE_256)) == 197 END

#endif

// to_list

BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(TUPLE), 2) == 2 END
BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(TUPLE_LARGE), 19) == 19 END
BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(TUPLE_VERY_LARGE), 62) == 62 END
BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(TUPLE_NONE), 0)) == 1 END
BEGIN BOOST_PP_LIST_SIZE(BOOST_PP_TUPLE_TO_LIST(TUPLE_NONE)) == 1 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(TUPLE_LARGE_128), 88) == 88 END
BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(TUPLE_VERY_LARGE_128), 113) == 113 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(TUPLE_LARGE_256), 137) == 137 END
BEGIN BOOST_PP_LIST_AT(BOOST_PP_TUPLE_TO_LIST(TUPLE_VERY_LARGE_256), 235) == 235 END

#endif

// to_seq

BEGIN BOOST_PP_SEQ_ELEM(4,BOOST_PP_TUPLE_TO_SEQ(TUPLE)) == 4 END
BEGIN BOOST_PP_SEQ_ELEM(31,BOOST_PP_TUPLE_TO_SEQ(TUPLE_LARGE)) == 31 END
BEGIN BOOST_PP_SEQ_ELEM(55,BOOST_PP_TUPLE_TO_SEQ(TUPLE_VERY_LARGE)) == 55 END
BEGIN BOOST_PP_IS_EMPTY(BOOST_PP_SEQ_ELEM(0,BOOST_PP_TUPLE_TO_SEQ(TUPLE_NONE))) == 1 END
BEGIN BOOST_PP_SEQ_SIZE(BOOST_PP_TUPLE_TO_SEQ(TUPLE_NONE)) == 1 END

#if BOOST_PP_LIMIT_TUPLE > 64

BEGIN BOOST_PP_SEQ_ELEM(97,BOOST_PP_TUPLE_TO_SEQ(TUPLE_LARGE_128)) == 97 END
BEGIN BOOST_PP_SEQ_ELEM(123,BOOST_PP_TUPLE_TO_SEQ(TUPLE_VERY_LARGE_128)) == 123 END

#endif

#if BOOST_PP_LIMIT_TUPLE > 128

BEGIN BOOST_PP_SEQ_ELEM(53,BOOST_PP_TUPLE_TO_SEQ(TUPLE_LARGE_256)) == 53 END
BEGIN BOOST_PP_SEQ_ELEM(181,BOOST_PP_TUPLE_TO_SEQ(TUPLE_VERY_LARGE_256)) == 181 END

#endif
