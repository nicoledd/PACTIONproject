# /* **************************************************************************
#  *                                                                          *
#  *     (C) Copyright Edward Diener 2019.
#  *     Distributed under the Boost Software License, Version 1.0. (See
#  *     accompanying file LICENSE_1_0.txt or copy at
#  *     http://www.boost.org/LICENSE_1_0.txt)
#  *                                                                          *
#  ************************************************************************** */
#
# /* See http://www.boost.org for most recent version. */
#
# include <libs/preprocessor/test/test.h>
# include <boost/preprocessor/variadic/has_opt.hpp>

# if BOOST_PP_VARIADIC_HAS_OPT()

# include <boost/preprocessor/facilities/empty.hpp>
# include <boost/preprocessor/facilities/va_opt.hpp>
# include <boost/preprocessor/variadic/elem.hpp>

#define DATA
#define OBJECT OBJECT2
#define OBJECT2
#define FUNC(x) FUNC2(x)
#define FUNC2(x)
#define FUNC_GEN() ()
#define FUNC_GEN2(x) ()
#define FUNC_GEN3() (&)
#define FUNC_GEN4(x) (y)
#define FUNC_GEN5() (y,z)
#define FUNC_GEN6() anything
#define FUNC_GEN7(x) anything
#define FUNC_GEN8(x,y) (1,2,3)
#define FUNC_GEN9(x,y,z) anything
#define FUNC_GEN10(x) (y) data
#define NAME &name
#define ATUPLE (atuple)
#define ATUPLE_PLUS (atuple) data

#define TRY_PASS_EMPTY_STD(...) BOOST_PP_VARIADIC_ELEM(0, __VA_OPT__ (0,) 1)
#define TRY_PASS_EMPTY_BOOST(...) BOOST_PP_VARIADIC_ELEM(0,BOOST_PP_VA_OPT((0),(1),__VA_ARGS__))

BEGIN TRY_PASS_EMPTY_STD(FUNC_GEN) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC_GEN) == 0 END
BEGIN TRY_PASS_EMPTY_STD(FUNC_GEN2) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC_GEN2) == 0 END
BEGIN TRY_PASS_EMPTY_STD(FUNC_GEN3) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC_GEN3) == 0 END
BEGIN TRY_PASS_EMPTY_STD(FUNC_GEN4) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC_GEN4) == 0 END
BEGIN TRY_PASS_EMPTY_STD(FUNC_GEN5) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC_GEN5) == 0 END
BEGIN TRY_PASS_EMPTY_STD(FUNC_GEN8) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC_GEN8) == 0 END
BEGIN TRY_PASS_EMPTY_STD(FUNC_GEN9) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC_GEN9) == 0 END
BEGIN TRY_PASS_EMPTY_STD(FUNC_GEN10) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC_GEN10) == 0 END
BEGIN TRY_PASS_EMPTY_STD(BOOST_PP_EMPTY()) == 1 END
BEGIN TRY_PASS_EMPTY_BOOST(BOOST_PP_EMPTY()) == 1 END
BEGIN TRY_PASS_EMPTY_STD(DATA BOOST_PP_EMPTY()) == 1 END
BEGIN TRY_PASS_EMPTY_BOOST(DATA BOOST_PP_EMPTY()) == 1 END
BEGIN TRY_PASS_EMPTY_STD() == 1 END
BEGIN TRY_PASS_EMPTY_BOOST() == 1 END
BEGIN TRY_PASS_EMPTY_STD(DATA) == 1 END
BEGIN TRY_PASS_EMPTY_BOOST(DATA) == 1 END
BEGIN TRY_PASS_EMPTY_STD(x BOOST_PP_EMPTY()) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(x BOOST_PP_EMPTY()) == 0 END
BEGIN TRY_PASS_EMPTY_STD(OBJECT BOOST_PP_EMPTY()) == 1 END
BEGIN TRY_PASS_EMPTY_BOOST(OBJECT BOOST_PP_EMPTY()) == 1 END
BEGIN TRY_PASS_EMPTY_STD(FUNC(z) BOOST_PP_EMPTY()) == 1 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC(z) BOOST_PP_EMPTY()) == 1 END
BEGIN TRY_PASS_EMPTY_STD(FUNC_GEN6) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC_GEN6) == 0 END
BEGIN TRY_PASS_EMPTY_STD(FUNC_GEN7) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(FUNC_GEN7) == 0 END
BEGIN TRY_PASS_EMPTY_STD(NAME) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(NAME) == 0 END
BEGIN TRY_PASS_EMPTY_STD(ATUPLE) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(ATUPLE) == 0 END
BEGIN TRY_PASS_EMPTY_STD(ATUPLE_PLUS) == 0 END
BEGIN TRY_PASS_EMPTY_BOOST(ATUPLE_PLUS) == 0 END

# else

BEGIN 1 == 1 END

# endif
