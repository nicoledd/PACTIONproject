// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007, 2009
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning(disable : 4756) // overflow in constant arithmetic
// Constants are too big for float case, but this doesn't matter for test.
#endif

#include <boost/math/concepts/real_concept.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/array.hpp>
#include "functor.hpp"

#include "handle_test_result.hpp"
#include "table_type.hpp"

#include <boost/math/special_functions/jacobi_elliptic.hpp>

#ifndef SC_
#define SC_(x) static_cast<typename table_type<T>::type>(BOOST_JOIN(x, L))
#endif

template <class Real, typename T>
void do_test_sn(T& data, const char* type_name, const char* test)
{
#if !(defined(ERROR_REPORTING_MODE) && !defined(SN_FUNCTION_TO_TEST))
   typedef Real                   value_type;

   std::cout << "Testing: " << test << std::endl;

#ifdef SN_FUNCTION_TO_TEST
   value_type(*fp2)(value_type, value_type) = SN_FUNCTION_TO_TEST;
#elif defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
    value_type (*fp2)(value_type, value_type) = boost::math::jacobi_sn<value_type, value_type>;
#else
    value_type (*fp2)(value_type, value_type) = boost::math::jacobi_sn;
#endif
    boost::math::tools::test_result<value_type> result;

    result = boost::math::tools::test_hetero<Real>(
      data,
      bind_func<Real>(fp2, 1, 0),
      extract_result<Real>(2));
    handle_test_result(result, data[result.worst()], result.worst(),
      type_name, "jacobi_sn", test);

#ifdef CN_FUNCTION_TO_TEST
    fp2 = CN_FUNCTION_TO_TEST;
#elif defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
    fp2 = boost::math::jacobi_cn<value_type, value_type>;
#else
    fp2 = boost::math::jacobi_cn;
#endif
    result = boost::math::tools::test_hetero<Real>(
      data,
      bind_func<Real>(fp2, 1, 0),
      extract_result<Real>(3));
    handle_test_result(result, data[result.worst()], result.worst(),
      type_name, "jacobi_cn", test);

#ifdef SN_FUNCTION_TO_TEST
    fp2 = DN_FUNCTION_TO_TEST;
#elif defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
    fp2 = boost::math::jacobi_dn<value_type, value_type>;
#else
    fp2 = boost::math::jacobi_dn;
#endif
    result = boost::math::tools::test_hetero<Real>(
      data,
      bind_func<Real>(fp2, 1, 0),
      extract_result<Real>(4));
    handle_test_result(result, data[result.worst()], result.worst(),
      type_name, "jacobi_dn", test);

   std::cout << std::endl;
#endif
}

template <typename T>
void test_spots(T, const char* type_name)
{
    BOOST_MATH_STD_USING
    // Function values calculated on http://functions.wolfram.com/
    // Note that Mathematica's Sn/Cn/Dn accepts k^2 as the second parameter.
    // Arguments here are theta, k, sn, cn, dn
    static const std::array<std::array<T, 5>, 36> data1 = {{
        {{ SC_(0.0), SC_(0.0), SC_(0.0), SC_(1.0), SC_(1.0) }},
        {{ ldexp(T(1), -25), ldexp(T(1), -25), SC_(2.98023223876953080883700663838486782870427050521881839342311e-8), SC_(0.99999999999999955591079014993741669975171697261290223678373), SC_(0.99999999999999999999999999999960556954738949421406900774443) }},
        {{ -ldexp(T(1), -25), ldexp(T(1), -25), SC_(-2.98023223876953080883700663838486782870427050521881839342311e-8), SC_(0.99999999999999955591079014993741669975171697261290223678373), SC_(0.99999999999999999999999999999960556954738949421406900774443) }},
        {{ SC_(0.25), ldexp(T(1), -25), SC_(0.247403959254522927383635623557663763268693729825996390997241), SC_(0.968912421710644784709721529742747886950140086772629513814665), SC_(0.99999999999999997281786831901333837240938011109848356555885) }},
        {{ SC_(-0.25), ldexp(T(1), -25), SC_(-0.247403959254522927383635623557663763268693729825996390997241), SC_(0.968912421710644784709721529742747886950140086772629513814665), SC_(0.99999999999999997281786831901333837240938011109848356555885) }},
        {{ SC_(1.25), ldexp(T(1), -25), SC_(0.948984619355586147780156037971989352776684194861616269831136), SC_(0.315322362395268865789580233344649598639316847638615703458263), SC_(0.99999999999999960006577747263860127231780811081154547949983) }},
        {{ SC_(-1.25), ldexp(T(1), -25), SC_(-0.948984619355586147780156037971989352776684194861616269831136), SC_(0.315322362395268865789580233344649598639316847638615703458263), SC_(0.99999999999999960006577747263860127231780811081154547949983) }},
        {{ SC_(25.0), ldexp(T(1), -25), SC_(-0.132351750097778560056127137329035522219365438979106560464704), SC_(0.991202811863472859528158119981178957382802975691690722810123), SC_(0.99999999999999999222089563757583834413059580275315226870704) }},
        {{ SC_(-25.0), ldexp(T(1), -25), SC_(0.132351750097778560056127137329035522219365438979106560464704), SC_(0.991202811863472859528158119981178957382802975691690722810123), SC_(0.99999999999999999222089563757583834413059580275315226870704) }},
        {{ ldexp(T(1), -25), SC_(0.5), SC_(2.9802322387695306985462582979816394722771323981194501394601e-8), SC_(0.99999999999999955591079014993744956895610118130967536624417), SC_(0.99999999999999988897769753748438088116649141278818704012037) }},
        {{ -ldexp(T(1), -25), SC_(0.5), SC_(-2.9802322387695306985462582979816394722771323981194501394601e-8), SC_(0.99999999999999955591079014993744956895610118130967536624417), SC_(0.99999999999999988897769753748438088116649141278818704012037) }},
        {{ SC_(0.25), SC_(0.5), SC_(0.246781405136141600483623741101255389743847413013817188632739), SC_(0.969071172865559727608777289021929824625726812182428398055476), SC_(0.992358168465276394946615469032829292963938826683866720698130) }},
        {{ SC_(-0.25), SC_(0.5), SC_(-0.246781405136141600483623741101255389743847413013817188632739), SC_(0.969071172865559727608777289021929824625726812182428398055476), SC_(0.992358168465276394946615469032829292963938826683866720698130) }},
        {{ SC_(1.25), SC_(0.5), SC_(0.928561236426319775700204388269999130782711902203415239399579), SC_(0.371179242693370810357222594552131893184749696381729988511999), SC_(0.885688154799196841458565445994481097477880319663264816077719) }},
        {{ SC_(-1.25), SC_(0.5), SC_(-0.928561236426319775700204388269999130782711902203415239399579), SC_(0.371179242693370810357222594552131893184749696381729988511999), SC_(0.885688154799196841458565445994481097477880319663264816077719) }},
        {{ SC_(25.0), SC_(0.5), SC_(-0.969223071486651608400225080456020493867827336842041561445359), SC_(-0.246184154035106463351874891855925292474628176040625311168501), SC_(0.874729477852721764836147376110255133761608728373832418508248) }},
        {{ SC_(-25.0), SC_(0.5), SC_(0.969223071486651608400225080456020493867827336842041561445359), SC_(-0.246184154035106463351874891855925292474628176040625311168501), SC_(0.874729477852721764836147376110255133761608728373832418508248) }},
        {{ ldexp(T(1), -25), 1 - ldexp(T(1), -9), SC_(2.98023223876953036939562331632512854347569015560128614888589e-8), SC_(0.99999999999999955591079014993754766348947956082687878223721), SC_(0.99999999999999955764381956001984590118394542979655101564079) }},
        {{ -ldexp(T(1), -25), 1 - ldexp(T(1), -9), SC_(-2.98023223876953036939562331632512854347569015560128614888589e-8), SC_(0.99999999999999955591079014993754766348947956082687878223721), SC_(0.99999999999999955764381956001984590118394542979655101564079) }},
        {{ SC_(0.25), 1 - ldexp(T(1), -9), SC_(0.244928335616519632082236089277654937383208524525331032303082), SC_(0.969541185516180906431546524888118346090913555188425579774305), SC_(0.969661908643964623248878987955178702010392829596222190545649) }},
        {{ SC_(-0.25), 1 - ldexp(T(1), -9), SC_(-0.244928335616519632082236089277654937383208524525331032303082), SC_(0.969541185516180906431546524888118346090913555188425579774305), SC_(0.969661908643964623248878987955178702010392829596222190545649) }},
        {{ SC_(1.25), 1 - ldexp(T(1), -9), SC_(0.848768940045053312079390719205939167551169094157365783446523), SC_(0.528763923140371497228677918580246099580380684604621321430057), SC_(0.531415689278260818860813380561526095359692710060403584603095) }},
        {{ SC_(-1.25), 1 - ldexp(T(1), -9), SC_(-0.848768940045053312079390719205939167551169094157365783446523), SC_(0.528763923140371497228677918580246099580380684604621321430057), SC_(0.531415689278260818860813380561526095359692710060403584603095) }},
        {{ SC_(25.0), 1 - ldexp(T(1), -9), SC_(-0.0252326124525503880903568715488227138184083895871544015366337), SC_(-0.999681606947341709011836635135181960590782564534371631099332), SC_(0.999682849652724146508471774051629114156076052044812654903417) }},
        {{ SC_(-25.0), 1 - ldexp(T(1), -9), SC_(0.0252326124525503880903568715488227138184083895871544015366337), SC_(-0.999681606947341709011836635135181960590782564534371631099332), SC_(0.999682849652724146508471774051629114156076052044812654903417) }},

        // Try modulus > 1
        {{ ldexp(T(1), -25), SC_(1.5), SC_(2.98023223876952981622027157475276613133414644789222481971590e-8), SC_(0.999999999999999555910790149937712522591174851747994454928040), SC_(0.999999999999999000799277837359575841918151654603571877092161) }},
        {{ -ldexp(T(1), -25), SC_(1.5), SC_(-2.98023223876952981622027157475276613133414644789222481971590e-8), SC_(0.999999999999999555910790149937712522591174851747994454928040), SC_(0.999999999999999000799277837359575841918151654603571877092161) }},
        {{ SC_(0.25), SC_(1.5), SC_(0.241830488135945315134822478837394038661484435596992059686086), SC_(0.970318512143270619246031961334217540099946232418710982266812), SC_(0.931888155181641649031244632258710371461078255228024421800363) }},
        {{ SC_(-0.25), SC_(1.5), SC_(-0.241830488135945315134822478837394038661484435596992059686086), SC_(0.970318512143270619246031961334217540099946232418710982266812), SC_(0.931888155181641649031244632258710371461078255228024421800363) }},
        {{ SC_(1.25), SC_(1.5), SC_(0.665875890711922169121186264316618499018039094009893317545462), SC_(0.746062529663971452521312655373498959968622875614588791642250), SC_(-0.0486921028438866868299166778939466685768843580182675008164949) }},
        {{ SC_(-1.25), SC_(1.5), SC_(-0.665875890711922169121186264316618499018039094009893317545462), SC_(0.746062529663971452521312655373498959968622875614588791642250), SC_(-0.0486921028438866868299166778939466685768843580182675008164949) }},
        {{ SC_(25.0), SC_(1.5), SC_(0.618665338981368217712277210270169521641154921220796362724248), SC_(0.785654630447163313102421517325310755764805805534154371583941), SC_(0.372585153048138377269609818284480926623056458773704266654150) }},
        {{ SC_(-25.0), SC_(1.5), SC_(-0.618665338981368217712277210270169521641154921220796362724248), SC_(0.785654630447163313102421517325310755764805805534154371583941), SC_(0.372585153048138377269609818284480926623056458773704266654150) }},

        // Special Values:
        {{ SC_(0.0), SC_(0.5), SC_(0.0), SC_(1.0), SC_(1.0) }},
        {{ SC_(5.0), SC_(0.0), SC_(-0.958924274663138468893154406155993973352461543964601778131672), SC_(0.283662185463226264466639171513557308334422592252215944930359), SC_(1.0) }},
        {{ SC_(5.0), SC_(1.0), SC_(0.999909204262595131210990447534473021089812615990547862736429), SC_(0.0134752822213045573055191382448821552908373539417006868332819), SC_(0.0134752822213045573055191382448821552908373539417006868332819) }},
    }};
    do_test_sn<T>(data1, type_name, "Jacobi Elliptic: Mathworld Data");

#include "jacobi_elliptic.ipp"
    do_test_sn<T>(jacobi_elliptic, type_name, "Jacobi Elliptic: Random Data");
#include "jacobi_elliptic_small.ipp"
    do_test_sn<T>(jacobi_elliptic_small, type_name, "Jacobi Elliptic: Random Small Values");
#include "jacobi_near_1.ipp"
    do_test_sn<T>(jacobi_near_1, type_name, "Jacobi Elliptic: Modulus near 1");
#include "jacobi_large_phi.ipp"
    do_test_sn<T>(jacobi_large_phi, type_name, "Jacobi Elliptic: Large Phi");

    //
    // Sanity checks for all the various derived functions - these are all
    // trivial wrappers around the main three that are tested above - so just
    // use a simple sanity check for each one.
    // Test values are from functions.wolfram.com:
    //
    T tol = boost::math::tools::epsilon<T>() * 100;
    boost::math::policies::policy<> pol;
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_cd(T(0.5), T(0.5)), static_cast<T>(0.905869360370352996327275878479104183407762212476128499788493L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_cd(T(0.5), T(0.5), pol), static_cast<T>(0.905869360370352996327275878479104183407762212476128499788493L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_cn(T(0.5), T(0.5)), static_cast<T>(0.879941022963758342138211939938800035594045353539382810624647L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_cn(T(0.5), T(0.5), pol), static_cast<T>(0.879941022963758342138211939938800035594045353539382810624647L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_cs(T(0.5), T(0.5)), static_cast<T>(1.85218402142505803268146025319200184620073865036147924150565L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_cs(T(0.5), T(0.5), pol), static_cast<T>(1.85218402142505803268146025319200184620073865036147924150565L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_dc(T(0.5), T(0.5)), static_cast<T>(1.10391193669599654696698383614539220889596741980833071370343L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_dc(T(0.5), T(0.5), pol), static_cast<T>(1.10391193669599654696698383614539220889596741980833071370343L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_dn(T(0.5), T(0.5)), static_cast<T>(0.971377398838178842823315157470233933307542433588855341182382L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_dn(T(0.5), T(0.5), pol), static_cast<T>(0.971377398838178842823315157470233933307542433588855341182382L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_ds(T(0.5), T(0.5)), static_cast<T>(2.04464805020871497502900445828888632133468724223115900866414L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_ds(T(0.5), T(0.5), pol), static_cast<T>(2.04464805020871497502900445828888632133468724223115900866414L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_nc(T(0.5), T(0.5)), static_cast<T>(1.13643979983097851593855424992691981204889859711476187519109L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_nc(T(0.5), T(0.5), pol), static_cast<T>(1.13643979983097851593855424992691981204889859711476187519109L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_nd(T(0.5), T(0.5)), static_cast<T>(1.02946599457230050141998045852435702297405263760707971258676L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_nd(T(0.5), T(0.5), pol), static_cast<T>(1.02946599457230050141998045852435702297405263760707971258676L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_ns(T(0.5), T(0.5)), static_cast<T>(2.10489563855842977359275221390569031239706339764770047142101L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_ns(T(0.5), T(0.5), pol), static_cast<T>(2.10489563855842977359275221390569031239706339764770047142101L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_sc(T(0.5), T(0.5)), static_cast<T>(0.539903156723383602910722041849329275299051877814755451255071L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_sc(T(0.5), T(0.5), pol), static_cast<T>(0.539903156723383602910722041849329275299051877814755451255071L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_sd(T(0.5), T(0.5)), static_cast<T>(0.489081727242945953222289853693492188561192086497066116267160L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_sd(T(0.5), T(0.5), pol), static_cast<T>(0.489081727242945953222289853693492188561192086497066116267160L), tol);

    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_sn(T(0.5), T(0.5)), static_cast<T>(0.475082936028536510082218324703870258745078171807428948028252L), tol);
    BOOST_CHECK_CLOSE_FRACTION(boost::math::jacobi_sn(T(0.5), T(0.5), pol), static_cast<T>(0.475082936028536510082218324703870258745078171807428948028252L), tol);

}

