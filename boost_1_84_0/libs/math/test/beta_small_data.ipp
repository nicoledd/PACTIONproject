//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using std::ldexp;

   static const std::array<std::array<typename table_type<T>::type, 3>, 23> beta_small_data = { {
      {{ SC_(0.1730655412757187150418758392333984375e-5), SC_(0.1730655412757187150418758392333984375e-5), SC_(1155631.551635027016649268884796909927277) }},
      {{ SC_(0.216575062950141727924346923828125e-5), SC_(0.1730655412757187150418758392333984375e-5), SC_(1039549.452063747329381617654200841254652) }},
      {{ SC_(0.216575062950141727924346923828125e-5), SC_(0.216575062950141727924346923828125e-5), SC_(923467.3524924676425690820378921903570447) }},
      {{ SC_(0.72700195232755504548549652099609375e-5), SC_(0.1730655412757187150418758392333984375e-5), SC_(715366.9882608199489156088500706918884474) }},
      {{ SC_(0.72700195232755504548549652099609375e-5), SC_(0.216575062950141727924346923828125e-5), SC_(599284.8886895402674421924477675861806454) }},
      {{ SC_(0.72700195232755504548549652099609375e-5), SC_(0.72700195232755504548549652099609375e-5), SC_(275102.4248866129549503632384850636365051) }},
      {{ SC_(0.14000004739500582218170166015625e-4), SC_(0.1730655412757187150418758392333984375e-5), SC_(649244.3230419389091055494874477610228104) }},
      {{ SC_(0.14000004739500582218170166015625e-4), SC_(0.216575062950141727924346923828125e-5), SC_(533162.2234706592346717218748633348329255) }},
      {{ SC_(0.14000004739500582218170166015625e-4), SC_(0.72700195232755504548549652099609375e-5), SC_(208979.7596677320047637517666971035319141) }},
      {{ SC_(0.14000004739500582218170166015625e-4), SC_(0.14000004739500582218170166015625e-4), SC_(142857.0944488511634633414470885024340382) }},
      {{ SC_(0.17196454791701398789882659912109375e-4), SC_(0.1730655412757187150418758392333984375e-5), SC_(635967.2966761120408405544283252854891051) }},
      {{ SC_(0.17196454791701398789882659912109375e-4), SC_(0.216575062950141727924346923828125e-5), SC_(519885.1971048323697502064131013612391462) }},
      {{ SC_(0.17196454791701398789882659912109375e-4), SC_(0.72700195232755504548549652099609375e-5), SC_(195702.7333019051790657557774869592595619) }},
      {{ SC_(0.17196454791701398789882659912109375e-4), SC_(0.14000004739500582218170166015625e-4), SC_(129580.0680830243894812627999923649303964) }},
      {{ SC_(0.17196454791701398789882659912109375e-4), SC_(0.17196454791701398789882659912109375e-4), SC_(116303.0417171976400618567245671670493297) }},
      {{ SC_(0.60085076256655156612396240234375e-4), SC_(0.1730655412757187150418758392333984375e-5), SC_(594458.8435535046961781833034439822035691) }},
      {{ SC_(0.60085076256655156612396240234375e-4), SC_(0.216575062950141727924346923828125e-5), SC_(478376.7439822250699480739314338811953296) }},
      {{ SC_(0.60085076256655156612396240234375e-4), SC_(0.72700195232755504548549652099609375e-5), SC_(154194.2801792984055346510144704788292969) }},
      {{ SC_(0.60085076256655156612396240234375e-4), SC_(0.14000004739500582218170166015625e-4), SC_(88071.61496041830983457563965247905141874) }},
      {{ SC_(0.60085076256655156612396240234375e-4), SC_(0.17196454791701398789882659912109375e-4), SC_(74794.58859459188997822531683455967562115) }},
      {{ SC_(0.60085076256655156612396240234375e-4), SC_(0.60085076256655156612396240234375e-4), SC_(33286.13547199056171829186221434157408122) }},
      {{ std::numeric_limits<T>::max_exponent > 154 ? SC_(2.9833362924800826973163861261851735349505886138402e-154) : SC_(1.0), std::numeric_limits<T>::max_exponent > 154 ? SC_(3.7291703656001033716454826577314669186882357673002e-155) : SC_(1.0), std::numeric_limits<T>::max_exponent > 154 ? SC_(3.0167567842370843474041556245963153786828573096333e154) : SC_(1.0) }},
      {{ std::numeric_limits<T>::max_exponent > 154 ? SC_(2.9833362924800826973163861261851735349505886138402e-154) : SC_(1.0), std::numeric_limits<T>::max_exponent > 154 ? SC_(1.8645851828000516858227413288657334593441178836501e-155) : SC_(1.0), std::numeric_limits<T>::max_exponent > 154 ? SC_(5.6983183702256037673189606242374846041787304737518e154) : SC_(1.0) }}
   } };


// static_cast<T>(BOOST_JOIN(x, L)) // changes type to T and adds suffix L to decimal digit string.


