// Copyright John Maddock 2008.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define SC_(x) static_cast<typename table_type<T>::type>(BOOST_JOIN(x, L))
   static const std::array<std::array<typename table_type<T>::type, 2>, 260> atanh_data = {{
      {{ SC_(-0.9999983310699462890625), SC_(-6.998237084679026894944012639589359039154) }}, 
      {{ SC_(-0.9999978542327880859375), SC_(-6.872579751329170618373487147465365112113) }}, 
      {{ SC_(-0.999992847442626953125), SC_(-6.270592097465752658938563627507298840894) }}, 
      {{ SC_(-0.999986171722412109375), SC_(-5.94096761408481303343713262973790932529) }}, 
      {{ SC_(-0.9999828338623046875), SC_(-5.832855225378502201556044619092381320613) }}, 
      {{ SC_(-0.99993991851806640625), SC_(-5.206463012087559958576821359878019831698) }}, 
      {{ SC_(-0.9998834133148193359375), SC_(-4.874982184181078415951979525893507005727) }}, 
      {{ SC_(-0.9998509883880615234375), SC_(-4.752279497280337973604777810710354182427) }}, 
      {{ SC_(-0.9996016025543212890625), SC_(-4.260504202858904159948511594799808172686) }}, 
      {{ SC_(-0.9993612766265869140625), SC_(-4.024433435320311666223352267171267188441) }}, 
      {{ SC_(-0.9989283084869384765625), SC_(-3.765564108299923386710499391156429825334) }}, 
      {{ SC_(-0.996978282928466796875), SC_(-3.246782610980921221976121268600437559073) }}, 
      {{ SC_(-0.9950058460235595703125), SC_(-2.995067117994093910327278839368162513171) }}, 
      {{ SC_(-0.9942638874053955078125), SC_(-2.925624274960953689217686505461794222524) }}, 
      {{ SC_(-0.9907157421112060546875), SC_(-2.683964628330836423435761449812536808992) }}, 
      {{ SC_(-0.990334033966064453125), SC_(-2.663723350226517894297444398633783948406) }}, 
      {{ SC_(-0.9760982990264892578125), SC_(-2.207464998348322224447911930881940831481) }}, 
      {{ SC_(-0.975830078125), SC_(-2.201817459680555978821333155035293326746) }}, 
      {{ SC_(-0.972824573516845703125), SC_(-2.142454230829143642120980722087236052917) }}, 
      {{ SC_(-0.964355945587158203125), SC_(-2.004668675602091699628236313728557407114) }}, 
      {{ SC_(-0.9377224445343017578125), SC_(-1.718833734617706500332405835370646660471) }}, 
      {{ SC_(-0.936240673065185546875), SC_(-1.706694048256515431631815885672233674251) }}, 
      {{ SC_(-0.9310147762298583984375), SC_(-1.665954300553314591750124550883856364694) }}, 
      {{ SC_(-0.9284839630126953125), SC_(-1.647283871876074775612635961977408037529) }}, 
      {{ SC_(-0.92702484130859375), SC_(-1.636806734088156154435563850211913111601) }}, 
      {{ SC_(-0.907566547393798828125), SC_(-1.513547347731111539787447924356304439238) }}, 
      {{ SC_(-0.8974773883819580078125), SC_(-1.45909860863314960721683880733144126731) }}, 
      {{ SC_(-0.89201068878173828125), SC_(-1.431681573516303182182852281763191690617) }}, 
      {{ SC_(-0.87765598297119140625), SC_(-1.365471286049011032989426113780315638737) }}, 
      {{ SC_(-0.864722728729248046875), SC_(-1.31177055834445394659795622633703969998) }}, 
      {{ SC_(-0.8482067584991455078125), SC_(-1.249725893334944048400811976134680437019) }}, 
      {{ SC_(-0.805655956268310546875), SC_(-1.114524602859225810626832366034868507092) }}, 
      {{ SC_(-0.8048388957977294921875), SC_(-1.112200609756455010839566277837214844325) }}, 
      {{ SC_(-0.780198574066162109375), SC_(-1.045877833082255648218230662475846475501) }}, 
      {{ SC_(-0.774993419647216796875), SC_(-1.032711173436253036051634900441573817364) }}, 
      {{ SC_(-0.761928558349609375), SC_(-1.000796728136218508274740551519944314059) }}, 
      {{ SC_(-0.7504425048828125), SC_(-0.9739672824457071545212558891424185128762) }}, 
      {{ SC_(-0.7495596408843994140625), SC_(-0.9719492983286864197920761015664283271976) }}, 
      {{ SC_(-0.7481319904327392578125), SC_(-0.9686989420144869980652193388237565424316) }}, 
      {{ SC_(-0.7459518909454345703125), SC_(-0.9637657636705832060477647484529718099027) }}, 
      {{ SC_(-0.740113735198974609375), SC_(-0.9507308314464193446531653906278264445993) }}, 
      {{ SC_(-0.7289731502532958984375), SC_(-0.926532531986765318247140551363667818701) }}, 
      {{ SC_(-0.7226788997650146484375), SC_(-0.9132299082876395943553155758462776294309) }}, 
      {{ SC_(-0.7161557674407958984375), SC_(-0.8997082193533088241452821117641146130759) }}, 
      {{ SC_(-0.7017018795013427734375), SC_(-0.8706453720344795324137713806233471643566) }}, 
      {{ SC_(-0.7013418674468994140625), SC_(-0.869936501309450037837255615527206388543) }}, 
      {{ SC_(-0.6910541057586669921875), SC_(-0.8499705913361888274674607430873352657984) }}, 
      {{ SC_(-0.6847054958343505859375), SC_(-0.8379194558420050243680161679293990119402) }}, 
      {{ SC_(-0.683816432952880859375), SC_(-0.8362476144993315210191634049830735329595) }}, 
      {{ SC_(-0.6747090816497802734375), SC_(-0.8193374156276963181623971413078624527231) }}, 
      {{ SC_(-0.6575610637664794921875), SC_(-0.7885046044142132636244565971148537710756) }}, 
      {{ SC_(-0.6522045135498046875), SC_(-0.7791255597799838714478027675042210151575) }}, 
      {{ SC_(-0.6261923313140869140625), SC_(-0.7351275788820003190814567583601297864281) }}, 
      {{ SC_(-0.62317371368408203125), SC_(-0.7301771459970386661436038123301221899899) }}, 
      {{ SC_(-0.6067488193511962890625), SC_(-0.703759752613062694606340098389939027562) }}, 
      {{ SC_(-0.583805561065673828125), SC_(-0.6682166303197607887197770868026216259119) }}, 
      {{ SC_(-0.57952404022216796875), SC_(-0.6617457665293066615721543987931600748293) }}, 
      {{ SC_(-0.5760939121246337890625), SC_(-0.6565964588573980673410012745962406194225) }}, 
      {{ SC_(-0.56546783447265625), SC_(-0.6408350116350283247012343252539212893332) }}, 
      {{ SC_(-0.557876110076904296875), SC_(-0.6297442839791668262390512924219530985196) }}, 
      {{ SC_(-0.552320957183837890625), SC_(-0.6217149641475686188387285842099747744544) }}, 
      {{ SC_(-0.5396339893341064453125), SC_(-0.6036390747171697770925487100621985766607) }}, 
      {{ SC_(-0.512898921966552734375), SC_(-0.566655625606477066963814924576009526094) }}, 
      {{ SC_(-0.5087778568267822265625), SC_(-0.5610793900942041010450981003401240686441) }}, 
      {{ SC_(-0.49778258800506591796875), SC_(-0.5463539505715040435095828997881628670016) }}, 
      {{ SC_(-0.49138653278350830078125), SC_(-0.5378865967606702825309467710240853062652) }}, 
      {{ SC_(-0.48976075649261474609375), SC_(-0.5357455496477737185521806219611537015446) }}, 
      {{ SC_(-0.4849350452423095703125), SC_(-0.5294166456244710849865087993139911130393) }}, 
      {{ SC_(-0.447905063629150390625), SC_(-0.4820764946679978939347292302681479739796) }}, 
      {{ SC_(-0.4461095333099365234375), SC_(-0.4798325976916710897104638200930521631286) }}, 
      {{ SC_(-0.442959308624267578125), SC_(-0.4759065337156127733780647024545262009008) }}, 
      {{ SC_(-0.4282791614532470703125), SC_(-0.4577873936293679104616687959332340265949) }}, 
      {{ SC_(-0.40590059757232666015625), SC_(-0.4306933608076878743806295938648362982359) }}, 
      {{ SC_(-0.40029656887054443359375), SC_(-0.4240020382545707218942432038948084774496) }}, 
      {{ SC_(-0.3961341381072998046875), SC_(-0.4190551379319939121143095387765817810309) }}, 
      {{ SC_(-0.38362753391265869140625), SC_(-0.4043062717590873607565891081374512836717) }}, 
      {{ SC_(-0.3668625354766845703125), SC_(-0.3847928551425506989844364335805100378916) }}, 
      {{ SC_(-0.36576449871063232421875), SC_(-0.383524642274593405044354917396013672621) }}, 
      {{ SC_(-0.33507001399993896484375), SC_(-0.3485286317501441759395903809820783149922) }}, 
      {{ SC_(-0.325722217559814453125), SC_(-0.3380352468276521988337991305273476357259) }}, 
      {{ SC_(-0.3191967010498046875), SC_(-0.3307524237890151189189282578009246778706) }}, 
      {{ SC_(-0.300002574920654296875), SC_(-0.3095224337886502772553678318463405679109) }}, 
      {{ SC_(-0.296651363372802734375), SC_(-0.3058438250228024929192737081479944125667) }}, 
      {{ SC_(-0.2944457530975341796875), SC_(-0.3034271164344304639271762476267134487006) }}, 
      {{ SC_(-0.287281036376953125), SC_(-0.2956001834724682534366346933530868053292) }}, 
      {{ SC_(-0.277384281158447265625), SC_(-0.2848460820316943671684268151981802822583) }}, 
      {{ SC_(-0.239084422588348388671875), SC_(-0.2438028008332661157406055005411967471601) }}, 
      {{ SC_(-0.23685944080352783203125), SC_(-0.2414442516939151847702918137086915005926) }}, 
      {{ SC_(-0.2253856658935546875), SC_(-0.2293228153248167922436892209184552796457) }}, 
      {{ SC_(-0.222838103771209716796875), SC_(-0.2266405306474514216815783549721039495306) }}, 
      {{ SC_(-0.2155244350433349609375), SC_(-0.2189577360114399508892266906786141753798) }}, 
      {{ SC_(-0.215337574481964111328125), SC_(-0.2187617810795299603977983273974692894269) }}, 
      {{ SC_(-0.210162580013275146484375), SC_(-0.213341433207717363756996669992890635922) }}, 
      {{ SC_(-0.202502727508544921875), SC_(-0.2053409277979887156428422481476167713845) }}, 
      {{ SC_(-0.1915638446807861328125), SC_(-0.1939600847413307517360830799803393191118) }}, 
      {{ SC_(-0.182519435882568359375), SC_(-0.1845877143932293893050297021219967479454) }}, 
      {{ SC_(-0.1746494770050048828125), SC_(-0.1764584460861806768647912839605350112421) }}, 
      {{ SC_(-0.15646183490753173828125), SC_(-0.1577576667718915543508871956214699612414) }}, 
      {{ SC_(-0.155809104442596435546875), SC_(-0.1570886262196417664166148374546893502577) }}, 
      {{ SC_(-0.15365445613861083984375), SC_(-0.1548811251554959173309605263848373765859) }}, 
      {{ SC_(-0.1224990189075469970703125), SC_(-0.1231173360990485154969499505998906050858) }}, 
      {{ SC_(-0.1088167130947113037109375), SC_(-0.1092492929673783678150776545155488162838) }}, 
      {{ SC_(-0.0879255831241607666015625), SC_(-0.08815322150790301444935652020860591175639) }}, 
      {{ SC_(-0.084013283252716064453125), SC_(-0.08421178632314607700273933837777975783955) }}, 
      {{ SC_(-0.06121261417865753173828125), SC_(-0.06128924075509795947903082912344529497771) }}, 
      {{ SC_(-0.0534169971942901611328125), SC_(-0.05346789060550385741248960046699294299987) }}, 
      {{ SC_(-0.0504775941371917724609375), SC_(-0.0505205318923802872873054309617566051014) }}, 
      {{ SC_(-0.029245793819427490234375), SC_(-0.02925413623733265658295603323356445573956) }}, 
      {{ SC_(-0.024859689176082611083984375), SC_(-0.02486481220617492008704879514284669118942) }}, 
      {{ SC_(-0.02046917378902435302734375), SC_(-0.02047203328100153078698678836073163016822) }}, 
      {{ SC_(-0.0188200175762176513671875), SC_(-0.01882224002175634892911631496693598606876) }}, 
      {{ SC_(-0.01615250110626220703125), SC_(-0.01615390607310920464304141890191065645265) }}, 
      {{ SC_(-0.003271550871431827545166015625), SC_(-0.003271562543358962122360602166863449874156) }}, 
      {{ SC_(0.165048140085555239409131900174543261528e-11), SC_(0.1650481400855552394091320500431428863902e-11) }}, 
      {{ SC_(0.2065420751096169738048047292977571487427e-11), SC_(0.206542075109616973804805022998032275867e-11) }}, 
      {{ SC_(0.6933230031758164102484442992135882377625e-11), SC_(0.6933230031758164102484554084849147212117e-11) }}, 
      {{ SC_(0.1335144494962747785393730737268924713135e-10), SC_(0.1335144494962747785393810072036465631048e-10) }}, 
      {{ SC_(0.1639981206391638579589198343455791473389e-10), SC_(0.16399812063916385795893453698677871597e-10) }}, 
      {{ SC_(0.5730159402528300915946601890027523040771e-10), SC_(0.573015940252830091595287349730680977303e-10) }}, 
      {{ SC_(0.1113731329382972035091370344161987304688e-9), SC_(0.1113731329382972035095975242587774151074e-9) }}, 
      {{ SC_(0.1421470718909745301061775535345077514648e-9), SC_(0.1421470718909745301071349514979911287802e-9) }}, 
      {{ SC_(0.3800632031314421510614920407533645629883e-9), SC_(0.3800632031314421510797918354702572472565e-9) }}, 
      {{ SC_(0.6091627202664540163823403418064117431641e-9), SC_(0.6091627202664540164576895507880226730245e-9) }}, 
      {{ SC_(0.1022164131114777774200774729251861572266e-8), SC_(0.1022164131114777774556767071774970972449e-8) }}, 
      {{ SC_(0.2881922256392499548383057117462158203125e-8), SC_(0.288192225639249955636163572505228913693e-8) }}, 
      {{ SC_(0.4762776839584148547146469354629516601563e-8), SC_(0.4762776839584148583159481252584483554515e-8) }}, 
      {{ SC_(0.8854133426439148024655878543853759765625e-8), SC_(0.8854133426439148256031145063819131862293e-8) }}, 
      {{ SC_(0.2305032609228874207474291324615478515625e-7), SC_(0.2305032609228874615709037767878933144975e-7) }}, 
      {{ SC_(0.593924909253473742865025997161865234375e-7), SC_(0.5939249092534744412153923027426683753111e-7) }}, 
      {{ SC_(0.116676488914890796877443790435791015625e-6), SC_(0.1166764889148913263321344126152128598244e-6) }}, 
      {{ SC_(0.23799674409019644372165203094482421875e-6), SC_(0.2379967440902009372945601325325051016678e-6) }}, 
      {{ SC_(0.468465941594331525266170501708984375e-6), SC_(0.4684659415943657951641995164188851895468e-6) }}, 
      {{ SC_(0.938269977268646471202373504638671875e-6), SC_(0.9382699772689218066992955176687134214945e-6) }}, 
      {{ SC_(0.11039855962735600769519805908203125e-5), SC_(0.1103985596274008583684717717140945872849e-5) }}, 
      {{ SC_(0.3291776010883040726184844970703125e-5), SC_(0.3291776010894930389950221585822451540922e-5) }}, 
      {{ SC_(0.7517213816754519939422607421875e-5), SC_(0.7517213816896115440686244581232090380807e-5) }}, 
      {{ SC_(0.1511466689407825469970703125e-4), SC_(0.151146668952292524810471268542864819142e-4) }}, 
      {{ SC_(0.2986399340443313121795654296875e-4), SC_(0.2986399341331127938191098992807520418566e-4) }}, 
      {{ SC_(0.338702811859548091888427734375e-4), SC_(0.3387028119890675897146774685716560577438e-4) }}, 
      {{ SC_(0.906601198948919773101806640625e-4), SC_(0.9066012014327826381277876686817592943471e-4) }}, 
      {{ SC_(0.0002194953267462551593780517578125), SC_(0.0002194953302712183992004987083863068286575) }}, 
      {{ SC_(0.000439521507360041141510009765625), SC_(0.0004395215356621756172832388067699815859163) }}, 
      {{ SC_(0.000633315183222293853759765625), SC_(0.0006333152678940465733692213551976187060063) }}, 
      {{ SC_(0.0011151232756674289703369140625), SC_(0.001115123737886341835357675971473032730801) }}, 
      {{ SC_(0.001962467096745967864990234375), SC_(0.001962469616086656343101697360448994447248) }}, 
      {{ SC_(0.00555375404655933380126953125), SC_(0.005553811147953337251608673095795926379418) }}, 
      {{ SC_(0.0073246769607067108154296875), SC_(0.007324807956742500088707032291078962161983) }}, 
      {{ SC_(0.0086911283433437347412109375), SC_(0.008691347183450786003338484767259005750621) }}, 
      {{ SC_(0.01191294193267822265625), SC_(0.01191350553503790551712224267804100577312) }}, 
      {{ SC_(0.0299333631992340087890625), SC_(0.02994230816857020273948279830730524489901) }}, 
      {{ SC_(0.05124260485172271728515625), SC_(0.05128752666822782052568477249433836427776) }}, 
      {{ SC_(0.05473744869232177734375), SC_(0.05479221508125443603584511761951841280462) }}, 
      {{ SC_(0.0615889132022857666015625), SC_(0.0616669638195183068123715085691638305692) }}, 
      {{ SC_(0.0937536060810089111328125), SC_(0.09402975380882211958737150286272764887037) }}, 
      {{ SC_(0.0944215953350067138671875), SC_(0.09470370926367391764663232920446630580688) }}, 
      {{ SC_(0.0944317281246185302734375), SC_(0.09471393321406025987803236345499696649722) }}, 
      {{ SC_(0.099437296390533447265625), SC_(0.09976699249016486016395021218834536220752) }}, 
      {{ SC_(0.1120129525661468505859375), SC_(0.1124849830355889429846170463794901835017) }}, 
      {{ SC_(0.123102605342864990234375), SC_(0.1237301640233916802510730091978889755734) }}, 
      {{ SC_(0.1356296539306640625), SC_(0.136470609508612479912140462032722727509) }}, 
      {{ SC_(0.137633502483367919921875), SC_(0.1385125786609474520950316182245783874744) }}, 
      {{ SC_(0.1474945545196533203125), SC_(0.1485782998046483417449330959482946265992) }}, 
      {{ SC_(0.161897182464599609375), SC_(0.163334331667904492581671816308390067676) }}, 
      {{ SC_(0.170511066913604736328125), SC_(0.1721929869363735526203806645402133592389) }}, 
      {{ SC_(0.170518338680267333984375), SC_(0.1722004764629990854468820758107275101776) }}, 
      {{ SC_(0.1856291294097900390625), SC_(0.1878064731815008603331954720151705472375) }}, 
      {{ SC_(0.188988208770751953125), SC_(0.1912876932893582208057658208318612269009) }}, 
      {{ SC_(0.23206615447998046875), SC_(0.2363721243391452317792593455314737006561) }}, 
      {{ SC_(0.23480379581451416015625), SC_(0.2392675448267426979141683097039145007706) }}, 
      {{ SC_(0.2646920680999755859375), SC_(0.271147290330230077990690324884981633411) }}, 
      {{ SC_(0.27949869632720947265625), SC_(0.287138205934443267371748751232576707068) }}, 
      {{ SC_(0.2878930568695068359375), SC_(0.2962673858386818760743191664695461734079) }}, 
      {{ SC_(0.29259669780731201171875), SC_(0.3014037366523923189237135841752647308831) }}, 
      {{ SC_(0.310164928436279296875), SC_(0.3207278827697849622886565511333038417846) }}, 
      {{ SC_(0.31092464923858642578125), SC_(0.321568689394455844291345719208687599694) }}, 
      {{ SC_(0.31145012378692626953125), SC_(0.3221505056451928603025272541357270447997) }}, 
      {{ SC_(0.3271782398223876953125), SC_(0.339664946169947805312598648269899355103) }}, 
      {{ SC_(0.3574345111846923828125), SC_(0.3739415343654542284372692080523616239196) }}, 
      {{ SC_(0.35936939716339111328125), SC_(0.3761615922309022376646457727379955025932) }}, 
      {{ SC_(0.35960352420806884765625), SC_(0.376430465969337149289307260724149960882) }}, 
      {{ SC_(0.36268270015716552734375), SC_(0.3799714809649667171023188371749734970087) }}, 
      {{ SC_(0.38961827754974365234375), SC_(0.4113499159905352551364249131399712396774) }}, 
      {{ SC_(0.3904266357421875), SC_(0.4123033008021400244918917567469642157756) }}, 
      {{ SC_(0.39811360836029052734375), SC_(0.4214052375603138698610724190400452927171) }}, 
      {{ SC_(0.41150724887847900390625), SC_(0.4374243870957909704688484065341332398408) }}, 
      {{ SC_(0.4120509624481201171875), SC_(0.438079118237434936389702798117975965435) }}, 
      {{ SC_(0.41868770122528076171875), SC_(0.4460997186945703025969092767121680940872) }}, 
      {{ SC_(0.4213654994964599609375), SC_(0.4493511447897728936819156706558562191803) }}, 
      {{ SC_(0.4516327381134033203125), SC_(0.4867494899047367852559605895765155264024) }}, 
      {{ SC_(0.45386397838592529296875), SC_(0.4895560176112374916678981118235718808999) }}, 
      {{ SC_(0.46555078029632568359375), SC_(0.5043748446613432644265368658689527269923) }}, 
      {{ SC_(0.48124635219573974609375), SC_(0.5246050193978662791545681419778026197093) }}, 
      {{ SC_(0.4862163066864013671875), SC_(0.5310932154891663069921221983411478381793) }}, 
      {{ SC_(0.48987305164337158203125), SC_(0.5358932909903700749055224508487508967699) }}, 
      {{ SC_(0.502483844757080078125), SC_(0.5526234425942533333497748944931154989089) }}, 
      {{ SC_(0.5074074268341064453125), SC_(0.5592320547729961905967162368329222398251) }}, 
      {{ SC_(0.50932216644287109375), SC_(0.5618140818296767894194272698942620318047) }}, 
      {{ SC_(0.5143489837646484375), SC_(0.5686253097655145706303693099066490153432) }}, 
      {{ SC_(0.5154285430908203125), SC_(0.5700943191671631767428541315644447456688) }}, 
      {{ SC_(0.5234100818634033203125), SC_(0.5810250825991417575052326759884079418476) }}, 
      {{ SC_(0.527447223663330078125), SC_(0.5866018515043636631396476583453040565353) }}, 
      {{ SC_(0.5309803485870361328125), SC_(0.5915094458340507867953465738492206518901) }}, 
      {{ SC_(0.5477793216705322265625), SC_(0.615203099922968789392563855123133294432) }}, 
      {{ SC_(0.5577394962310791015625), SC_(0.6295459624918965701898701179409871117547) }}, 
      {{ SC_(0.558278560638427734375), SC_(0.6303287742357745455208912482493484545533) }}, 
      {{ SC_(0.5843560695648193359375), SC_(0.6690521906099504661825293948827417437714) }}, 
      {{ SC_(0.58713626861572265625), SC_(0.6732844960442398402790025978545962104976) }}, 
      {{ SC_(0.587891101837158203125), SC_(0.6744372167164567880432735474026414000937) }}, 
      {{ SC_(0.59034061431884765625), SC_(0.6781887236623533927703862563260315607999) }}, 
      {{ SC_(0.5945003032684326171875), SC_(0.684597775489552027208298442397307157518) }}, 
      {{ SC_(0.59579753875732421875), SC_(0.686606510213166504449707209584681039618) }}, 
      {{ SC_(0.5961520671844482421875), SC_(0.6871563252400655469505525776500166368732) }}, 
      {{ SC_(0.6005008220672607421875), SC_(0.6939300827887145041443406444376899260605) }}, 
      {{ SC_(0.6150004863739013671875), SC_(0.7169242329194352408292786592670047478091) }}, 
      {{ SC_(0.6162893772125244140625), SC_(0.718999805549710753361347867549953319929) }}, 
      {{ SC_(0.6194069385528564453125), SC_(0.7240422748778544828852658840299672671419) }}, 
      {{ SC_(0.62850666046142578125), SC_(0.7389438896054792083692956718844417119266) }}, 
      {{ SC_(0.6293842792510986328125), SC_(0.7403958734869583175052115753160537014235) }}, 
      {{ SC_(0.641617298126220703125), SC_(0.7609178886018203867768570854585116358618) }}, 
      {{ SC_(0.6424276828765869140625), SC_(0.7622965466812235050838313606820415418465) }}, 
      {{ SC_(0.643742084503173828125), SC_(0.7645378650643100878873860589743586646119) }}, 
      {{ SC_(0.6468508243560791015625), SC_(0.7698647951781609863982541324817513695536) }}, 
      {{ SC_(0.661591053009033203125), SC_(0.7956379107512945487304413199186620675388) }}, 
      {{ SC_(0.669950008392333984375), SC_(0.8106524185805045490744310224054265386895) }}, 
      {{ SC_(0.6813662052154541015625), SC_(0.8316597473423232100463095038096703191227) }}, 
      {{ SC_(0.6968657970428466796875), SC_(0.861181279065929560893743209045814929335) }}, 
      {{ SC_(0.69818878173828125), SC_(0.8637579113749143394578018115274414809977) }}, 
      {{ SC_(0.7447831630706787109375), SC_(0.9611360201710216528199065175645895870797) }}, 
      {{ SC_(0.7518312931060791015625), SC_(0.9771540941752986840434538256366757954859) }}, 
      {{ SC_(0.753439426422119140625), SC_(0.9808634133542228689124290623167881332843) }}, 
      {{ SC_(0.7567856311798095703125), SC_(0.9886489208209698781284720459176833206972) }}, 
      {{ SC_(0.781728267669677734375), SC_(1.049799171982895646434536064092181561482) }}, 
      {{ SC_(0.8115026950836181640625), SC_(1.13141414441875853996497925019619325203) }}, 
      {{ SC_(0.8146479129791259765625), SC_(1.140694775558441751141132778566016009721) }}, 
      {{ SC_(0.8266689777374267578125), SC_(1.177523083369968036920622084233473757412) }}, 
      {{ SC_(0.8313877582550048828125), SC_(1.192613822570143323026054454909849845367) }}, 
      {{ SC_(0.83430385589599609375), SC_(1.202132332303961254992584113869002097916) }}, 
      {{ SC_(0.8416652679443359375), SC_(1.226857064433516189319978030947873360052) }}, 
      {{ SC_(0.85844135284423828125), SC_(1.287389667157365295856183109440507044615) }}, 
      {{ SC_(0.8678996562957763671875), SC_(1.324504043392939851303939943290664403741) }}, 
      {{ SC_(0.8679344654083251953125), SC_(1.324645130926160707039340948233535357264) }}, 
      {{ SC_(0.8800599575042724609375), SC_(1.376033487778217701076449442272672558944) }}, 
      {{ SC_(0.900353908538818359375), SC_(1.474085296119410724208965851686464033313) }}, 
      {{ SC_(0.909944057464599609375), SC_(1.527199085186199482742984118973093246969) }}, 
      {{ SC_(0.9142425060272216796875), SC_(1.552776894827300414876536030721521100862) }}, 
      {{ SC_(0.9149219989776611328125), SC_(1.556931837197935820864845994901316963069) }}, 
      {{ SC_(0.918490886688232421875), SC_(1.579289662838161287932020406856633071199) }}, 
      {{ SC_(0.91889286041259765625), SC_(1.581866335942762753838726342973068309432) }}, 
      {{ SC_(0.919395923614501953125), SC_(1.585108284332000716017087123313008639504) }}, 
      {{ SC_(0.9296839237213134765625), SC_(1.656055522329536863771852381599659592259) }}, 
      {{ SC_(0.929839611053466796875), SC_(1.6572041418041492174589693922205395004) }}, 
      {{ SC_(0.9352962970733642578125), SC_(1.69909864336192661532416961180476119912) }}, 
      {{ SC_(0.937641620635986328125), SC_(1.718164398807964846646482962730886595432) }}, 
      {{ SC_(0.9410912990570068359375), SC_(1.747508407724663243809927907107491941803) }}, 
      {{ SC_(0.96212291717529296875), SC_(1.973718016345510160586126866935192466646) }}, 
      {{ SC_(0.974821567535400390625), SC_(2.181122777108378116176876272951870175355) }}, 
      {{ SC_(0.976945400238037109375), SC_(2.225721449969825368029395993270411091677) }}, 
      {{ SC_(0.985663890838623046875), SC_(2.465463560165053850170474882317767473537) }}, 
      {{ SC_(0.98803806304931640625), SC_(2.556586922814200202245299755886851508354) }}, 
      {{ SC_(0.9928233623504638671875), SC_(2.813238353909419376753157236188061588174) }}
   }};
#undef SC_

