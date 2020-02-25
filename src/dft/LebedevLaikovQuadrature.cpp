//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "LebedevLaikovQuadrature.hpp"

#include <cmath>

CLebedevLaikovQuadrature::CLebedevLaikovQuadrature(const int32_t nAngularPoints)

    : _nAngularPoints(nAngularPoints)
{
}

CLebedevLaikovQuadrature::~CLebedevLaikovQuadrature()
{
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::generate() const
{
    if (_nAngularPoints == 6) return _generateQuadratureWith6Points();

    if (_nAngularPoints == 50) return _generateQuadratureWith50Points();

    if (_nAngularPoints == 110) return _generateQuadratureWith110Points();

    if (_nAngularPoints == 194) return _generateQuadratureWith194Points();

    if (_nAngularPoints == 302) return _generateQuadratureWith302Points();

    if (_nAngularPoints == 434) return _generateQuadratureWith434Points();

    if (_nAngularPoints == 590) return _generateQuadratureWith590Points();

    if (_nAngularPoints == 770) return _generateQuadratureWith770Points();

    if (_nAngularPoints == 974) return _generateQuadratureWith974Points();

    if (_nAngularPoints == 2030) return _generateQuadratureWith2030Points();

    return CMemBlock2D<double>();
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::_generateQuadratureWith6Points() const
{
    CMemBlock2D<double> qpoints(_nAngularPoints, 4);

    _generateCaseOne(qpoints, 0, 0.1666666666666667);

    return qpoints;
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::_generateQuadratureWith50Points() const
{
    CMemBlock2D<double> qpoints(_nAngularPoints, 4);

    _generateCaseOne(qpoints, 0, 0.1269841269841270e-1);

    _generateCaseTwo(qpoints, 6, 0.2257495590828924e-1);

    _generateCaseThree(qpoints, 18, 0.2109375000000000e-1);

    _generateCaseFour(qpoints, 26, 0.3015113445777636, 0.2017333553791887e-1);

    return qpoints;
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::_generateQuadratureWith110Points() const
{
    CMemBlock2D<double> qpoints(_nAngularPoints, 4);

    _generateCaseOne(qpoints, 0, 0.3828270494937162e-2);

    _generateCaseThree(qpoints, 6, 0.9793737512487512e-2);

    _generateCaseFour(qpoints, 14, 0.1851156353447362, 0.8211737283191111e-2);

    _generateCaseFour(qpoints, 38, 0.6904210483822922, 0.9942814891178103e-2);

    _generateCaseFour(qpoints, 62, 0.3956894730559419, 0.9595471336070963e-2);

    _generateCaseFive(qpoints, 86, 0.4783690288121502, 0.9694996361663028e-2);

    return qpoints;
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::_generateQuadratureWith194Points() const
{
    CMemBlock2D<double> qpoints(_nAngularPoints, 4);

    _generateCaseOne(qpoints, 0, 0.1782340447244611e-2);

    _generateCaseTwo(qpoints, 6, 0.5716905949977102e-2);

    _generateCaseThree(qpoints, 18, 0.5573383178848738e-2);

    _generateCaseFour(qpoints, 26, 0.6712973442695226, 0.5608704082587997e-2);

    _generateCaseFour(qpoints, 50, 0.2892465627575439, 0.5158237711805383e-2);

    _generateCaseFour(qpoints, 74, 0.4446933178717437, 0.5518771467273614e-2);

    _generateCaseFour(qpoints, 98, 0.1299335447650067, 0.4106777028169394e-2);

    _generateCaseFive(qpoints, 122, 0.3457702197611283, 0.5051846064614808e-2);

    _generateCaseSix(qpoints, 146, 0.1590417105383530, 0.8360360154824589, 0.5530248916233094e-2);

    return qpoints;
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::_generateQuadratureWith302Points() const
{
    CMemBlock2D<double> qpoints(_nAngularPoints, 4);

    _generateCaseOne(qpoints, 0, 0.8545911725128148e-3);

    _generateCaseThree(qpoints, 6, 0.3599119285025571e-2);

    _generateCaseFour(qpoints, 14, 0.3515640345570105, 0.3449788424305883e-2);

    _generateCaseFour(qpoints, 38, 0.6566329410219612, 0.3604822601419882e-2);

    _generateCaseFour(qpoints, 62, 0.4729054132581005, 0.3576729661743367e-2);

    _generateCaseFour(qpoints, 86, 0.9618308522614784e-1, 0.2352101413689164e-2);

    _generateCaseFour(qpoints, 110, 0.2219645236294178, 0.3108953122413675e-2);

    _generateCaseFour(qpoints, 134, 0.7011766416089545, 0.3650045807677255e-2);

    _generateCaseFive(qpoints, 158, 0.2644152887060663, 0.2982344963171804e-2);

    _generateCaseFive(qpoints, 182, 0.5718955891878961, 0.3600820932216460e-2);

    _generateCaseSix(qpoints, 206, 0.2510034751770465, 0.8000727494073952, 0.3571540554273387e-2);

    _generateCaseSix(qpoints, 254, 0.1233548532583327, 0.4127724083168531, 0.3392312205006170e-2);

    return qpoints;
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::_generateQuadratureWith434Points() const
{
    CMemBlock2D<double> qpoints(_nAngularPoints, 4);

    _generateCaseOne(qpoints, 0, 0.5265897968224436e-3);

    _generateCaseTwo(qpoints, 6, 0.2548219972002607e-2);

    _generateCaseThree(qpoints, 18, 0.2512317418927307e-2);

    _generateCaseFour(qpoints, 26, 0.6909346307509111, 0.2530403801186355e-2);

    _generateCaseFour(qpoints, 50, 0.1774836054609158, 0.2014279020918528e-2);

    _generateCaseFour(qpoints, 74, 0.4914342637784746, 0.2501725168402936e-2);

    _generateCaseFour(qpoints, 98, 0.6456664707424256, 0.2513267174597564e-2);

    _generateCaseFour(qpoints, 122, 0.2861289010307638, 0.2302694782227416e-2);

    _generateCaseFour(qpoints, 146, 0.7568084367178018e-1, 0.1462495621594614e-2);

    _generateCaseFour(qpoints, 170, 0.3927259763368002, 0.2445373437312980e-2);

    _generateCaseFive(qpoints, 194, 0.8818132877794288, 0.2417442375638981e-2);

    _generateCaseFive(qpoints, 218, 0.9776428111182649, 0.1910951282179532e-2);

    _generateCaseSix(qpoints, 242, 0.2054823696403044, 0.8689460322872412, 0.2416930044324775e-2);

    _generateCaseSix(qpoints, 290, 0.5905157048925271, 0.7999278543857286, 0.2512236854563495e-2);

    _generateCaseSix(qpoints, 338, 0.5550152361076807, 0.7717462626915901, 0.2496644054553086e-2);

    _generateCaseSix(qpoints, 386, 0.9371809858553722, 0.3344363145343455, 0.2236607760437849e-2);

    return qpoints;
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::_generateQuadratureWith590Points() const
{
    CMemBlock2D<double> qpoints(_nAngularPoints, 4);

    _generateCaseOne(qpoints, 0, 0.3095121295306187e-3);

    _generateCaseThree(qpoints, 6, 0.1852379698597489e-2);

    _generateCaseFour(qpoints, 14, 0.7040954938227469, 0.1871790639277744e-2);

    _generateCaseFour(qpoints, 38, 0.6807744066455243, 0.1858812585438317e-2);

    _generateCaseFour(qpoints, 62, 0.6372546939258752, 0.1852028828296213e-2);

    _generateCaseFour(qpoints, 86, 0.5044419707800358, 0.1846715956151242e-2);

    _generateCaseFour(qpoints, 110, 0.4215761784010967, 0.1818471778162769e-2);

    _generateCaseFour(qpoints, 134, 0.3317920736472123, 0.1749564657281154e-2);

    _generateCaseFour(qpoints, 158, 0.2384736701421887, 0.1617210647254411e-2);

    _generateCaseFour(qpoints, 182, 0.1459036449157763, 0.1384737234851692e-2);

    _generateCaseFour(qpoints, 206, 0.6095034115507196e-1, 0.9764331165051050e-3);

    _generateCaseFive(qpoints, 230, 0.6116843442009876, 0.1857161196774078e-2);

    _generateCaseFive(qpoints, 254, 0.3964755348199858, 0.1705153996395864e-2);

    _generateCaseFive(qpoints, 278, 0.1724782009907724, 0.1300321685886048e-2);

    _generateCaseSix(qpoints, 302, 0.5610263808622060, 0.3518280927733519, 0.1842866472905286e-2);

    _generateCaseSix(qpoints, 350, 0.4742392842551980, 0.2634716655937950, 0.1802658934377451e-2);

    _generateCaseSix(qpoints, 398, 0.5984126497885380, 0.1816640840360209, 0.1849830560443660e-2);

    _generateCaseSix(qpoints, 446, 0.3791035407695563, 0.1720795225656878, 0.1713904507106709e-2);

    _generateCaseSix(qpoints, 494, 0.2778673190586244, 0.8213021581932511e-1, 0.1555213603396808e-2);

    _generateCaseSix(qpoints, 542, 0.5033564271075117, 0.8999205842074875e-1, 0.1802239128008525e-2);

    return qpoints;
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::_generateQuadratureWith770Points() const
{
    CMemBlock2D<double> qpoints(_nAngularPoints, 4);

    _generateCaseOne(qpoints, 0, 0.2192942088181184e-3);

    _generateCaseTwo(qpoints, 6, 0.1436433617319080e-2);

    _generateCaseThree(qpoints, 18, 0.1421940344335877e-2);

    _generateCaseFour(qpoints, 26, 0.5087204410502360e-1, 0.6798123511050502e-3);

    _generateCaseFour(qpoints, 50, 0.1228198790178831, 0.9913184235294912e-3);

    _generateCaseFour(qpoints, 74, 0.2026890814408786, 0.1180207833238949e-2);

    _generateCaseFour(qpoints, 98, 0.2847745156464294, 0.1296599602080921e-2);

    _generateCaseFour(qpoints, 122, 0.3656719078978026, 0.1365871427428316e-2);

    _generateCaseFour(qpoints, 146, 0.4428264886713469, 0.1402988604775325e-2);

    _generateCaseFour(qpoints, 170, 0.5140619627249735, 0.1418645563595609e-2);

    _generateCaseFour(qpoints, 194, 0.6306401219166803, 0.1421376741851662e-2);

    _generateCaseFour(qpoints, 218, 0.6716883332022612, 0.1423996475490962e-2);

    _generateCaseFour(qpoints, 242, 0.6979792685336881, 0.1431554042178567e-2);

    _generateCaseFive(qpoints, 266, 0.1446865674195309, 0.9254401499865368e-3);

    _generateCaseFive(qpoints, 290, 0.3390263475411216, 0.1250239995053509e-2);

    _generateCaseFive(qpoints, 314, 0.5335804651263506, 0.1394365843329230e-2);

    _generateCaseSix(qpoints, 338, 0.6944024393349413e-1, 0.2355187894242326, 0.1127089094671749e-2);

    _generateCaseSix(qpoints, 386, 0.2269004109529460, 0.4102182474045730, 0.1345753760910670e-2);

    _generateCaseSix(qpoints, 434, 0.8025574607775339e-1, 0.6214302417481605, 0.1424957283316783e-2);

    _generateCaseSix(qpoints, 482, 0.1467999527896572, 0.3245284345717394, 0.1261523341237750e-2);

    _generateCaseSix(qpoints, 530, 0.1571507769824727, 0.5224482189696630, 0.1392547106052696e-2);

    _generateCaseSix(qpoints, 578, 0.2365702993157246, 0.6017546634089558, 0.1418761677877656e-2);

    _generateCaseSix(qpoints, 626, 0.7714815866765732e-1, 0.4346575516141163, 0.1338366684479554e-2);

    _generateCaseSix(qpoints, 674, 0.3062936666210730, 0.4908826589037616, 0.1393700862676131e-2);

    _generateCaseSix(qpoints, 722, 0.3822477379524787, 0.5648768149099500, 0.1415914757466932e-2);

    return qpoints;
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::_generateQuadratureWith974Points() const
{
    CMemBlock2D<double> qpoints(_nAngularPoints, 4);

    _generateCaseOne(qpoints, 0, 0.1438294190527431e-3);

    _generateCaseThree(qpoints, 6, 0.1125772288287004e-2);

    _generateCaseFour(qpoints, 14, 0.4292963545341347e-1, 0.4948029341949241e-3);

    _generateCaseFour(qpoints, 38, 0.1051426854086404, 0.7357990109125470e-3);

    _generateCaseFour(qpoints, 62, 0.1750024867623087, 0.8889132771304384e-3);

    _generateCaseFour(qpoints, 86, 0.2477653379650257, 0.9888347838921435e-3);

    _generateCaseFour(qpoints, 110, 0.3206567123955957, 0.1053299681709471e-2);

    _generateCaseFour(qpoints, 134, 0.3916520749849983, 0.1092778807014578e-2);

    _generateCaseFour(qpoints, 158, 0.4590825874187624, 0.1114389394063227e-2);

    _generateCaseFour(qpoints, 182, 0.5214563888415861, 0.1123724788051555e-2);

    _generateCaseFour(qpoints, 206, 0.6253170244654199, 0.1125239325243814e-2);

    _generateCaseFour(qpoints, 230, 0.6637926744523170, 0.1126153271815905e-2);

    _generateCaseFour(qpoints, 254, 0.6910410398498301, 0.1130286931123841e-2);

    _generateCaseFour(qpoints, 278, 0.7052907007457760, 0.1134986534363955e-2);

    _generateCaseFive(qpoints, 302, 0.1236686762657990, 0.6823367927109931e-3);

    _generateCaseFive(qpoints, 326, 0.2940777114468387, 0.9454158160447096e-3);

    _generateCaseFive(qpoints, 350, 0.4697753849207649, 0.1074429975385679e-2);

    _generateCaseFive(qpoints, 374, 0.6334563241139567, 0.1129300086569132e-2);

    _generateCaseSix(qpoints, 398, 0.5974048614181342e-1, 0.2029128752777523, 0.8436884500901954e-3);

    _generateCaseSix(qpoints, 446, 0.1375760408473636, 0.4602621942484054, 0.1075255720448885e-2);

    _generateCaseSix(qpoints, 494, 0.3391016526336286, 0.5030673999662036, 0.1108577236864462e-2);

    _generateCaseSix(qpoints, 542, 0.1271675191439820, 0.2817606422442134, 0.9566475323783357e-3);

    _generateCaseSix(qpoints, 590, 0.2693120740413512, 0.4331561291720157, 0.1080663250717391e-2);

    _generateCaseSix(qpoints, 638, 0.1419786452601918, 0.6256167358580814, 0.1126797131196295e-2);

    _generateCaseSix(qpoints, 686, 0.6709284600738255e-1, 0.3798395216859157, 0.1022568715358061e-2);

    _generateCaseSix(qpoints, 734, 0.7057738183256172e-1, 0.5517505421423520, 0.1108960267713108e-2);

    _generateCaseSix(qpoints, 782, 0.2783888477882155, 0.6029619156159187, 0.1122790653435766e-2);

    _generateCaseSix(qpoints, 830, 0.1979578938917407, 0.3589606329589096, 0.1032401847117460e-2);

    _generateCaseSix(qpoints, 878, 0.2087307061103274, 0.5348666438135476, 0.1107249382283854e-2);

    _generateCaseSix(qpoints, 926, 0.4055122137872836, 0.5674997546074373, 0.1121780048519972e-2);

    return qpoints;
}

CMemBlock2D<double>
CLebedevLaikovQuadrature::_generateQuadratureWith2030Points() const
{
    CMemBlock2D<double> qpoints(_nAngularPoints, 4);

    _generateCaseOne(qpoints, 0, 0.4656031899197431e-4);

    _generateCaseThree(qpoints, 6, 0.5421549195295507e-3);

    _generateCaseFour(qpoints, 14, 0.2540835336814348e-1, 0.1778522133346553e-3);

    _generateCaseFour(qpoints, 38, 0.6399322800504915e-1, 0.2811325405682796e-3);

    _generateCaseFour(qpoints, 62, 0.1088269469804125, 0.3548896312631459e-3);

    _generateCaseFour(qpoints, 86, 0.1570670798818287, 0.4090310897173364e-3);

    _generateCaseFour(qpoints, 110, 0.2071163932282514, 0.4493286134169965e-3);

    _generateCaseFour(qpoints, 134, 0.2578914044450844, 0.4793728447962723e-3);

    _generateCaseFour(qpoints, 158, 0.3085687558169623, 0.5015415319164265e-3);

    _generateCaseFour(qpoints, 182, 0.3584719706267024, 0.5175127372677937e-3);

    _generateCaseFour(qpoints, 206, 0.4070135594428709, 0.5285522262081019e-3);

    _generateCaseFour(qpoints, 230, 0.4536618626222638, 0.5356832703713962e-3);

    _generateCaseFour(qpoints, 254, 0.4979195686463577, 0.5397914736175170e-3);

    _generateCaseFour(qpoints, 278, 0.5393075111126999, 0.5416899441599930e-3);

    _generateCaseFour(qpoints, 302, 0.6115617676843916, 0.5419308476889938e-3);

    _generateCaseFour(qpoints, 326, 0.6414308435160159, 0.5416936902030596e-3);

    _generateCaseFour(qpoints, 350, 0.6664099412721607, 0.5419544338703164e-3);

    _generateCaseFour(qpoints, 374, 0.6859161771214913, 0.5428983656630975e-3);

    _generateCaseFour(qpoints, 398, 0.6993625593503890, 0.5442286500098193e-3);

    _generateCaseFour(qpoints, 422, 0.7062393387719380, 0.5452250345057301e-3);

    _generateCaseFive(qpoints, 446, 0.7479028168349763e-1, 0.2568002497728530e-3);

    _generateCaseFive(qpoints, 470, 0.1848951153969366, 0.3827211700292145e-3);

    _generateCaseFive(qpoints, 494, 0.3059529066581305, 0.4579491561917824e-3);

    _generateCaseFive(qpoints, 518, 0.4285556101021362, 0.5042003969083574e-3);

    _generateCaseFive(qpoints, 542, 0.5468758653496526, 0.5312708889976025e-3);

    _generateCaseFive(qpoints, 566, 0.6565821978343439, 0.5438401790747117e-3);

    _generateCaseSix(qpoints, 590, 0.1253901572367117, 0.3681917226439641e-1, 0.3316041873197344e-3);

    _generateCaseSix(qpoints, 638, 0.1775721510383941, 0.7982487607213301e-1, 0.3899113567153771e-3);

    _generateCaseSix(qpoints, 686, 0.2305693358216114, 0.1264640966592335, 0.4343343327201309e-3);

    _generateCaseSix(qpoints, 734, 0.2836502845992063, 0.1751585683418957, 0.4679415262318919e-3);

    _generateCaseSix(qpoints, 782, 0.3361794746232590, 0.2247995907632670, 0.4930847981631031e-3);

    _generateCaseSix(qpoints, 830, 0.3875979172264824, 0.2745299257422246, 0.5115031867540091e-3);

    _generateCaseSix(qpoints, 878, 0.4374019316999074, 0.3236373482441118, 0.5245217148457367e-3);

    _generateCaseSix(qpoints, 926, 0.4851275843340022, 0.3714967859436741, 0.5332041499895321e-3);

    _generateCaseSix(qpoints, 974, 0.5303391803806868, 0.4175353646321745, 0.5384583126021542e-3);

    _generateCaseSix(qpoints, 1022, 0.5726197380596287, 0.4612084406355461, 0.5411067210798852e-3);

    _generateCaseSix(qpoints, 1070, 0.2431520732564863, 0.4258040133043952e-1, 0.4259797391468714e-3);

    _generateCaseSix(qpoints, 1118, 0.3002096800895869, 0.8869424306722721e-1, 0.4604931368460021e-3);

    _generateCaseSix(qpoints, 1166, 0.3558554457457432, 0.1368811706510655, 0.4871814878255202e-3);

    _generateCaseSix(qpoints, 1214, 0.4097782537048887, 0.1860739985015033, 0.5072242910074885e-3);

    _generateCaseSix(qpoints, 1262, 0.4616337666067458, 0.2354235077395853, 0.5217069845235350e-3);

    _generateCaseSix(qpoints, 1310, 0.5110707008417874, 0.2842074921347011, 0.5315785966280310e-3);

    _generateCaseSix(qpoints, 1358, 0.5577415286163795, 0.3317784414984102, 0.5376833708758905e-3);

    _generateCaseSix(qpoints, 1406, 0.6013060431366950, 0.3775299002040700, 0.5408032092069521e-3);

    _generateCaseSix(qpoints, 1454, 0.3661596767261781, 0.4599367887164592e-1, 0.4842744917904866e-3);

    _generateCaseSix(qpoints, 1502, 0.4237633153506581, 0.9404893773654421e-1, 0.5048926076188130e-3);

    _generateCaseSix(qpoints, 1550, 0.4786328454658452, 0.1431377109091971, 0.5202607980478373e-3);

    _generateCaseSix(qpoints, 1598, 0.5305702076789774, 0.1924186388843570, 0.5309932388325743e-3);

    _generateCaseSix(qpoints, 1646, 0.5793436224231788, 0.2411590944775190, 0.5377419770895208e-3);

    _generateCaseSix(qpoints, 1694, 0.6247069017094747, 0.2886871491583605, 0.5411696331677717e-3);

    _generateCaseSix(qpoints, 1742, 0.4874315552535204, 0.4804978774953206e-1, 0.5197996293282420e-3);

    _generateCaseSix(qpoints, 1790, 0.5427337322059053, 0.9716857199366665e-1, 0.5311120836622945e-3);

    _generateCaseSix(qpoints, 1838, 0.5943493747246700, 0.1465205839795055, 0.5384309319956951e-3);

    _generateCaseSix(qpoints, 1886, 0.6421314033564943, 0.1953579449803574, 0.5421859504051886e-3);

    _generateCaseSix(qpoints, 1934, 0.6020628374713980, 0.4916375015738108e-1, 0.5390948355046314e-3);

    _generateCaseSix(qpoints, 1982, 0.6529222529856881, 0.9861621540127005e-1, 0.5433312705027845e-3);

    return qpoints;
}

void
CLebedevLaikovQuadrature::_generateCaseOne(CMemBlock2D<double>& gridPoints, const int32_t offset, const double weight) const
{
    auto x = gridPoints.data(0, offset);

    auto y = gridPoints.data(1, offset);

    auto z = gridPoints.data(2, offset);

    auto w = gridPoints.data(3, offset);

    x[0] = 1.0;
    y[0] = 0.0;
    z[0] = 0.0;
    w[0] = weight;

    x[1] = -1.0;
    y[1] = 0.0;
    z[1] = 0.0;
    w[1] = weight;

    x[2] = 0.0;
    y[2] = 1.0;
    z[2] = 0.0;
    w[2] = weight;

    x[3] = 0.0;
    y[3] = -1.0;
    z[3] = 0.0;
    w[3] = weight;

    x[4] = 0.0;
    y[4] = 0.0;
    z[4] = 1.0;
    w[4] = weight;

    x[5] = 0.0;
    y[5] = 0.0;
    z[5] = -1.0;
    w[5] = weight;
}

void
CLebedevLaikovQuadrature::_generateCaseTwo(CMemBlock2D<double>& gridPoints, const int32_t offset, const double weight) const
{
    auto x = gridPoints.data(0, offset);

    auto y = gridPoints.data(1, offset);

    auto z = gridPoints.data(2, offset);

    auto w = gridPoints.data(3, offset);

    double a = 1.0 / std::sqrt(2.0);

    x[0] = 0.0;
    y[0] = a;
    z[0] = a;
    w[0] = weight;

    x[1] = 0.0;
    y[1] = a;
    z[1] = -a;
    w[1] = weight;

    x[2] = 0.0;
    y[2] = -a;
    z[2] = a;
    w[2] = weight;

    x[3] = 0.0;
    y[3] = -a;
    z[3] = -a;
    w[3] = weight;

    x[4] = a;
    y[4] = 0.0;
    z[4] = a;
    w[4] = weight;

    x[5] = a;
    y[5] = 0.0;
    z[5] = -a;
    w[5] = weight;

    x[6] = -a;
    y[6] = 0.0;
    z[6] = a;
    w[6] = weight;

    x[7] = -a;
    y[7] = 0.0;
    z[7] = -a;
    w[7] = weight;

    x[8] = a;
    y[8] = a;
    z[8] = 0.0;
    w[8] = weight;

    x[9] = a;
    y[9] = -a;
    z[9] = 0.0;
    w[9] = weight;

    x[10] = -a;
    y[10] = a;
    z[10] = 0.0;
    w[10] = weight;

    x[11] = -a;
    y[11] = -a;
    z[11] = 0.0;
    w[11] = weight;
}

void
CLebedevLaikovQuadrature::_generateCaseThree(CMemBlock2D<double>& gridPoints, const int32_t offset, const double weight) const
{
    auto x = gridPoints.data(0, offset);

    auto y = gridPoints.data(1, offset);

    auto z = gridPoints.data(2, offset);

    auto w = gridPoints.data(3, offset);

    double a = 1.0 / std::sqrt(3.0);

    x[0] = a;
    y[0] = a;
    z[0] = a;
    w[0] = weight;

    x[1] = a;
    y[1] = a;
    z[1] = -a;
    w[1] = weight;

    x[2] = a;
    y[2] = -a;
    z[2] = a;
    w[2] = weight;

    x[3] = a;
    y[3] = -a;
    z[3] = -a;
    w[3] = weight;

    x[4] = -a;
    y[4] = a;
    z[4] = a;
    w[4] = weight;

    x[5] = -a;
    y[5] = a;
    z[5] = -a;
    w[5] = weight;

    x[6] = -a;
    y[6] = -a;
    z[6] = a;
    w[6] = weight;

    x[7] = -a;
    y[7] = -a;
    z[7] = -a;
    w[7] = weight;
}

void
CLebedevLaikovQuadrature::_generateCaseFour(CMemBlock2D<double>& gridPoints, const int32_t offset, const double factor, const double weight) const
{
    auto x = gridPoints.data(0, offset);

    auto y = gridPoints.data(1, offset);

    auto z = gridPoints.data(2, offset);

    auto w = gridPoints.data(3, offset);

    auto a = factor;

    auto b = std::sqrt(1.0 - 2.0 * factor * factor);

    x[0] = a;
    y[0] = a;
    z[0] = b;
    w[0] = weight;

    x[1] = a;
    y[1] = a;
    z[1] = -b;
    w[1] = weight;

    x[2] = a;
    y[2] = -a;
    z[2] = b;
    w[2] = weight;

    x[3] = a;
    y[3] = -a;
    z[3] = -b;
    w[3] = weight;

    x[4] = -a;
    y[4] = a;
    z[4] = b;
    w[4] = weight;

    x[5] = -a;
    y[5] = a;
    z[5] = -b;
    w[5] = weight;

    x[6] = -a;
    y[6] = -a;
    z[6] = b;
    w[6] = weight;

    x[7] = -a;
    y[7] = -a;
    z[7] = -b;
    w[7] = weight;

    x[8] = a;
    y[8] = b;
    z[8] = a;
    w[8] = weight;

    x[9] = a;
    y[9] = -b;
    z[9] = a;
    w[9] = weight;

    x[10] = a;
    y[10] = b;
    z[10] = -a;
    w[10] = weight;

    x[11] = a;
    y[11] = -b;
    z[11] = -a;
    w[11] = weight;

    x[12] = -a;
    y[12] = b;
    z[12] = a;
    w[12] = weight;

    x[13] = -a;
    y[13] = -b;
    z[13] = a;
    w[13] = weight;

    x[14] = -a;
    y[14] = b;
    z[14] = -a;
    w[14] = weight;

    x[15] = -a;
    y[15] = -b;
    z[15] = -a;
    w[15] = weight;

    x[16] = b;
    y[16] = a;
    z[16] = a;
    w[16] = weight;

    x[17] = -b;
    y[17] = a;
    z[17] = a;
    w[17] = weight;

    x[18] = b;
    y[18] = a;
    z[18] = -a;
    w[18] = weight;

    x[19] = -b;
    y[19] = a;
    z[19] = -a;
    w[19] = weight;

    x[20] = b;
    y[20] = -a;
    z[20] = a;
    w[20] = weight;

    x[21] = -b;
    y[21] = -a;
    z[21] = a;
    w[21] = weight;

    x[22] = b;
    y[22] = -a;
    z[22] = -a;
    w[22] = weight;

    x[23] = -b;
    y[23] = -a;
    z[23] = -a;
    w[23] = weight;
}

void
CLebedevLaikovQuadrature::_generateCaseFive(CMemBlock2D<double>& gridPoints, const int32_t offset, const double factor, const double weight) const
{
    auto x = gridPoints.data(0, offset);

    auto y = gridPoints.data(1, offset);

    auto z = gridPoints.data(2, offset);

    auto w = gridPoints.data(3, offset);

    auto a = factor;

    auto b = std::sqrt(1.0 - factor * factor);

    x[0] = a;
    y[0] = b;
    z[0] = 0.0;
    w[0] = weight;

    x[1] = a;
    y[1] = -b;
    z[1] = 0.0;
    w[1] = weight;

    x[2] = -a;
    y[2] = b;
    z[2] = 0.0;
    w[2] = weight;

    x[3] = -a;
    y[3] = -b;
    z[3] = 0.0;
    w[3] = weight;

    x[4] = b;
    y[4] = a;
    z[4] = 0.0;
    w[4] = weight;

    x[5] = b;
    y[5] = -a;
    z[5] = 0.0;
    w[5] = weight;

    x[6] = -b;
    y[6] = a;
    z[6] = 0.0;
    w[6] = weight;

    x[7] = -b;
    y[7] = -a;
    z[7] = 0.0;
    w[7] = weight;

    x[8] = a;
    y[8] = 0.0;
    z[8] = b;
    w[8] = weight;

    x[9] = a;
    y[9] = 0.0;
    z[9] = -b;
    w[9] = weight;

    x[10] = -a;
    y[10] = 0.0;
    z[10] = b;
    w[10] = weight;

    x[11] = -a;
    y[11] = 0.0;
    z[11] = -b;
    w[11] = weight;

    x[12] = b;
    y[12] = 0.0;
    z[12] = a;
    w[12] = weight;

    x[13] = b;
    y[13] = 0.0;
    z[13] = -a;
    w[13] = weight;

    x[14] = -b;
    y[14] = 0.0;
    z[14] = a;
    w[14] = weight;

    x[15] = -b;
    y[15] = 0.0;
    z[15] = -a;
    w[15] = weight;

    x[16] = 0.0;
    y[16] = a;
    z[16] = b;
    w[16] = weight;

    x[17] = 0.0;
    y[17] = a;
    z[17] = -b;
    w[17] = weight;

    x[18] = 0.0;
    y[18] = -a;
    z[18] = b;
    w[18] = weight;

    x[19] = 0.0;
    y[19] = -a;
    z[19] = -b;
    w[19] = weight;

    x[20] = 0.0;
    y[20] = b;
    z[20] = a;
    w[20] = weight;

    x[21] = 0.0;
    y[21] = b;
    z[21] = -a;
    w[21] = weight;

    x[22] = 0.0;
    y[22] = -b;
    z[22] = a;
    w[22] = weight;

    x[23] = 0.0;
    y[23] = -b;
    z[23] = -a;
    w[23] = weight;
}

void
CLebedevLaikovQuadrature::_generateCaseSix(CMemBlock2D<double>& gridPoints,
                                           const int32_t        offset,
                                           const double         factorA,
                                           const double         factorB,
                                           const double         weight) const
{
    auto x = gridPoints.data(0, offset);

    auto y = gridPoints.data(1, offset);

    auto z = gridPoints.data(2, offset);

    auto w = gridPoints.data(3, offset);

    auto a = factorA;

    auto b = factorB;

    auto c = std::sqrt(1.0 - factorA * factorA - factorB * factorB);

    x[0] = a;
    y[0] = b;
    z[0] = c;
    w[0] = weight;

    x[1] = a;
    y[1] = b;
    z[1] = -c;
    w[1] = weight;

    x[2] = a;
    y[2] = -b;
    z[2] = c;
    w[2] = weight;

    x[3] = a;
    y[3] = -b;
    z[3] = -c;
    w[3] = weight;

    x[4] = -a;
    y[4] = b;
    z[4] = c;
    w[4] = weight;

    x[5] = -a;
    y[5] = b;
    z[5] = -c;
    w[5] = weight;

    x[6] = -a;
    y[6] = -b;
    z[6] = c;
    w[6] = weight;

    x[7] = -a;
    y[7] = -b;
    z[7] = -c;
    w[7] = weight;

    x[8] = a;
    y[8] = c;
    z[8] = b;
    w[8] = weight;

    x[9] = a;
    y[9] = c;
    z[9] = -b;
    w[9] = weight;

    x[10] = a;
    y[10] = -c;
    z[10] = b;
    w[10] = weight;

    x[11] = a;
    y[11] = -c;
    z[11] = -b;
    w[11] = weight;

    x[12] = -a;
    y[12] = c;
    z[12] = b;
    w[12] = weight;

    x[13] = -a;
    y[13] = c;
    z[13] = -b;
    w[13] = weight;

    x[14] = -a;
    y[14] = -c;
    z[14] = b;
    w[14] = weight;

    x[15] = -a;
    y[15] = -c;
    z[15] = -b;
    w[15] = weight;

    x[16] = b;
    y[16] = a;
    z[16] = c;
    w[16] = weight;

    x[17] = b;
    y[17] = a;
    z[17] = -c;
    w[17] = weight;

    x[18] = b;
    y[18] = -a;
    z[18] = c;
    w[18] = weight;

    x[19] = b;
    y[19] = -a;
    z[19] = -c;
    w[19] = weight;

    x[20] = -b;
    y[20] = a;
    z[20] = c;
    w[20] = weight;

    x[21] = -b;
    y[21] = a;
    z[21] = -c;
    w[21] = weight;

    x[22] = -b;
    y[22] = -a;
    z[22] = c;
    w[22] = weight;

    x[23] = -b;
    y[23] = -a;
    z[23] = -c;
    w[23] = weight;

    x[24] = b;
    y[24] = c;
    z[24] = a;
    w[24] = weight;

    x[25] = b;
    y[25] = c;
    z[25] = -a;
    w[25] = weight;

    x[26] = b;
    y[26] = -c;
    z[26] = a;
    w[26] = weight;

    x[27] = b;
    y[27] = -c;
    z[27] = -a;
    w[27] = weight;

    x[28] = -b;
    y[28] = c;
    z[28] = a;
    w[28] = weight;

    x[29] = -b;
    y[29] = c;
    z[29] = -a;
    w[29] = weight;

    x[30] = -b;
    y[30] = -c;
    z[30] = a;
    w[30] = weight;

    x[31] = -b;
    y[31] = -c;
    z[31] = -a;
    w[31] = weight;

    x[32] = c;
    y[32] = a;
    z[32] = b;
    w[32] = weight;

    x[33] = c;
    y[33] = a;
    z[33] = -b;
    w[33] = weight;

    x[34] = c;
    y[34] = -a;
    z[34] = b;
    w[34] = weight;

    x[35] = c;
    y[35] = -a;
    z[35] = -b;
    w[35] = weight;

    x[36] = -c;
    y[36] = a;
    z[36] = b;
    w[36] = weight;

    x[37] = -c;
    y[37] = a;
    z[37] = -b;
    w[37] = weight;

    x[38] = -c;
    y[38] = -a;
    z[38] = b;
    w[38] = weight;

    x[39] = -c;
    y[39] = -a;
    z[39] = -b;
    w[39] = weight;

    x[40] = c;
    y[40] = b;
    z[40] = a;
    w[40] = weight;

    x[41] = c;
    y[41] = b;
    z[41] = -a;
    w[41] = weight;

    x[42] = c;
    y[42] = -b;
    z[42] = a;
    w[42] = weight;

    x[43] = c;
    y[43] = -b;
    z[43] = -a;
    w[43] = weight;

    x[44] = -c;
    y[44] = b;
    z[44] = a;
    w[44] = weight;

    x[45] = -c;
    y[45] = b;
    z[45] = -a;
    w[45] = weight;

    x[46] = -c;
    y[46] = -b;
    z[46] = a;
    w[46] = weight;

    x[47] = -c;
    y[47] = -b;
    z[47] = -a;
    w[47] = weight;
}
