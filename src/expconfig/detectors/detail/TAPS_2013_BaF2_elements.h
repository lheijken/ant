#include "detectors/TAPS.h"
const std::vector<ant::expconfig::detector::TAPS_2013::BaF2_Element_t>
ant::expconfig::detector::TAPS_2013::BaF2_elements_init = {
  // element, {xy}, TAC, LG, SG, LGS, SGS, {neigbours}
  { 12, {-15.588,   9.000}, 20011, 20018, 20019, 20016, 20017, {5,379,383,15,16,13,7}},
  { 13, {-15.588,   3.000}, 20040, 20034, 20035, 20032, 20033, {9,7,5,12,16,17,14,11}},
  { 14, {-15.588,  -3.000}, 20041, 20038, 20039, 20036, 20037, {78,11,9,13,17,18,85,80}},
  { 15, {-20.785,  12.000}, 20042, 20046, 20047, 20044, 20045, {12,383,388,19,20,16}},
  { 16, {-20.785,   6.000}, 20043, 20050, 20051, 20048, 20049, {13,12,15,20,21,17}},
  { 17, {-20.785,   0.000}, 20072, 20066, 20067, 20064, 20065, {13,16,21,22,18,14}},
  { 18, {-20.785,  -6.000}, 20073, 20070, 20071, 20068, 20069, {88,85,14,17,22,23}},
  { 19, {-25.981,  15.000}, 20074, 20078, 20079, 20076, 20077, {20,15,388,394,24,25}},
  { 20, {-25.981,   9.000}, 20075, 20082, 20083, 20080, 20081, {21,16,15,19,25,26}},
  { 21, {-25.981,   3.000}, 20104, 20098, 20099, 20096, 20097, {22,17,16,20,26,27}},
  { 22, {-25.981,  -3.000}, 20105, 20102, 20103, 20100, 20101, {23,18,17,21,27,28}},
  { 23, {-25.981,  -9.000}, 20106, 20110, 20111, 20108, 20109, {92,88,18,22,28,29}},
  { 24, {-31.177,  18.000}, 20107, 20114, 20115, 20112, 20113, {25,19,394,401,30,31}},
  { 25, {-31.177,  12.000}, 20136, 20130, 20131, 20128, 20129, {26,20,19,24,31,32}},
  { 26, {-31.177,   6.000}, 20137, 20134, 20135, 20132, 20133, {27,21,20,25,32,33}},
  { 27, {-31.177,   0.000}, 20138, 20142, 20143, 20140, 20141, {28,22,21,26,33,34}},
  { 28, {-31.177,  -6.000}, 20139, 20146, 20147, 20144, 20145, {29,23,22,27,34,35}},
  { 29, {-31.177, -12.000}, 20168, 20162, 20163, 20160, 20161, {97,92,23,28,35,36}},
  { 30, {-36.373,  21.000}, 20169, 20166, 20167, 20164, 20165, {31,24,401,409,37,38}},
  { 31, {-36.373,  15.000}, 20170, 20174, 20175, 20172, 20173, {32,25,24,30,38,39}},
  { 32, {-36.373,   9.000}, 20171, 20178, 20179, 20176, 20177, {33,26,25,31,39,40}},
  { 33, {-36.373,   3.000}, 20200, 20194, 20195, 20192, 20193, {34,27,26,32,40,41}},
  { 34, {-36.373,  -3.000}, 20201, 20198, 20199, 20196, 20197, {35,28,27,33,41,42}},
  { 35, {-36.373,  -9.000}, 20202, 20206, 20207, 20204, 20205, {36,29,28,34,42,43}},
  { 36, {-36.373, -15.000}, 20203, 20210, 20211, 20208, 20209, {103,97,29,35,43,44}},
  { 37, {-41.569,  24.000}, 20232, 20226, 20227, 20224, 20225, {38,30,409,418,45,46}},
  { 38, {-41.569,  18.000}, 20233, 20230, 20231, 20228, 20229, {39,31,30,37,46,47}},
  { 39, {-41.569,  12.000}, 20234, 20238, 20239, 20236, 20237, {40,32,31,38,47,48}},
  { 40, {-41.569,   6.000}, 20235, 20242, 20243, 20240, 20241, {41,33,32,39,48,49}},
  { 41, {-41.569,   0.000}, 20264, 20258, 20259, 20256, 20257, {42,34,33,40,49,50}},
  { 42, {-41.569,  -6.000}, 20265, 20262, 20263, 20260, 20261, {43,35,34,41,50,51}},
  { 43, {-41.569, -12.000}, 20266, 20270, 20271, 20268, 20269, {44,36,35,42,51,52}},
  { 44, {-41.569, -18.000}, 20267, 20274, 20275, 20272, 20273, {110,103,36,43,52,53}},
  { 45, {-46.765,  27.000}, 20296, 20290, 20291, 20288, 20289, {46,37,418,428,54,55}},
  { 46, {-46.765,  21.000}, 20297, 20294, 20295, 20292, 20293, {47,38,37,45,55,56}},
  { 47, {-46.765,  15.000}, 20298, 20302, 20303, 20300, 20301, {48,39,38,46,56,57}},
  { 48, {-46.765,   9.000}, 20299, 20306, 20307, 20304, 20305, {49,40,39,47,57,58}},
  { 49, {-46.765,   3.000}, 20328, 20322, 20323, 20320, 20321, {50,41,40,48,58,59}},
  { 50, {-46.765,  -3.000}, 20329, 20326, 20327, 20324, 20325, {51,42,41,49,59,60}},
  { 51, {-46.765,  -9.000}, 20330, 20334, 20335, 20332, 20333, {52,43,42,50,60,61}},
  { 52, {-46.765, -15.000}, 20331, 20338, 20339, 20336, 20337, {53,44,43,51,61,62}},
  { 53, {-46.765, -21.000}, 20360, 20354, 20355, 20352, 20353, {118,110,44,52,62,63}},
  { 54, {-51.962,  30.000}, 20361, 20358, 20359, 20356, 20357, {55,45,428,437,64}},
  { 55, {-51.962,  24.000}, 20362, 20366, 20367, 20364, 20365, {56,46,45,54,64,65}},
  { 56, {-51.962,  18.000}, 20363, 20370, 20371, 20368, 20369, {57,47,46,55,65,66}},
  { 57, {-51.962,  12.000}, 20392, 20386, 20387, 20384, 20385, {58,48,47,56,66,67}},
  { 58, {-51.962,   6.000}, 20393, 20390, 20391, 20388, 20389, {59,49,48,57,67,68}},
  { 59, {-51.962,   0.000}, 20394, 20398, 20399, 20396, 20397, {60,50,49,58,68,69}},
  { 60, {-51.962,  -6.000}, 20395, 20402, 20403, 20400, 20401, {61,51,50,59,69,70}},
  { 61, {-51.962, -12.000}, 20424, 20418, 20419, 20416, 20417, {62,52,51,60,70,71}},
  { 62, {-51.962, -18.000}, 20425, 20422, 20423, 20420, 20421, {63,53,52,61,71,72}},
  { 63, {-51.962, -24.000}, 20426, 20430, 20431, 20428, 20429, {127,118,53,62,72}},
  { 64, {-57.158,  27.000}, 20427, 20434, 20435, 20432, 20433, {65,55,54}},
  { 65, {-57.158,  21.000}, 20456, 20450, 20451, 20448, 20449, {66,56,55,64}},
  { 66, {-57.158,  15.000}, 20457, 20454, 20455, 20452, 20453, {67,57,56,65}},
  { 67, {-57.158,   9.000}, 20458, 20462, 20463, 20460, 20461, {68,58,57,66}},
  { 68, {-57.158,   3.000}, 20459, 20466, 20467, 20464, 20465, {69,59,58,67}},
  { 69, {-57.158,  -3.000}, 20488, 20482, 20483, 20480, 20481, {70,60,59,68}},
  { 70, {-57.158,  -9.000}, 20489, 20486, 20487, 20484, 20485, {71,61,60,69}},
  { 71, {-57.158, -15.000}, 20490, 20494, 20495, 20492, 20493, {72,62,61,70}},
  { 72, {-57.158, -21.000}, 20491, 20498, 20499, 20496, 20497, {63,62,71}},
  { 85, {-15.588,  -9.000}, 21011, 21018, 21019, 21016, 21017, {80,78,14,18,88,89,86}},
  { 86, {-10.392, -12.000}, 21040, 21034, 21035, 21032, 21033, {84,82,79,80,85,89,90,87}},
  { 87, { -5.196, -15.000}, 21041, 21038, 21039, 21036, 21037, {153,151,83,84,86,90,91,158}},
  { 88, {-20.785, -12.000}, 21042, 21046, 21047, 21044, 21045, {93,89,85,18,23,92}},
  { 89, {-15.588, -15.000}, 21043, 21050, 21051, 21048, 21049, {94,90,86,85,88,93}},
  { 90, {-10.392, -18.000}, 21072, 21066, 21067, 21064, 21065, {95,91,87,86,89,94}},
  { 91, { -5.196, -21.000}, 21073, 21070, 21071, 21068, 21069, {96,161,158,87,90,95}},
  { 92, {-25.981, -15.000}, 21074, 21078, 21079, 21076, 21077, {98,93,88,23,29,97}},
  { 93, {-20.785, -18.000}, 21075, 21082, 21083, 21080, 21081, {99,94,89,88,92,98}},
  { 94, {-15.588, -21.000}, 21104, 21098, 21099, 21096, 21097, {100,95,90,89,93,99}},
  { 95, {-10.392, -24.000}, 21105, 21102, 21103, 21100, 21101, {101,96,91,90,94,100}},
  { 96, { -5.196, -27.000}, 21106, 21110, 21111, 21108, 21109, {102,165,161,91,95,101}},
  { 97, {-31.177, -18.000}, 21107, 21114, 21115, 21112, 21113, {104,98,92,29,36,103}},
  { 98, {-25.981, -21.000}, 21136, 21130, 21131, 21128, 21129, {105,99,93,92,97,104}},
  { 99, {-20.785, -24.000}, 21137, 21134, 21135, 21132, 21133, {106,100,94,93,98,105}},
  {100, {-15.588, -27.000}, 21138, 21142, 21143, 21140, 21141, {107,101,95,94,99,106}},
  {101, {-10.392, -30.000}, 21139, 21146, 21147, 21144, 21145, {108,102,96,95,100,107}},
  {102, { -5.196, -33.000}, 21168, 21162, 21163, 21160, 21161, {109,170,165,96,101,108}},
  {103, {-36.373, -21.000}, 21169, 21166, 21167, 21164, 21165, {111,104,97,36,44,110}},
  {104, {-31.177, -24.000}, 21170, 21174, 21175, 21172, 21173, {112,105,98,97,103,111}},
  {105, {-25.981, -27.000}, 21171, 21178, 21179, 21176, 21177, {113,106,99,98,104,112}},
  {106, {-20.785, -30.000}, 21200, 21194, 21195, 21192, 21193, {114,107,100,99,105,113}},
  {107, {-15.588, -33.000}, 21201, 21198, 21199, 21196, 21197, {115,108,101,100,106,114}},
  {108, {-10.392, -36.000}, 21202, 21206, 21207, 21204, 21205, {116,109,102,101,107,115}},
  {109, { -5.196, -39.000}, 21203, 21210, 21211, 21208, 21209, {117,176,170,102,108,116}},
  {110, {-41.569, -24.000}, 21232, 21226, 21227, 21224, 21225, {119,111,103,44,53,118}},
  {111, {-36.373, -27.000}, 21233, 21230, 21231, 21228, 21229, {120,112,104,103,110,119}},
  {112, {-31.177, -30.000}, 21234, 21238, 21239, 21236, 21237, {121,113,105,104,111,120}},
  {113, {-25.981, -33.000}, 21235, 21242, 21243, 21240, 21241, {122,114,106,105,112,121}},
  {114, {-20.785, -36.000}, 21264, 21258, 21259, 21256, 21257, {123,115,107,106,113,122}},
  {115, {-15.588, -39.000}, 21265, 21262, 21263, 21260, 21261, {124,116,108,107,114,123}},
  {116, {-10.392, -42.000}, 21266, 21270, 21271, 21268, 21269, {125,117,109,108,115,124}},
  {117, { -5.196, -45.000}, 21267, 21274, 21275, 21272, 21273, {126,183,176,109,116,125}},
  {118, {-46.765, -27.000}, 21296, 21290, 21291, 21288, 21289, {128,119,110,53,63,127}},
  {119, {-41.569, -30.000}, 21297, 21294, 21295, 21292, 21293, {129,120,111,110,118,128}},
  {120, {-36.373, -33.000}, 21298, 21302, 21303, 21300, 21301, {130,121,112,111,119,129}},
  {121, {-31.177, -36.000}, 21299, 21306, 21307, 21304, 21305, {131,122,113,112,120,130}},
  {122, {-25.981, -39.000}, 21328, 21322, 21323, 21320, 21321, {132,123,114,113,121,131}},
  {123, {-20.785, -42.000}, 21329, 21326, 21327, 21324, 21325, {133,124,115,114,122,132}},
  {124, {-15.588, -45.000}, 21330, 21334, 21335, 21332, 21333, {134,125,116,115,123,133}},
  {125, {-10.392, -48.000}, 21331, 21338, 21339, 21336, 21337, {135,126,117,116,124,134}},
  {126, { -5.196, -51.000}, 21360, 21354, 21355, 21352, 21353, {136,191,183,117,125,135}},
  {127, {-51.962, -30.000}, 21361, 21358, 21359, 21356, 21357, {137,128,118,63}},
  {128, {-46.765, -33.000}, 21362, 21366, 21367, 21364, 21365, {138,129,119,118,127,137}},
  {129, {-41.569, -36.000}, 21363, 21370, 21371, 21368, 21369, {139,130,120,119,128,138}},
  {130, {-36.373, -39.000}, 21392, 21386, 21387, 21384, 21385, {140,131,121,120,129,139}},
  {131, {-31.177, -42.000}, 21393, 21390, 21391, 21388, 21389, {141,132,122,121,130,140}},
  {132, {-25.981, -45.000}, 21394, 21398, 21399, 21396, 21397, {142,133,123,122,131,141}},
  {133, {-20.785, -48.000}, 21395, 21402, 21403, 21400, 21401, {143,134,124,123,132,142}},
  {134, {-15.588, -51.000}, 21424, 21418, 21419, 21416, 21417, {144,135,125,124,133,143}},
  {135, {-10.392, -54.000}, 21425, 21422, 21423, 21420, 21421, {145,136,126,125,134,144}},
  {136, { -5.196, -57.000}, 21426, 21430, 21431, 21428, 21429, {200,191,126,135,145}},
  {137, {-51.962, -36.000}, 21427, 21434, 21435, 21432, 21433, {138,128,127}},
  {138, {-46.765, -39.000}, 21456, 21450, 21451, 21448, 21449, {139,129,128,137}},
  {139, {-41.569, -42.000}, 21457, 21454, 21455, 21452, 21453, {140,130,129,138}},
  {140, {-36.373, -45.000}, 21458, 21462, 21463, 21460, 21461, {141,131,130,139}},
  {141, {-31.177, -48.000}, 21459, 21466, 21467, 21464, 21465, {142,132,131,140}},
  {142, {-25.981, -51.000}, 21488, 21482, 21483, 21480, 21481, {143,133,132,141}},
  {143, {-20.785, -54.000}, 21489, 21486, 21487, 21484, 21485, {144,134,133,142}},
  {144, {-15.588, -57.000}, 21490, 21494, 21495, 21492, 21493, {145,135,134,143}},
  {145, {-10.392, -60.000}, 21491, 21498, 21499, 21496, 21497, {136,135,144}},
  {158, {  0.000, -18.000}, 22011, 22018, 22019, 22016, 22017, {159,152,153,87,91,161,162}},
  {159, {  5.196, -15.000}, 22040, 22034, 22035, 22032, 22033, {160,156,157,150,152,158,162,163}},
  {160, { 10.392, -12.000}, 22041, 22038, 22039, 22036, 22037, {231,225,226,154,156,159,163,164}},
  {161, {  0.000, -24.000}, 22042, 22046, 22047, 22044, 22045, {165,166,162,158,91,96}},
  {162, {  5.196, -21.000}, 22043, 22050, 22051, 22048, 22049, {166,167,163,159,158,161}},
  {163, { 10.392, -18.000}, 22072, 22066, 22067, 22064, 22065, {167,168,164,160,159,162}},
  {164, { 15.588, -15.000}, 22073, 22070, 22071, 22068, 22069, {168,169,234,231,160,163}},
  {165, {  0.000, -30.000}, 22074, 22078, 22079, 22076, 22077, {170,171,166,161,96,102}},
  {166, {  5.196, -27.000}, 22075, 22082, 22083, 22080, 22081, {171,172,167,162,161,165}},
  {167, { 10.392, -24.000}, 22104, 22098, 22099, 22096, 22097, {172,173,168,163,162,166}},
  {168, { 15.588, -21.000}, 22105, 22102, 22103, 22100, 22101, {173,174,169,164,163,167}},
  {169, { 20.785, -18.000}, 22106, 22110, 22111, 22108, 22109, {174,175,238,234,164,168}},
  {170, {  0.000, -36.000}, 22107, 22114, 22115, 22112, 22113, {176,177,171,165,102,109}},
  {171, {  5.196, -33.000}, 22136, 22130, 22131, 22128, 22129, {177,178,172,166,165,170}},
  {172, { 10.392, -30.000}, 22137, 22134, 22135, 22132, 22133, {178,179,173,167,171,166}},
  {173, { 15.588, -27.000}, 22138, 22142, 22143, 22140, 22141, {179,180,174,168,172,167}},
  {174, { 20.785, -24.000}, 22139, 22146, 22147, 22144, 22145, {180,181,175,169,173,168}},
  {175, { 25.981, -21.000}, 22168, 22162, 22163, 22160, 22161, {181,182,243,238,169,174}},
  {176, {  0.000, -42.000}, 22169, 22166, 22167, 22164, 22165, {183,184,177,170,109,117}},
  {177, {  5.196, -39.000}, 22170, 22174, 22175, 22172, 22173, {184,185,178,171,170,176}},
  {178, { 10.392, -36.000}, 22171, 22178, 22179, 22176, 22177, {185,186,179,172,171,177}},
  {179, { 15.588, -33.000}, 22200, 22194, 22195, 22192, 22193, {186,187,180,173,172,178}},
  {180, { 20.785, -30.000}, 22201, 22198, 22199, 22196, 22197, {187,188,181,174,173,179}},
  {181, { 25.981, -27.000}, 22202, 22206, 22207, 22204, 22205, {188,189,182,175,174,180}},
  {182, { 31.177, -24.000}, 22203, 22210, 22211, 22208, 22209, {189,190,249,243,175,181}},
  {183, {  0.000, -48.000}, 22232, 22226, 22227, 22224, 22225, {191,192,184,176,117,126}},
  {184, {  5.196, -45.000}, 22233, 22230, 22231, 22228, 22229, {192,193,185,177,176,183}},
  {185, { 10.392, -42.000}, 22234, 22238, 22239, 22236, 22237, {193,194,186,178,177,184}},
  {186, { 15.588, -39.000}, 22235, 22242, 22243, 22240, 22241, {194,195,187,179,178,185}},
  {187, { 20.785, -36.000}, 22264, 22258, 22259, 22256, 22257, {195,196,188,180,179,186}},
  {188, { 25.981, -33.000}, 22265, 22262, 22263, 22260, 22261, {196,197,189,181,180,187}},
  {189, { 31.177, -30.000}, 22266, 22270, 22271, 22268, 22269, {197,198,190,182,181,188}},
  {190, { 36.373, -27.000}, 22267, 22274, 22275, 22272, 22273, {198,199,256,249,182,189}},
  {191, {  0.000, -54.000}, 22296, 22290, 22291, 22288, 22289, {200,201,192,183,126,136}},
  {192, {  5.196, -51.000}, 22297, 22294, 22295, 22292, 22293, {201,202,193,184,183,191}},
  {193, { 10.392, -48.000}, 22298, 22302, 22303, 22300, 22301, {202,203,194,185,184,192}},
  {194, { 15.588, -45.000}, 22299, 22306, 22307, 22304, 22305, {203,204,195,186,185,193}},
  {195, { 20.785, -42.000}, 22328, 22322, 22323, 22320, 22321, {204,205,196,187,186,194}},
  {196, { 25.981, -39.000}, 22329, 22326, 22327, 22324, 22325, {205,206,197,188,187,195}},
  {197, { 31.177, -36.000}, 22330, 22334, 22335, 22332, 22333, {206,207,198,189,188,196}},
  {198, { 36.373, -33.000}, 22331, 22338, 22339, 22336, 22337, {207,208,199,190,189,197}},
  {199, { 41.569, -30.000}, 22360, 22354, 22355, 22352, 22353, {208,209,264,256,190,198}},
  {200, {  0.000, -60.000}, 22361, 22358, 22359, 22356, 22357, {201,191,136}},
  {201, {  5.196, -57.000}, 22362, 22366, 22367, 22364, 22365, {210,202,192,191,200}},
  {202, { 10.392, -54.000}, 22363, 22370, 22371, 22368, 22369, {210,211,203,193,192,201}},
  {203, { 15.588, -51.000}, 22392, 22386, 22387, 22384, 22385, {211,212,204,194,193,202}},
  {204, { 20.785, -48.000}, 22393, 22390, 22391, 22388, 22389, {212,213,205,195,194,203}},
  {205, { 25.981, -45.000}, 22394, 22398, 22399, 22396, 22397, {213,214,206,196,195,204}},
  {206, { 31.177, -42.000}, 22395, 22402, 22403, 22400, 22401, {214,215,207,197,196,205}},
  {207, { 36.373, -39.000}, 22424, 22418, 22419, 22416, 22417, {215,216,208,198,197,206}},
  {208, { 41.569, -36.000}, 22425, 22422, 22423, 22420, 22421, {216,217,209,199,198,207}},
  {209, { 46.765, -33.000}, 22426, 22430, 22431, 22428, 22429, {217,218,273,264,199,208}},
  {210, { 10.392, -60.000}, 22427, 22434, 22435, 22432, 22433, {211,202,201}},
  {211, { 15.588, -57.000}, 22456, 22450, 22451, 22448, 22449, {212,203,202,210}},
  {212, { 20.785, -54.000}, 22457, 22454, 22455, 22452, 22453, {213,204,203,211}},
  {213, { 25.981, -51.000}, 22458, 22462, 22463, 22460, 22461, {214,205,204,212}},
  {214, { 31.177, -48.000}, 22459, 22466, 22467, 22464, 22465, {215,206,205,213}},
  {215, { 36.373, -45.000}, 22488, 22482, 22483, 22480, 22481, {216,207,206,214}},
  {216, { 41.569, -42.000}, 22489, 22486, 22487, 22484, 22485, {217,208,207,215}},
  {217, { 46.765, -39.000}, 22490, 22494, 22495, 22492, 22493, {218,209,208,216}},
  {218, { 51.962, -36.000}, 22491, 22498, 22499, 22496, 22497, {273,209,217}},
  {231, { 15.588,  -9.000}, 23011, 23018, 23019, 23016, 23017, {234,235,232,223,225,160,164}},
  {232, { 15.588,  -3.000}, 23040, 23034, 23035, 23032, 23033, {235,236,233,227,229,223,225,231}},
  {233, { 15.588,   3.000}, 23041, 23038, 23039, 23036, 23037, {236,237,304,296,298,227,229,232}},
  {234, { 20.785, -12.000}, 23042, 23046, 23047, 23044, 23045, {169,238,239,235,231,164}},
  {235, { 20.785,  -6.000}, 23043, 23050, 23051, 23048, 23049, {234,239,240,236,232,231}},
  {236, { 20.785,   0.000}, 23072, 23066, 23067, 23064, 23065, {235,240,241,237,233,232}},
  {237, { 20.785,   6.000}, 23073, 23070, 23071, 23068, 23069, {236,241,242,307,304,233}},
  {238, { 25.981, -15.000}, 23074, 23078, 23079, 23076, 23077, {175,243,244,239,234,169}},
  {239, { 25.981,  -9.000}, 23075, 23082, 23083, 23080, 23081, {238,244,245,240,235,234}},
  {240, { 25.981,  -3.000}, 23104, 23098, 23099, 23096, 23097, {239,245,246,241,236,235}},
  {241, { 25.981,   3.000}, 23105, 23102, 23103, 23100, 23101, {240,246,247,242,237,236}},
  {242, { 25.981,   9.000}, 23106, 23110, 23111, 23108, 23109, {241,247,248,311,307,237}},
  {243, { 31.177, -18.000}, 23107, 23114, 23115, 23112, 23113, {182,249,250,244,238,175}},
  {244, { 31.177, -12.000}, 23136, 23130, 23131, 23128, 23129, {243,250,251,245,239,238}},
  {245, { 31.177,  -6.000}, 23137, 23134, 23135, 23132, 23133, {244,251,252,246,240,239}},
  {246, { 31.177,   0.000}, 23138, 23142, 23143, 23140, 23141, {245,252,253,247,241,240}},
  {247, { 31.177,   6.000}, 23139, 23146, 23147, 23144, 23145, {246,253,254,248,242,241}},
  {248, { 31.177,  12.000}, 23168, 23162, 23163, 23160, 23161, {247,254,255,316,311,242}},
  {249, { 36.373, -21.000}, 23169, 23166, 23167, 23164, 23165, {190,256,257,250,243,182}},
  {250, { 36.373, -15.000}, 23170, 23174, 23175, 23172, 23173, {249,257,258,251,244,243}},
  {251, { 36.373,  -9.000}, 23171, 23178, 23179, 23176, 23177, {250,258,259,252,245,244}},
  {252, { 36.373,  -3.000}, 23200, 23194, 23195, 23192, 23193, {251,259,260,253,246,245}},
  {253, { 36.373,   3.000}, 23201, 23198, 23199, 23196, 23197, {252,260,261,254,247,246}},
  {254, { 36.373,   9.000}, 23202, 23206, 23207, 23204, 23205, {253,261,262,255,248,247}},
  {255, { 36.373,  15.000}, 23203, 23210, 23211, 23208, 23209, {254,262,263,322,316,248}},
  {256, { 41.569, -24.000}, 23232, 23226, 23227, 23224, 23225, {199,264,265,257,249,190}},
  {257, { 41.569, -18.000}, 23233, 23230, 23231, 23228, 23229, {256,265,266,258,250,249}},
  {258, { 41.569, -12.000}, 23234, 23238, 23239, 23236, 23237, {257,266,267,259,251,250}},
  {259, { 41.569,  -6.000}, 23235, 23242, 23243, 23240, 23241, {258,267,268,260,252,251}},
  {260, { 41.569,   0.000}, 23264, 23258, 23259, 23256, 23257, {259,268,269,261,253,252}},
  {261, { 41.569,   6.000}, 23265, 23262, 23263, 23260, 23261, {260,269,270,262,254,253}},
  {262, { 41.569,  12.000}, 23266, 23270, 23271, 23268, 23269, {261,270,271,263,255,254}},
  {263, { 41.569,  18.000}, 23267, 23274, 23275, 23272, 23273, {262,271,272,329,322,255}},
  {264, { 46.765, -27.000}, 23296, 23290, 23291, 23288, 23289, {209,273,274,265,256,199}},
  {265, { 46.765, -21.000}, 23297, 23294, 23295, 23292, 23293, {264,274,275,266,257,256}},
  {266, { 46.765, -15.000}, 23298, 23302, 23303, 23300, 23301, {265,275,276,267,258,257}},
  {267, { 46.765,  -9.000}, 23299, 23306, 23307, 23304, 23305, {266,276,277,268,259,258}},
  {268, { 46.765,  -3.000}, 23328, 23322, 23323, 23320, 23321, {267,277,278,269,260,259}},
  {269, { 46.765,   3.000}, 23329, 23326, 23327, 23324, 23325, {268,278,279,270,261,260}},
  {270, { 46.765,   9.000}, 23330, 23334, 23335, 23332, 23333, {269,279,280,271,262,261}},
  {271, { 46.765,  15.000}, 23331, 23338, 23339, 23336, 23337, {270,280,281,272,263,262}},
  {272, { 46.765,  21.000}, 23360, 23354, 23355, 23352, 23353, {271,281,282,337,329,263}},
  {273, { 51.962, -30.000}, 23361, 23358, 23359, 23356, 23357, {283,274,264,209,218}},
  {274, { 51.962, -24.000}, 23362, 23366, 23367, 23364, 23365, {273,283,284,275,265,264}},
  {275, { 51.962, -18.000}, 23363, 23370, 23371, 23368, 23369, {274,284,285,276,266,265}},
  {276, { 51.962, -12.000}, 23392, 23386, 23387, 23384, 23385, {275,285,286,277,267,266}},
  {277, { 51.962,  -6.000}, 23393, 23390, 23391, 23388, 23389, {276,286,287,278,268,267}},
  {278, { 51.962,   0.000}, 23394, 23398, 23399, 23396, 23397, {277,287,288,279,269,268}},
  {279, { 51.962,   6.000}, 23395, 23402, 23403, 23400, 23401, {278,288,289,280,270,269}},
  {280, { 51.962,  12.000}, 23424, 23418, 23419, 23416, 23417, {279,289,290,281,271,270}},
  {281, { 51.962,  18.000}, 23425, 23422, 23423, 23420, 23421, {280,290,291,282,272,271}},
  {282, { 51.962,  24.000}, 23426, 23430, 23431, 23428, 23429, {281,291,346,337,272}},
  {283, { 57.158, -27.000}, 23427, 23434, 23435, 23432, 23433, {273,274,284}},
  {284, { 57.158, -21.000}, 23456, 23450, 23451, 23448, 23449, {283,285,275,274}},
  {285, { 57.158, -15.000}, 23457, 23454, 23455, 23452, 23453, {284,286,276,275}},
  {286, { 57.158,  -9.000}, 23458, 23462, 23463, 23460, 23461, {285,287,277,276}},
  {287, { 57.158,  -3.000}, 23459, 23466, 23467, 23464, 23465, {286,288,278,277}},
  {288, { 57.158,   3.000}, 23488, 23482, 23483, 23480, 23481, {287,289,279,278}},
  {289, { 57.158,   9.000}, 23489, 23486, 23487, 23484, 23485, {288,290,280,279}},
  {290, { 57.158,  15.000}, 23490, 23494, 23495, 23492, 23493, {289,291,281,280}},
  {291, { 57.158,  21.000}, 23491, 23498, 23499, 23496, 23497, {290,282,281}},
  {304, { 15.588,   9.000}, 24011, 24018, 24019, 24016, 24017, {237,307,308,305,296,298,233}},
  {305, { 10.392,  12.000}, 24040, 24034, 24035, 24032, 24033, {304,308,309,306,300,302,297,296}},
  {306, {  5.196,  15.000}, 24041, 24038, 24039, 24036, 24037, {305,309,310,377,369,371,301,300}},
  {307, { 20.785,  12.000}, 24042, 24046, 24047, 24044, 24045, {237,242,311,312,308,304}},
  {308, { 15.588,  15.000}, 24043, 24050, 24051, 24048, 24049, {304,307,312,313,309,305}},
  {309, { 10.392,  18.000}, 24072, 24066, 24067, 24064, 24065, {305,308,313,314,310,306}},
  {310, {  5.196,  21.000}, 24073, 24070, 24071, 24068, 24069, {306,309,314,315,380,377}},
  {311, { 25.981,  15.000}, 24074, 24078, 24079, 24076, 24077, {242,248,316,317,312,307}},
  {312, { 20.785,  18.000}, 24075, 24082, 24083, 24080, 24081, {307,311,317,318,313,308}},
  {313, { 15.588,  21.000}, 24104, 24098, 24099, 24096, 24097, {308,312,318,319,314,309}},
  {314, { 10.392,  24.000}, 24105, 24102, 24103, 24100, 24101, {309,313,319,320,315,310}},
  {315, {  5.196,  27.000}, 24106, 24110, 24111, 24108, 24109, {310,314,320,321,384,380}},
  {316, { 31.177,  18.000}, 24107, 24114, 24115, 24112, 24113, {248,255,322,323,317,311}},
  {317, { 25.981,  21.000}, 24136, 24130, 24131, 24128, 24129, {311,316,323,324,318,312}},
  {318, { 20.785,  24.000}, 24137, 24134, 24135, 24132, 24133, {312,317,324,325,319,313}},
  {319, { 15.588,  27.000}, 24138, 24142, 24143, 24140, 24141, {313,318,325,326,320,314}},
  {320, { 10.392,  30.000}, 24139, 24146, 24147, 24144, 24145, {314,319,326,327,321,315}},
  {321, {  5.196,  33.000}, 24168, 24162, 24163, 24160, 24161, {315,320,327,328,389,384}},
  {322, { 36.373,  21.000}, 24169, 24166, 24167, 24164, 24165, {255,263,329,330,323,316}},
  {323, { 31.177,  24.000}, 24170, 24174, 24175, 24172, 24173, {316,322,330,331,324,317}},
  {324, { 25.981,  27.000}, 24171, 24178, 24179, 24176, 24177, {317,323,331,332,325,318}},
  {325, { 20.785,  30.000}, 24200, 24194, 24195, 24192, 24193, {318,324,332,333,326,319}},
  {326, { 15.588,  33.000}, 24201, 24198, 24199, 24196, 24197, {319,325,333,334,327,320}},
  {327, { 10.392,  36.000}, 24202, 24206, 24207, 24204, 24205, {320,326,334,335,328,321}},
  {328, {  5.196,  39.000}, 24203, 24210, 24211, 24208, 24209, {321,327,335,336,395,389}},
  {329, { 41.569,  24.000}, 24232, 24226, 24227, 24224, 24225, {263,272,337,338,330,322}},
  {330, { 36.373,  27.000}, 24233, 24230, 24231, 24228, 24229, {322,329,338,339,331,323}},
  {331, { 31.177,  30.000}, 24234, 24238, 24239, 24236, 24237, {323,330,339,340,332,324}},
  {332, { 25.981,  33.000}, 24235, 24242, 24243, 24240, 24241, {324,331,340,341,333,325}},
  {333, { 20.785,  36.000}, 24264, 24258, 24259, 24256, 24257, {325,332,341,342,334,326}},
  {334, { 15.588,  39.000}, 24265, 24262, 24263, 24260, 24261, {326,333,342,343,335,327}},
  {335, { 10.392,  42.000}, 24266, 24270, 24271, 24268, 24269, {327,334,343,344,336,328}},
  {336, {  5.196,  45.000}, 24267, 24274, 24275, 24272, 24273, {328,335,344,345,402,395}},
  {337, { 46.765,  27.000}, 24296, 24290, 24291, 24288, 24289, {272,282,346,347,338,329}},
  {338, { 41.569,  30.000}, 24297, 24294, 24295, 24292, 24293, {329,337,347,348,339,330}},
  {339, { 36.373,  33.000}, 24298, 24302, 24303, 24300, 24301, {330,338,348,349,340,331}},
  {340, { 31.177,  36.000}, 24299, 24306, 24307, 24304, 24305, {331,339,349,350,341,332}},
  {341, { 25.981,  39.000}, 24328, 24322, 24323, 24320, 24321, {332,340,350,351,342,333}},
  {342, { 20.785,  42.000}, 24329, 24326, 24327, 24324, 24325, {333,341,351,352,343,334}},
  {343, { 15.588,  45.000}, 24330, 24334, 24335, 24332, 24333, {334,342,352,353,344,335}},
  {344, { 10.392,  48.000}, 24331, 24338, 24339, 24336, 24337, {335,343,353,354,345,336}},
  {345, {  5.196,  51.000}, 24360, 24354, 24355, 24352, 24353, {336,344,354,355,410,402}},
  {346, { 51.962,  30.000}, 24361, 24358, 24359, 24356, 24357, {282,356,347,337}},
  {347, { 46.765,  33.000}, 24362, 24366, 24367, 24364, 24365, {337,346,356,357,348,338}},
  {348, { 41.569,  36.000}, 24363, 24370, 24371, 24368, 24369, {338,347,357,358,349,339}},
  {349, { 36.373,  39.000}, 24392, 24386, 24387, 24384, 24385, {339,348,358,359,350,340}},
  {350, { 31.177,  42.000}, 24393, 24390, 24391, 24388, 24389, {340,349,359,360,351,341}},
  {351, { 25.981,  45.000}, 24394, 24398, 24399, 24396, 24397, {341,350,360,361,352,342}},
  {352, { 20.785,  48.000}, 24395, 24402, 24403, 24400, 24401, {342,351,361,362,353,343}},
  {353, { 15.588,  51.000}, 24424, 24418, 24419, 24416, 24417, {343,352,362,363,354,344}},
  {354, { 10.392,  54.000}, 24425, 24422, 24423, 24420, 24421, {344,353,363,364,355,345}},
  {355, {  5.196,  57.000}, 24426, 24430, 24431, 24428, 24429, {345,354,364,419,410}},
  {356, { 51.962,  36.000}, 24427, 24434, 24435, 24432, 24433, {346,357,347}},
  {357, { 46.765,  39.000}, 24456, 24450, 24451, 24448, 24449, {356,347,348,358}},
  {358, { 41.569,  42.000}, 24457, 24454, 24455, 24452, 24453, {357,348,349,359}},
  {359, { 36.373,  45.000}, 24458, 24462, 24463, 24460, 24461, {358,349,350,360}},
  {360, { 31.177,  48.000}, 24459, 24466, 24467, 24464, 24465, {359,350,351,361}},
  {361, { 25.981,  51.000}, 24488, 24482, 24483, 24480, 24481, {360,351,352,362}},
  {362, { 20.785,  54.000}, 24489, 24486, 24487, 24484, 24485, {361,352,353,363}},
  {363, { 15.588,  57.000}, 24490, 24494, 24495, 24492, 24493, {362,353,354,364}},
  {364, { 10.392,  60.000}, 24491, 24498, 24499, 24496, 24497, {363,354,355}},
  {377, {  0.000,  18.000}, 25011, 25018, 25019, 25016, 25017, {306,310,380,381,378,370,369}},
  {378, { -5.196,  15.000}, 25040, 25034, 25035, 25032, 25033, {370,377,381,382,379,374,373,372}},
  {379, {-10.392,  12.000}, 25041, 25038, 25039, 25036, 25037, {374,378,382,383,12,5,4,376}},
  {380, {  0.000,  24.000}, 25042, 25046, 25047, 25044, 25045, {377,310,315,384,385,381}},
  {381, { -5.196,  21.000}, 25043, 25050, 25051, 25048, 25049, {378,377,380,385,386,382}},
  {382, {-10.392,  18.000}, 25072, 25066, 25067, 25064, 25065, {379,378,381,386,387,383}},
  {383, {-15.588,  15.000}, 25073, 25070, 25071, 25068, 25069, {12,379,382,387,388,15}},
  {384, {  0.000,  30.000}, 25074, 25078, 25079, 25076, 25077, {380,315,321,389,390,385}},
  {385, { -5.196,  27.000}, 25075, 25082, 25083, 25080, 25081, {381,380,384,390,391,386}},
  {386, {-10.392,  24.000}, 25104, 25098, 25099, 25096, 25097, {382,381,385,391,392,387}},
  {387, {-15.588,  21.000}, 25105, 25102, 25103, 25100, 25101, {383,382,386,392,393,388}},
  {388, {-20.785,  18.000}, 25106, 25110, 25111, 25108, 25109, {15,383,387,393,394,19}},
  {389, {  0.000,  36.000}, 25107, 25114, 25115, 25112, 25113, {384,321,328,395,396,390}},
  {390, { -5.196,  33.000}, 25136, 25130, 25131, 25128, 25129, {385,384,389,396,397,391}},
  {391, {-10.392,  30.000}, 25137, 25134, 25135, 25132, 25133, {386,385,390,397,398,392}},
  {392, {-15.588,  27.000}, 25138, 25142, 25143, 25140, 25141, {387,386,391,398,399,393}},
  {393, {-20.785,  24.000}, 25139, 25146, 25147, 25144, 25145, {388,387,392,399,400,394}},
  {394, {-25.981,  21.000}, 25168, 25162, 25163, 25160, 25161, {19,388,393,400,401,24}},
  {395, {  0.000,  42.000}, 25169, 25166, 25167, 25164, 25165, {389,328,336,402,403,396}},
  {396, { -5.196,  39.000}, 25170, 25174, 25175, 25172, 25173, {390,389,395,403,404,397}},
  {397, {-10.392,  36.000}, 25171, 25178, 25179, 25176, 25177, {391,390,396,404,405,398}},
  {398, {-15.588,  33.000}, 25200, 25194, 25195, 25192, 25193, {392,391,397,405,406,399}},
  {399, {-20.785,  30.000}, 25201, 25198, 25199, 25196, 25197, {393,392,398,406,407,400}},
  {400, {-25.981,  27.000}, 25202, 25206, 25207, 25204, 25205, {394,393,399,407,408,401}},
  {401, {-31.177,  24.000}, 25203, 25210, 25211, 25208, 25209, {24,394,400,408,409,30}},
  {402, {  0.000,  48.000}, 25232, 25226, 25227, 25224, 25225, {395,336,345,410,411,403}},
  {403, { -5.196,  45.000}, 25233, 25230, 25231, 25228, 25229, {396,395,402,411,412,404}},
  {404, {-10.392,  42.000}, 25234, 25238, 25239, 25236, 25237, {397,396,403,412,413,405}},
  {405, {-15.588,  39.000}, 25235, 25242, 25243, 25240, 25241, {398,397,404,413,414,406}},
  {406, {-20.785,  36.000}, 25264, 25258, 25259, 25256, 25257, {399,398,405,414,415,407}},
  {407, {-25.981,  33.000}, 25265, 25262, 25263, 25260, 25261, {400,399,406,415,416,408}},
  {408, {-31.177,  30.000}, 25266, 25270, 25271, 25268, 25269, {401,400,407,416,417,409}},
  {409, {-36.373,  27.000}, 25267, 25274, 25275, 25272, 25273, {30,401,408,417,418,37}},
  {410, {  0.000,  54.000}, 25296, 25290, 25291, 25288, 25289, {402,345,355,419,420,411}},
  {411, { -5.196,  51.000}, 25297, 25294, 25295, 25292, 25293, {403,402,410,420,421,412}},
  {412, {-10.392,  48.000}, 25298, 25302, 25303, 25300, 25301, {404,403,411,421,422,413}},
  {413, {-15.588,  45.000}, 25299, 25306, 25307, 25304, 25305, {405,404,412,422,423,414}},
  {414, {-20.785,  42.000}, 25328, 25322, 25323, 25320, 25321, {406,405,413,423,424,415}},
  {415, {-25.981,  39.000}, 25329, 25326, 25327, 25324, 25325, {407,406,414,424,425,416}},
  {416, {-31.177,  36.000}, 25330, 25334, 25335, 25332, 25333, {408,407,415,425,426,417}},
  {417, {-36.373,  33.000}, 25331, 25338, 25339, 25336, 25337, {409,408,416,426,427,418}},
  {418, {-41.569,  30.000}, 25360, 25354, 25355, 25352, 25353, {37,409,417,427,428,45}},
  {419, {  0.000,  60.000}, 25361, 25358, 25359, 25356, 25357, {355,410,420}},
  {420, { -5.196,  57.000}, 25362, 25366, 25367, 25364, 25365, {419,410,411,421,429}},
  {421, {-10.392,  54.000}, 25363, 25370, 25371, 25368, 25369, {412,411,420,429,430,422}},
  {422, {-15.588,  51.000}, 25392, 25386, 25387, 25384, 25385, {413,412,421,430,431,423}},
  {423, {-20.785,  48.000}, 25393, 25390, 25391, 25388, 25389, {414,413,422,431,432,424}},
  {424, {-25.981,  45.000}, 25394, 25398, 25399, 25396, 25397, {415,414,423,432,433,425}},
  {425, {-31.177,  42.000}, 25395, 25402, 25403, 25400, 25401, {416,415,424,433,434,426}},
  {426, {-36.373,  39.000}, 25424, 25418, 25419, 25416, 25417, {417,416,425,434,435,427}},
  {427, {-41.569,  36.000}, 25425, 25422, 25423, 25420, 25421, {418,417,426,435,436,428}},
  {428, {-46.765,  33.000}, 25426, 25430, 25431, 25428, 25429, {45,418,427,436,437,54}},
  {429, {-10.392,  60.000}, 25427, 25434, 25435, 25432, 25433, {420,421,430}},
  {430, {-15.588,  57.000}, 25456, 25450, 25451, 25448, 25449, {429,421,422,431}},
  {431, {-20.785,  54.000}, 25457, 25454, 25455, 25452, 25453, {430,422,423,432}},
  {432, {-25.981,  51.000}, 25458, 25462, 25463, 25460, 25461, {431,423,424,433}},
  {433, {-31.177,  48.000}, 25459, 25466, 25467, 25464, 25465, {432,424,425,434}},
  {434, {-36.373,  45.000}, 25488, 25482, 25483, 25480, 25481, {433,425,426,435}},
  {435, {-41.569,  42.000}, 25489, 25486, 25487, 25484, 25485, {434,426,427,436}},
  {436, {-46.765,  39.000}, 25490, 25494, 25495, 25492, 25493, {435,427,428,437}},
  {437, {-51.962,  36.000}, 25491, 25498, 25499, 25496, 25497, {436,428,54}},
};