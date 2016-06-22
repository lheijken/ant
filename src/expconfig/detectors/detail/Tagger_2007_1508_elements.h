#include <detectors/Tagger.h>
const std::vector<ant::expconfig::detector::Tagger::Element_t>
ant::expconfig::detector::Tagger_2007::elements_init =
{
    // channel, TDC, Scaler, ADC, electronEnergy
    {  0, 1438,  164,  816, 105.500},
    {  1, 1432,  152,  800, 107.500},
    {  2, 1454,  165,  817, 109.500},
    {  3, 1448,  153,  801, 111.600},
    {  4, 1439,  166,  818, 113.800},
    {  5, 1433,  154,  802, 116.000},
    {  6, 1455,  167,  819, 118.400},
    {  7, 1449,  155,  803, 120.800},
    {  8, 1440,  168,  820, 123.400},
    {  9, 1434,  156,  804, 126.000},
    { 10, 1456,  169,  821, 128.700},
    { 11, 1450,  157,  805, 131.400},
    { 12, 1441,  170,  822, 134.300},
    { 13, 1435,  158,  806, 137.200},
    { 14, 1457,  171,  823, 140.100},
    { 15, 1451,  159,  807, 143.200},
    { 16, 1442,  172,  824, 146.200},
    { 17, 1436,  160,  808, 149.300},
    { 18, 1458,  173,  825, 152.500},
    { 19, 1452,  161,  809, 155.600},
    { 20, 1443,  174,  826, 158.800},
    { 21, 1437,  162,  810, 161.900},
    { 22, 1459,  175,  827, 165.100},
    { 23, 1453,  163,  811, 168.200},
    { 24, 1444,  176,  828, 171.300},
    { 25, 1460,  177,  829, 174.400},
    { 26, 1445,  178,  830, 177.400},
    { 27, 1461,  179,  831, 180.300},
    { 28, 1446,  180,  812, 183.300},
    { 29, 1462,  181,  813, 186.100},
    { 30, 1447,  182,  814, 188.900},
    { 31, 1463,  183,  815, 191.700},
    { 32, 1464,  184,  832, 194.400},
    { 33, 1480,  185,  833, 197.100},
    { 34, 1465,  186,  834, 199.800},
    { 35, 1481,  187,  835, 202.500},
    { 36, 1466,  188,  836, 205.300},
    { 37, 1482,  189,  837, 208.000},
    { 38, 1467,  190,  838, 210.800},
    { 39, 1483,  191,  839, 213.600},
    { 40, 1468,  192,  840, 216.400},
    { 41, 1484,  193,  841, 219.200},
    { 42, 1469,  194,  842, 222.000},
    { 43, 1485,  195,  843, 224.800},
    { 44, 1470,  196,  844, 227.700},
    { 45, 1486,  197,  845, 230.500},
    { 46, 1471,  198,  846, 233.400},
    { 47, 1487,  199,  847, 236.200},
    { 48, 1472,  200,  848, 239.100},
    { 49, 1488,  201,  849, 242.000},
    { 50, 1473,  202,  850, 244.900},
    { 51, 1489,  203,  851, 247.800},
    { 52, 1474,  204,  852, 250.700},
    { 53, 1490,  205,  853, 253.700},
    { 54, 1475,  206,  854, 256.600},
    { 55, 1491,  207,  855, 259.600},
    { 56, 1476,  208,  856, 262.600},
    { 57, 1492,  209,  857, 265.500},
    { 58, 1477,  210,  858, 268.500},
    { 59, 1493,  211,  859, 271.500},
    { 60, 1478,  212,  860, 274.500},
    { 61, 1494,  213,  861, 277.500},
    { 62, 1479,  214,  862, 280.600},
    { 63, 1495,  215,  863, 283.600},
    { 64, 1496,  216,  864, 286.700},
    { 65, 1512,  217,  865, 289.800},
    { 66, 1497,  218,  866, 292.800},
    { 67, 1513,  219,  867, 296.000},
    { 68, 1498,  220,  868, 299.000},
    { 69, 1514,  221,  869, 302.200},
    { 70, 1499,  222,  870, 305.300},
    { 71, 1515,  223,  871, 308.400},
    { 72, 1500,  224,  872, 311.600},
    { 73, 1516,  225,  873, 314.700},
    { 74, 1501,  226,  874, 317.900},
    { 75, 1517,  227,  875, 321.100},
    { 76, 1502,  228,  876, 324.300},
    { 77, 1518,  229,  877, 327.500},
    { 78, 1503,  230,  878, 330.700},
    { 79, 1519,  231,  879, 334.000},
    { 80, 1504,  232,  880, 337.200},
    { 81, 1520,  233,  881, 340.500},
    { 82, 1505,  234,  882, 343.700},
    { 83, 1521,  235,  883, 347.000},
    { 84, 1506,  236,  884, 350.300},
    { 85, 1522,  237,  885, 353.600},
    { 86, 1507,  238,  886, 356.900},
    { 87, 1523,  239,  887, 360.200},
    { 88, 1508,  240,  888, 363.600},
    { 89, 1524,  241,  889, 366.900},
    { 90, 1509,  242,  890, 370.200},
    { 91, 1525,  243,  891, 373.600},
    { 92, 1510,  244,  892, 377.000},
    { 93, 1526,  245,  893, 380.400},
    { 94, 1511,  246,  894, 383.800},
    { 95, 1527,  247,  895, 387.200},
    { 96, 1528,  248,  896, 390.600},
    { 97, 1544,  249,  897, 394.000},
    { 98, 1529,  250,  898, 397.500},
    { 99, 1545,  251,  899, 400.900},
    {100, 1530,  252,  900, 404.400},
    {101, 1546,  253,  901, 407.900},
    {102, 1531,  254,  902, 411.300},
    {103, 1547,  255,  903, 414.900},
    {104, 1532,  256,  904, 418.300},
    {105, 1548,  257,  905, 421.900},
    {106, 1533,  258,  906, 425.400},
    {107, 1549,  259,  907, 428.900},
    {108, 1534,  260,  908, 432.500},
    {109, 1550,  261,  909, 436.000},
    {110, 1535,  262,  910, 439.600},
    {111, 1551,  263,  911, 443.200},
    {112, 1536,  264,  912, 446.700},
    {113, 1552,  265,  913, 450.400},
    {114, 1537,  266,  914, 453.900},
    {115, 1553,  267,  915, 457.600},
    {116, 1538,  268,  916, 461.200},
    {117, 1554,  269,  917, 464.800},
    {118, 1539,  270,  918, 468.500},
    {119, 1555,  271,  919, 472.100},
    {120, 1540,  272,  920, 475.800},
    {121, 1556,  273,  921, 479.500},
    {122, 1541,  274,  922, 483.100},
    {123, 1557,  275,  923, 486.900},
    {124, 1542,  276,  924, 490.500},
    {125, 1558,  277,  925, 494.200},
    {126, 1543,  278,  926, 498.000},
    {127, 1559,  279,  927, 501.700},
    {128, 1560,  280,  928, 505.400},
    {129, 1576,  281,  929, 509.200},
    {130, 1561,  282,  930, 512.900},
    {131, 1577,  283,  931, 516.700},
    {132, 1562,  284,  932, 520.400},
    {133, 1578,  285,  933, 524.300},
    {134, 1563,  286,  934, 528.000},
    {135, 1579,  287,  935, 531.800},
    {136, 1564,  288,  936, 535.600},
    {137, 1580,  289,  937, 539.400},
    {138, 1565,  290,  938, 543.300},
    {139, 1581,  291,  939, 547.100},
    {140, 1566,  292,  940, 550.900},
    {141, 1582,  293,  941, 554.800},
    {142, 1567,  294,  942, 558.600},
    {143, 1583,  295,  943, 562.500},
    {144, 1568,  296,  944, 566.400},
    {145, 1584,  297,  945, 570.200},
    {146, 1569,  298,  946, 574.100},
    {147, 1585,  299,  947, 578.000},
    {148, 1570,  300,  948, 581.900},
    {149, 1586,  301,  949, 585.800},
    {150, 1571,  302,  950, 589.700},
    {151, 1587,  303,  951, 593.700},
    {152, 1572,  304,  952, 597.600},
    {153, 1588,  305,  953, 601.500},
    {154, 1573,  306,  954, 605.500},
    {155, 1589,  307,  955, 609.500},
    {156, 1574,  308,  956, 613.400},
    {157, 1590,  309,  957, 617.400},
    {158, 1575,  310,  958, 621.300},
    {159, 1591,  311,  959, 625.300},
    {160, 1592,  312,  960, 629.300},
    {161, 1608,  313,  961, 633.300},
    {162, 1593,  314,  962, 637.300},
    {163, 1609,  315,  963, 641.300},
    {164, 1594,  316,  964, 645.300},
    {165, 1610,  317,  965, 649.300},
    {166, 1595,  318,  966, 653.400},
    {167, 1611,  319,  967, 657.400},
    {168, 1596,  320,  968, 661.400},
    {169, 1612,  321,  969, 665.500},
    {170, 1597,  322,  970, 669.500},
    {171, 1613,  323,  971, 673.600},
    {172, 1598,  324,  972, 677.600},
    {173, 1614,  325,  973, 681.700},
    {174, 1599,  326,  974, 685.800},
    {175, 1615,  327,  975, 689.900},
    {176, 1600,  328,  976, 693.900},
    {177, 1616,  329,  977, 698.100},
    {178, 1601,  330,  978, 702.100},
    {179, 1617,  331,  979, 706.200},
    {180, 1602,  332,  980, 710.300},
    {181, 1618,  333,  981, 714.400},
    {182, 1603,  334,  982, 718.600},
    {183, 1619,  335,  983, 722.700},
    {184, 1604,  336,  984, 726.800},
    {185, 1620,  337,  985, 730.900},
    {186, 1605,  338,  986, 735.100},
    {187, 1621,  339,  987, 739.200},
    {188, 1606,  340,  988, 743.400},
    {189, 1622,  341,  989, 747.500},
    {190, 1607,  342,  990, 751.700},
    {191, 1623,  343,  991, 755.800},
    {192, 1624,  344,  992, 760.000},
    {193, 1640,  345,  993, 764.100},
    {194, 1625,  346,  994, 768.300},
    {195, 1641,  347,  995, 772.500},
    {196, 1626,  348,  996, 776.700},
    {197, 1642,  349,  997, 780.800},
    {198, 1627,  350,  998, 785.000},
    {199, 1643,  351,  999, 789.200},
    {200, 1628,  352, 1000, 793.500},
    {201, 1644,  353, 1001, 797.600},
    {202, 1629,  354, 1002, 801.800},
    {203, 1645,  355, 1003, 806.000},
    {204, 1630,  356, 1004, 810.200},
    {205, 1646,  357, 1005, 814.400},
    {206, 1631,  358, 1006, 818.700},
    {207, 1647,  359, 1007, 822.900},
    {208, 1632,  360, 1008, 827.100},
    {209, 1648,  361, 1009, 831.300},
    {210, 1633,  362, 1010, 835.500},
    {211, 1649,  363, 1011, 839.800},
    {212, 1634,  364, 1012, 844.100},
    {213, 1650,  365, 1013, 848.200},
    {214, 1635,  366, 1014, 852.500},
    {215, 1651,  367, 1015, 856.700},
    {216, 1636,  368, 1016, 861.000},
    {217, 1652,  369, 1017, 865.200},
    {218, 1637,  370, 1018, 869.500},
    {219, 1653,  371, 1019, 873.700},
    {220, 1638,  372, 1020, 878.000},
    {221, 1654,  373, 1021, 882.200},
    {222, 1639,  374, 1022, 886.500},
    {223, 1655,  375, 1023, 890.700},
    {224, 1656,  376, 1024, 895.100},
    {225, 1672,  377, 1025, 899.300},
    {226, 1657,  378, 1026, 903.500},
    {227, 1673,  379, 1027, 907.800},
    {228, 1658,  380, 1028, 912.100},
    {229, 1674,  381, 1029, 916.300},
    {230, 1659,  382, 1030, 920.600},
    {231, 1675,  383, 1031, 924.900},
    {232, 1660,  384, 1032, 929.200},
    {233, 1676,  385, 1033, 933.400},
    {234, 1661,  386, 1034, 937.700},
    {235, 1677,  387, 1035, 942.000},
    {236, 1662,  388, 1036, 946.400},
    {237, 1678,  389, 1037, 950.600},
    {238, 1663,  390, 1038, 954.900},
    {239, 1679,  391, 1039, 959.100},
    {240, 1664,  392, 1040, 963.400},
    {241, 1680,  393, 1041, 967.700},
    {242, 1665,  394, 1042, 972.000},
    {243, 1681,  395, 1043, 976.300},
    {244, 1666,  396, 1044, 980.600},
    {245, 1682,  397, 1045, 984.900},
    {246, 1667,  398, 1046, 989.200},
    {247, 1683,  399, 1047, 993.500},
    {248, 1668,  400, 1048, 997.800},
    {249, 1684,  401, 1049, 1002.100},
    {250, 1669,  402, 1050, 1006.300},
    {251, 1685,  403, 1051, 1010.600},
    {252, 1670,  404, 1052, 1014.900},
    {253, 1686,  405, 1053, 1019.200},
    {254, 1671,  406, 1054, 1023.500},
    {255, 1687,  407, 1055, 1027.800},
    {256, 1688,  408, 1056, 1032.100},
    {257, 1704,  409, 1057, 1036.400},
    {258, 1689,  410, 1058, 1040.700},
    {259, 1705,  411, 1059, 1045.000},
    {260, 1690,  412, 1060, 1049.400},
    {261, 1706,  413, 1061, 1053.600},
    {262, 1691,  414, 1062, 1057.900},
    {263, 1707,  415, 1063, 1062.200},
    {264, 1692,  416, 1064, 1066.500},
    {265, 1708,  417, 1065, 1070.800},
    {266, 1693,  418, 1066, 1075.100},
    {267, 1709,  419, 1067, 1079.400},
    {268, 1694,  420, 1068, 1083.700},
    {269, 1710,  421, 1069, 1087.900},
    {270, 1695,  422, 1070, 1092.200},
    {271, 1711,  423, 1071, 1096.500},
    {272, 1696,  424, 1072, 1100.900},
    {273, 1712,  425, 1073, 1105.100},
    {274, 1697,  426, 1074, 1109.400},
    {275, 1713,  427, 1075, 1113.700},
    {276, 1698,  428, 1076, 1118.000},
    {277, 1714,  429, 1077, 1122.300},
    {278, 1699,  430, 1078, 1126.500},
    {279, 1715,  431, 1079, 1130.800},
    {280, 1700,  432, 1080, 1135.100},
    {281, 1716,  433, 1081, 1139.400},
    {282, 1701,  434, 1082, 1143.700},
    {283, 1717,  435, 1083, 1148.000},
    {284, 1702,  436, 1084, 1152.300},
    {285, 1718,  437, 1085, 1156.500},
    {286, 1703,  438, 1086, 1160.800},
    {287, 1719,  439, 1087, 1165.100},
    {288, 1720,  440, 1088, 1169.300},
    {289, 1736,  441, 1089, 1173.600},
    {290, 1721,  442, 1090, 1177.900},
    {291, 1737,  443, 1091, 1182.100},
    {292, 1722,  444, 1092, 1186.400},
    {293, 1738,  445, 1093, 1190.700},
    {294, 1723,  446, 1094, 1194.900},
    {295, 1739,  447, 1095, 1199.200},
    {296, 1724,  448, 1096, 1203.500},
    {297, 1740,  449, 1097, 1207.700},
    {298, 1725,  450, 1098, 1212.000},
    {299, 1741,  451, 1099, 1216.200},
    {300, 1726,  452, 1100, 1220.500},
    {301, 1742,  453, 1101, 1224.700},
    {302, 1727,  454, 1102, 1229.000},
    {303, 1743,  455, 1103, 1233.200},
    {304, 1728,  456, 1104, 1237.500},
    {305, 1744,  457, 1105, 1241.700},
    {306, 1729,  458, 1106, 1245.900},
    {307, 1745,  459, 1107, 1250.200},
    {308, 1730,  460, 1108, 1254.500},
    {309, 1746,  461, 1109, 1258.700},
    {310, 1731,  462, 1110, 1262.900},
    {311, 1747,  463, 1111, 1267.100},
    {312, 1732,  464, 1112, 1271.300},
    {313, 1748,  465, 1113, 1275.600},
    {314, 1733,  466, 1114, 1279.800},
    {315, 1749,  467, 1115, 1284.000},
    {316, 1734,  468, 1116, 1288.200},
    {317, 1750,  469, 1117, 1292.400},
    {318, 1735,  470, 1118, 1296.600},
    {319, 1751,  471, 1119, 1300.900},
    {320, 1752,  472, 1120, 1305.100},
    {321, 1768,  473, 1121, 1309.300},
    {322, 1753,  474, 1122, 1313.500},
    {323, 1769,  475, 1123, 1317.700},
    {324, 1754,  476, 1124, 1321.900},
    {325, 1770,  477, 1125, 1326.000},
    {326, 1755,  478, 1126, 1330.200},
    {327, 1771,  479, 1127, 1334.400},
    {328, 1756,  480, 1128, 1338.600},
    {329, 1772,  481, 1129, 1342.800},
    {330, 1757,  482, 1130, 1347.000},
    {331, 1773,  483, 1131, 1351.100},
    {332, 1758,  484, 1132, 1355.300},
    {333, 1774,  485, 1133, 1359.500},
    {334, 1759,  486, 1134, 1363.700},
    {335, 1775,  487, 1135, 1367.800},
    {336, 1760,  488, 1136, 1372.000},
    {337, 1776,  489, 1137, 1376.100},
    {338, 1761,  490, 1138, 1380.300},
    {339, 1777,  491, 1139, 1384.400},
    {340, 1762,  492, 1140, 1388.600},
    {341, 1778,  493, 1141, 1392.700},
    {342, 1763,  494, 1142, 1396.900},
    {343, 1779,  495, 1143, 1401.000},
    {344, 1764,  496, 1144, 1405.100},
    {345, 1780,  497, 1145, 1409.300},
    {346, 1765,  498, 1146, 1413.400},
    {347, 1781,  499, 1147, 1417.600},
    {348, 1766,  500, 1148, 1421.600},
    {349, 1782,  501, 1149, 1425.800},
    {350, 1767,  502, 1150, 1429.900},
    {351, 1783,  503, 1151, 1434.000}
};