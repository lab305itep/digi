	const struct {
		char name[32];
		int first;
		int last;
		int period;
		double bgnd;
	} positions[] = {
		{ "down_21.04.16",   2307, 2361, 1, 0.05},
		{ "up_22.04.16",     2366, 2387, 1, 0.05},
		{ "down_23.04.16",   2399, 2445, 1, 0.05},
		{ "mid_24.04.16",    2449, 2512, 1, 0.05},
		{ "up_25.04.16",     2514, 2562, 1, 0.05},
		{ "down_26.04.16",   2564, 2609, 1, 0.05},
		{ "mid_27.04.16",    2620, 2685, 1, 0.05},
		{ "up_28.04.16",     2687, 2730, 1, 0.05},
		{ "down_29.04.16",   2733, 2788, 1, 0.05},
		{ "mid_30.04.16",    2791, 2832, 1, 0.05},
		{ "stuck_01.05.16",  2836, 2952, 1, 0.05},
		{ "stuck_03.05.16",  2954, 3082, 1, 0.05},
		{ "stuck_05.05.16",  3086, 3197, 1, 0.05},
		{ "stuck_07.05.16",  3199, 3317, 1, 0.05},
		{ "stuck_09.05.16",  3319, 3442, 1, 0.05},
		{ "stuck_11.05.16",  3444, 3549, 1, 0.05},
		{ "stuck_13.05.16",  3560, 3682, 1, 0.05},
		{ "stuck_15.05.16",  3683, 3778, 1, 0.05},
		{ "up_03.06.16",     4261, 4400, 1, 0.1},
		{ "down_06.06.16",   4403, 4439, 1, 0.1},
//		{ "up_08.06.16",     4508, 4601, 1, 0.1},	// empty
		{ "raised_30.09.16", 5540, 5807, 0, 0.05},
		{ "raised_04.10.16", 5808, 5903, 0, 0.025},	// veto corners on
		{ "mid_05.10.16",    5907, 5995, 2, 0.025},
		{ "up_10.10.16",     6007, 6130, 2, 0.025},
		{ "mid_12.10.16",    6136, 6275, 0, 0.025},
		{ "down_14.10.16",   6278, 6467, 0, 0.025},
		{ "up_17.10.16",     6469, 6570, 0, 0.025},
		{ "mid_21.10.16",    6573, 6582, 0, 0.025},
		{ "down_21.10.16",   6587, 6745, 0, 0.025},
		{ "up_24.10.16",     6757, 6815, 0, 0.025},
		{ "mid_27.10.16",    6842, 6909, 0, 0.025},	// was 6923
		{ "down_28.10.16",   6926, 7095, 0, 0.025},
		{ "up_31.10.16",     7106, 7364, 2, 0.025},
		{ "mid_11.11.16",    7387, 7406, 2, 0.025},
		{ "down_11.11.16",   7418, 7458, 2, 0.025},
		{ "up_14.11.16",     7478, 7579, 2, 0.025},
		{ "mid_16.11.16",    7581, 7717, 0, 0.025},
		{ "down_18.11.16",   7727, 7913, 2, 0.025},
		{ "up_21.11.16",     7922, 8042, 2, 0.025},
		{ "mid_23.11.16",    8048, 8179, 2, 0.025},
		{ "down_25.11.16",   8185, 8353, 2, 0.025},
		{ "up_28.11.16",     8357, 8430, 2, 0.025},
		{ "mid_01.12.16",    8470, 8571, 2, 0.025},
		{ "down_02.12.16",   8574, 8738, 2, 0.025},
		{ "up_05.12.16",     8741, 8869, 2, 0.025},
		{ "mid_07.12.16",    8873, 9009, 2, 0.025},
		{ "up_12.12.16",     9012, 9112, 2, 0.025},
		{ "mid_14.12.16",    9116, 9245, 2, 0.025},
		{ "down_16.12.16",   9253, 9470, 2, 0.025},
		{ "up_19.12.16",     9475, 9600, 0, 0.025},
		{ "mid_21.12.16",    9603, 9712, 0, 0.025},
		{ "down_23.12.16",   9715, 9869, 0, 0.025},
		{ "up_26.12.16",     9871, 10019, 0, 0.025},
		{ "mid_28.12.16",    10021, 10171, 0, 0.025},
		{ "down_30.12.16",   10175, 10307, 0, 0.025},
		{ "down_02.01.17",   10308, 10356, 0, 0.025},
		{ "down_03.01.17",   10357, 10424, 0, 0.025},
		{ "up_04.01.17",     10433, 10832, 2, 0.025},
		{ "mid_11.01.17",    10834, 10973, 2, 0.025},
		{ "down_13.01.17",   10979, 11147, 2, 0.025},
		{ "up_16.01.17",     11150, 11267, 2, 0.025},
		{ "mid_18.01.17",    11271, 11401, 2, 0.025},
		{ "down_20.01.17",   11404, 11563, 2, 0.025},
		{ "up_23.01.17",     11570, 11694, 2, 0.025},
		{ "mid_25.01.17",    11695, 11816, 2, 0.025},
		{ "down_27.01.17",   11820, 12079, 2, 0.025},
		{ "down_01.02.17",   12080, 12421, 2, 0.025},	// there are sources inside this period
		{ "up_06.02.17",     12425, 12547, 2, 0.025},
		{ "mid_08.02.17",    12548, 12671, 2, 0.025},
		{ "down_10.02.17",   12672, 12936, 2, 0.025},	// undefined position inside
		{ "mid_14.02.17",    12938, 12984, 2, 0.025},
		{ "down_15.02.17",   12985, 12994, 2, 0.033},
		{ "up_16.02.17",     13050, 13115, 2, 0.033},
		{ "down_17.02.17",   13117, 13367, 2, 0.033},
		{ "up_20.02.17",     13368, 13504, 2, 0.033},
		{ "mid_22.02.17",    13505, 13663, 2, 0.033},
		{ "down_24.02.17",   13665, 13852, 2, 0.033},
		{ "up_27.02.17",     13854, 13901, 2, 0.033},
		{ "mid_28.02.17",    13902, 14145, 2, 0.033},
		{ "mid_03.03.17",    14147, 14359, 3, 0.033},
		{ "down_06.03.17",   14360, 14426, 3, 0.033},
		{ "down_10.03.17",   14636, 14831, 3, 0.033},
		{ "up_13.03.17",     14833, 14968, 3, 0.033},
		{ "mid_15.03.17",    14971, 15130, 3, 0.025},
		{ "down_17.03.17",   15132, 15306, 3, 0.025},
		{ "up_20.03.17",     15308, 15412, 3, 0.025},
		{ "mid_22.03.17",    15413, 15436, 3, 0.025},
		{ "mid_29.03.17",    17314, 17436, 3, 0.025},
		{ "down_31.03.17",   17438, 17612, 3, 0.025},
		{ "up_03.04.17",     17613, 17715, 3, 0.025},
		{ "mid_05.04.17",    17717, 17848, 3, 0.025},
		{ "down_07.04.17",   17849, 18020, 3, 0.025},
		{ "up_10.04.17",     18021, 18131, 3, 0.025},
		{ "mid_12.04.17",    18133, 18263, 3, 0.025},
		{ "down_14.04.17",   18264, 18444, 3, 0.025},
		{ "up_17.04.17",     18445, 18552, 3, 0.025},
		{ "mid_19.04.17",    18553, 18611, 3, 0.025},
		{ "mid_20.04.17",    18612, 18682, 0, 0.025},	// reactor off and getting back
		{ "down_21.04.17",   18684, 18858, 3, 0.025},
		{ "up_24.04.17",     18859, 18921, 3, 0.025},
		{ "mid_26.04.17",    18922, 19040, 3, 0.025},
		{ "down_28.04.17",   19041, 19221, 3, 0.025},
		{ "up_01.05.17",     19222, 19336, 3, 0.025},
		{ "down_03.05.17",   19337, 19339, 3, 0.025},
		{ "up_03.05.17",     19340, 19448, 3, 0.025},
		{ "mid_05.05.17",    19449, 19456, 3, 0.025},
//		{ "down_05.05.17",   19457, 19638, 3, 0.025},	// Empty
		{ "up_08.05.17",     19639, 19738, 3, 0.025},
		{ "mid_10.05.17",    19739, 19856, 3, 0.025},
		{ "down_12.05.17",   19857, 20041, 3, 0.025},
		{ "up_15.05.17",     20042, 20151, 3, 0.025},
		{ "down_17.05.17",   20195, 20474, 3, 0.025},	// 22Na and 60Co inside
		{ "up_22.05.17",     20475, 20606, 3, 0.025},
		{ "mid_24.05.17",    20612, 20752, 3, 0.025},
		{ "down_26.05.17",   20754, 20911, 3, 0.025},
		{ "up_29.05.17",     20915, 21027, 3, 0.025},
		{ "mid_31.05.17",    21028, 21147, 3, 0.025},	// position is identified as unknown by Ira.
		{ "down_02.06.17",   21148, 21313, 3, 0.025},
		{ "up_05.06.17",     21315, 21411, 3, 0.025},
		{ "mid_07.06.17",    21412, 21566, 3, 0.025},
		{ "up_13.06.17",     21616, 21672, 3, 0.025},
		{ "mid_14.06.17",    21673, 21800, 3, 0.025},
		{ "down_16.06.17",   21801, 21959, 3, 0.025},
		{ "up_19.06.17",     21960, 22077, 3, 0.025},
		{ "mid_21.06.17",    22079, 22200, 3, 0.025},
		{ "down_23.06.17",   22201, 22380, 3, 0.025},
		{ "up_26.06.17",     22381, 22497, 3, 0.025},
		{ "mid_28.06.17",    22499, 22637, 3, 0.025},
		{ "down_30.06.17",   22638, 22812, 3, 0.025},
		{ "up_03.07.17",     22813, 22924, 3, 0.025},
		{ "mid_05.07.17",    22925, 23041, 3, 0.025},
		{ "down_07.07.17",   23042, 23081, 0, 0.025},	// starting power off
		{ "down_08.07.17",   23082, 23185, 0, 0.025},	// OFF
		{ "up_10.07.17",     23188, 23281, 0, 0.025},
		{ "mid_12.07.17",    23283, 23380, 0, 0.025},
		{ "down_14.07.17",   23382, 23523, 0, 0.025},
		{ "up_17.07.17",     23525, 23624, 0, 0.025},
		{ "mid_19.07.17",    23625, 23724, 0, 0.025},
		{ "down_21.07.17",   23725, 23861, 0, 0.025},
		{ "up_24.07.17",     23863, 23957, 0, 0.025},
		{ "mid_26.07.17",    23958, 24060, 0, 0.025},
		{ "down_28.07.17",   24061, 24203, 0, 0.025},
		{ "up_31.07.17",     24204, 24209, 0, 0.025},
		{ "mid_02.08.17",    24212, 24343, 0, 0.025},
		{ "down_04.08.17",   24347, 24480, 0, 0.025},
		{ "up_07.08.17",     24482, 24575, 0, 0.025},
		{ "mid_09.08.17",    24576, 24670, 0, 0.025},
		{ "down_11.08.17",   24672, 24706, 0, 0.025},
		{ "down_18.08.17",   24727, 24782, 0, 0.025},
		{ "down_19.08.17",   24783, 24892, 4, 0.025},
		{ "up_21.08.17",     24893, 24908, 4, 0.025},
		{ "up_22.08.17",     24909, 24980, 0, 0.025},
		{ "mid_23.08.17",    24985, 25100, 0, 0.025},
		{ "down_25.08.17",   25101, 25212, 0, 0.025},
		{ "down_28.08.17",   25213, 25251, 4, 0.025},
		{ "up_28.08.17",     25252, 25356, 4, 0.025},
		{ "mid_30.08.17",    25357, 25473, 4, 0.025},
		{ "down_01.09.17",   25475, 25659, 4, 0.025},
		{ "up_04.09.17",     25660, 25761, 4, 0.025},
		{ "mid_06.09.17",    25762, 25875, 4, 0.025},
		{ "down_08.09.17",   25877, 25940, 4, 0.025},
		{ "down_09.09.17",   25941, 25993, 0, 0.025},
		{ "down_10.09.17",   25994, 26155, 4, 0.025},
		{ "up_13.09.17",     26157, 26263, 4, 0.025},
		{ "down_15.09.17",   26265, 26436, 4, 0.025},
		{ "up_18.09.17",     26437, 26551, 4, 0.025},
		{ "mid_20.09.17",    26552, 26676, 4, 0.025},
		{ "down_22.09.17",   26678, 26852, 4, 0.025},
		{ "up_25.09.17",     26853, 26958, 4, 0.025},
		{ "mid_27.09.17",    26959, 27061, 4, 0.025},
		{ "down_29.09.17",   27062, 27244, 4, 0.025},
		{ "up_02.10.17",     27245, 27349, 4, 0.025},
		{ "mid_04.10.17",    27350, 27469, 4, 0.025},
		{ "down_06.10.17",   27470, 27643, 4, 0.025},
		{ "up_09.10.17",     27644, 27750, 4, 0.025},
		{ "mid_11.10.17",    27751, 27861, 4, 0.025},
		{ "down_13.10.17",   27862, 28036, 4, 0.025},
		{ "up_16.10.17",     28037, 28145, 4, 0.025},
		{ "mid_18.10.17",    28146, 28283, 4, 0.025},
		{ "down_20.10.17",   28284, 28470, 4, 0.025},
		{ "up_23.10.17",     28471, 28532, 4, 0.025},
		{ "mid_25.10.17",    28533, 28658, 4, 0.025},
		{ "down_27.10.17",   28659, 28858, 4, 0.025},
		{ "up_30.10.17",     28859, 28964, 4, 0.025},
		{ "mid_01.11.17",    28965, 29093, 4, 0.025},
		{ "down_03.11.17",   29097, 29279, 4, 0.025},
		{ "up_06.11.17",     29280, 29392, 4, 0.025},
		{ "mid_08.11.17",    29393, 29522, 4, 0.025},
		{ "down_10.11.17",   29523, 29628, 4, 0.025},
		{ "up_13.11.17",     29701, 29803, 4, 0.025},
		{ "mid_15.11.17",    29804, 29933, 4, 0.025},
		{ "down_17.11.17",   29934, 30107, 4, 0.025},
		{ "up_20.11.17",     30108, 30130, 4, 0.025},
		{ "down_24.11.17",   30331, 30498, 4, 0.025},
		{ "up_27.11.17",     30499, 30609, 4, 0.025},
		{ "mid_29.11.17",    30610, 30723, 4, 0.025},
		{ "down_01.12.17",   30724, 30895, 4, 0.025},
		{ "up_04.12.17",     30896, 30998, 4, 0.025},
		{ "mid_06.12.17",    30999, 31117, 4, 0.025},
		{ "down_08.12.17",   31118, 31296, 4, 0.025},
		{ "up_11.12.17",     31297, 31403, 4, 0.025},
		{ "mid_13.12.17",    31404, 31513, 4, 0.025},
		{ "down_15.12.17",   31514, 31697, 4, 0.025},
		{ "mid_18.12.17",    31698, 31802, 4, 0.025},
		{ "up_20.12.17",     31803, 31912, 4, 0.025},
		{ "down_22.12.17",   31913, 32080, 4, 0.025},
		{ "mid_25.12.17",    32081, 32185, 4, 0.025},
		{ "up_27.12.17",     32186, 32945, 4, 0.025},
		{ "mid_10.01.18",    32946, 33063, 4, 0.025},
		{ "down_12.01.18",   33064, 33255, 4, 0.025},
		{ "up_15.01.18",     33256, 33370, 4, 0.025},
		{ "mid_17.01.18",    33371, 33483, 4, 0.025},
		{ "down_19.01.18",   33484, 33667, 4, 0.025},
		{ "up_22.01.18",     33668, 33766, 4, 0.025},
		{ "mid_24.01.18",    33767, 33799, 4, 0.025},
		{ "down_26.01.18",   33800, 33976, 4, 0.025},
		{ "up_29.01.18",     33977, 34095, 4, 0.025},
		{ "mid_31.01.18",    34096, 34209, 4, 0.025},
		{ "down_02.02.18",   34211, 34385, 4, 0.025},
		{ "up_05.02.18",     34386, 34496, 4, 0.025},
		{ "mid_07.02.18",    34497, 34608, 4, 0.025},
		{ "down_09.02.18",   34609, 34784, 4, 0.025},
		{ "up_12.02.18",     34785, 34894, 4, 0.025},
		{ "mid_14.02.18",    34895, 35011, 4, 0.025},
		{ "down_16.02.18",   35012, 35204, 4, 0.025},
		{ "up_19.02.18",     35205, 35308, 4, 0.025},
		{ "down_23.02.18",   35449, 35643, 4, 0.025},
		{ "up_26.02.18",     35644, 35755, 4, 0.025},
		{ "mid_28.02.18",    35756, 35877, 4, 0.025},
		{ "down_02.03.18",   35878, 36102, 4, 0.025},
		{ "up_05.03.18",     36103, 36231, 4, 0.025},
		{ "mid_07.03.18",    36232, 36402, 4, 0.025},
		{ "down_10.03.18",   36403, 36408, 4, 0.025},
		{ "down_04.05.18",   38048, 38265, 5, 0.025},
		{ "up_07.05.18",     38266, 38398, 5, 0.025},
		{ "mid_09.05.18",    38399, 38535, 5, 0.025},
		{ "down_11.05.18",   38536, 38807, 5, 0.025},
		{ "up_15.05.18",     38808, 38941, 5, 0.025},
		{ "mid_17.05.18",    38942, 39067, 5, 0.025},
		{ "down_19.05.18",   39068, 39299, 5, 0.025},
		{ "up_22.05.18",     39300, 39435, 5, 0.025},
		{ "mid_24.05.18",    39436, 39505, 5, 0.025},
		{ "down_26.05.18",   39506, 39720, 5, 0.025},
		{ "up_28.05.18",     39721, 39848, 5, 0.025},
		{ "mid_30.05.18",    39849, 39982, 5, 0.025},
		{ "down_01.06.18",   39983, 40189, 5, 0.025}
	};
	const char periods[][30] = {"April-June 16", "October 16 - February 17", "March-July 17", "August 17 - January 18"};
