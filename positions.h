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
		{ "up_08.06.16",     4508, 4601, 1, 0.1},
		{ "raised_30.09.16", 5540, 5807, 2, 0.05},
		{ "raised_04.10.16", 5808, 5903, 2, 0.025},	// veto corners on
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
		{ "up_04.01.17",     10433, 10832, 3, 0.025},
		{ "mid_11.01.17",    10834, 10973, 3, 0.025},
		{ "down_13.01.17",   10979, 11147, 3, 0.025},
		{ "up_16.01.17",     11150, 11267, 3, 0.025},
		{ "mid_18.01.17",    11271, 11401, 3, 0.025},
		{ "down_20.01.17",   11404, 11563, 3, 0.025},
		{ "up_23.01.17",     11570, 11694, 3, 0.025},
		{ "mid_25.01.17",    11695, 11816, 3, 0.025},
		{ "down_27.01.17",   11820, 12079, 3, 0.025},
		{ "down_01.02.17",   12080, 12421, 3, 0.025},	// there are sources inside this period
		{ "up_06.02.17",     12425, 12547, 3, 0.025},
		{ "mid_08.02.17",    12548, 12671, 3, 0.025},
		{ "down_10.02.17",   12672, 12936, 3, 0.025},	// undefined position inside
		{ "mid_14.02.17",    12938, 12984, 3, 0.025},
		{ "down_15.02.17",   12985, 12994, 3, 0.033},
		{ "up_16.02.17",     13050, 13115, 3, 0.033},
		{ "down_17.02.17",   13117, 13367, 3, 0.033},
		{ "up_20.02.17",     13368, 13504, 3, 0.033},
		{ "mid_22.02.17",    13505, 13663, 3, 0.033},
		{ "down_24.02.17",   13665, 13852, 3, 0.033},
		{ "up_27.02.17",     13854, 13901, 3, 0.033},
		{ "mid_28.02.17",    13902, 14145, 3, 0.033},
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
		{ "up_03.04.17",     17613, 17715, 4, 0.025},
		{ "mid_05.04.17",    17717, 17848, 4, 0.025},
		{ "down_07.04.17",   17849, 18020, 4, 0.025},
		{ "up_10.04.17",     18021, 18131, 4, 0.025},
		{ "mid_12.04.17",    18133, 18263, 4, 0.025},
		{ "down_14.04.17",   18264, 18444, 4, 0.025},
		{ "up_17.04.17",     18445, 18552, 4, 0.025},
		{ "mid_19.04.17",    18553, 18611, 4, 0.025},
		{ "mid_20.04.17",    18612, 18682, 0, 0.025},	// reactor off and getting back
		{ "down_21.04.17",   18684, 18858, 4, 0.025},
		{ "up_24.04.17",     18859, 18921, 4, 0.025},
		{ "mid_26.04.17",    18922, 19040, 4, 0.025},
		{ "down_28.04.17",   19041, 19221, 4, 0.025},
		{ "up_01.05.17",     19222, 19336, 4, 0.025},
//		{ "down_03.05.17",   19337, 19339, 4, 0.025},
		{ "up_03.05.17",     19340, 19448, 4, 0.025},
//		{ "mid_05.05.17",    19449, 19456, 4, 0.025},
		{ "down_05.05.17",   19457, 19638, 4, 0.025},
		{ "up_08.05.17",     19639, 19738, 4, 0.025},
		{ "mid_10.05.17",    19739, 19856, 4, 0.025},
		{ "down_12.05.17",   19857, 20041, 4, 0.025},
		{ "up_15.05.17",     20042, 20151, 4, 0.025},
		{ "down_17.05.17",   20195, 20474, 4, 0.025},	// 22Na and 60Co inside
		{ "up_22.05.17",     20475, 20606, 4, 0.025},
		{ "mid_24.05.17",    20612, 20752, 4, 0.025},
		{ "down_26.05.17",   20754, 20911, 4, 0.025},
		{ "up_29.05.17",     20915, 21027, 4, 0.025},
		{ "mid_31.05.17",    21028, 21147, 4, 0.025},	// position is identified as unknown by Ira.
		{ "down_02.06.17",   21148, 21313, 4, 0.025},
		{ "up_05.06.17",     21315, 21411, 4, 0.025},
		{ "mid_07.06.17",    21412, 21566, 4, 0.025},
		{ "up_12.06.17",     21567, 21672, 4, 0.025},
		{ "mid_14.06.17",    21673, 21800, 4, 0.025},
		{ "down_16.06.17",   21801, 21959, 4, 0.025},
		{ "up_19.06.17",     21960, 22077, 4, 0.025},
		{ "mid_21.06.17",    22079, 22200, 4, 0.025},
		{ "down_23.06.17",   22201, 22380, 4, 0.025},
		{ "up_26.06.17",     22381, 22497, 4, 0.025},
		{ "mid_28.06.17",    22499, 22637, 4, 0.025},
		{ "down_30.06.17",   22638, 22812, 4, 0.025},
		{ "up_03.07.17",     22813, 22924, 4, 0.025},
		{ "mid_05.07.17",    22925, 23041, 4, 0.025},
		{ "down_07.07.17",   23042, 23081, 0, 0.025},	// starting power off
		{ "down_08.07.17",   23082, 23185, 0, 0.025}	// OFF
	};
	const char periods[][30] = {"April-June 16", "October-December 16", "January-March 17", "April-June 17"};

