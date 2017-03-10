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
		{ "stuck_03.05.16",  2954, 3788, 1, 0.05},
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
		{ "mid_12.10.16",    6136, 6275, 2, 0.025},
		{ "down_14.10.16",   6278, 6467, 2, 0.025},
		{ "up_17.10.16",     6469, 6570, 2, 0.025},
		{ "mid_21.10.16",    6573, 6582, 2, 0.025},
		{ "down_21.10.16",   6587, 6745, 2, 0.025},
		{ "up_24.10.16",     6757, 6815, 2, 0.025},
		{ "mid_27.10.16",    6842, 6909, 2, 0.025},	// was 6923
		{ "down_28.10.16",   6926, 7095, 2, 0.025},
		{ "up_31.10.16",     7106, 7364, 2, 0.025},
		{ "mid_11.11.16",    7387, 7406, 2, 0.025},
		{ "down_11.11.16",   7418, 7458, 2, 0.025},
		{ "up_14.11.16",     7478, 7579, 2, 0.025},
		{ "mid_16.11.16",    7581, 7717, 2, 0.025},
		{ "down_18.11.16",   7727, 7913, 2, 0.025},
		{ "up_21.11.16",     7922, 8042, 2, 0.025},
		{ "mid_23.11.16",    8048, 8179, 2, 0.025},
		{ "down_25.11.16",   8185, 8353, 2, 0.025},
		{ "up_28.11.16",     8357, 8430, 2, 0.025},
		{ "mid_01.12.16",    8470, 8571, 3, 0.025},
		{ "down_02.12.16",   8574, 8738, 3, 0.025},
		{ "up_05.12.16",     8741, 8869, 3, 0.025},
		{ "mid_07.12.16",    8873, 9009, 3, 0.025},
		{ "up_12.12.16",     9012, 9112, 3, 0.025},
		{ "mid_14.12.16",    9116, 9245, 3, 0.025},
		{ "down_16.12.16",   9253, 9470, 3, 0.025},
		{ "up_19.12.16",     9475, 9600, 3, 0.025},
		{ "mid_21.12.16",    9603, 9712, 3, 0.025},
		{ "down_23.12.16",   9715, 9869, 3, 0.025},
		{ "up_26.12.16",     9871, 10019, 3, 0.025},
		{ "mid_28.12.16",    10021, 10171, 3, 0.025},
		{ "down_30.12.16",   10175, 10307, 3, 0.025},
		{ "down_02.01.17",   10308, 10356, 3, 0.025},
		{ "down_03.01.17",   10357, 10424, 3, 0.025},
		{ "up_04.01.17",     10433, 10832, 3, 0.025},
		{ "mid_11.01.17",    10834, 10973, 3, 0.025},
		{ "down_13.01.17",   10979, 11147, 3, 0.025},
		{ "up_16.01.17",     11150, 11267, 3, 0.025},
		{ "mid_18.01.17",    11271, 11401, 3, 0.025},
		{ "down_20.01.17",   11404, 11563, 3, 0.025},
		{ "up_23.01.17",     11570, 11688, 3, 0.025}
	};
	const char periods[][30] = {"April-June", "October-November", "December-January"};
