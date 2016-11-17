{
	TCut cVeto("gtFromVeto > 60");
	TCut cIso("gtFromPrevious > 45 && gtToNext > 80 && EventsBetween == 0");
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
	TCut cR("Distance < 100 && DistanceZ > -40 && DistanceZ < 40");
	TCut c20("gtDiff > 2");
        TCut cGamma("AnnihilationEnergy < 1.5 && AnnihilationGammas < 9");
        TCut cPe("PositronEnergy > 1");
	TCut cSel = cX && cY && cZ && cR && c20 && cGamma && cPe;
	TCut cSig = cSel && cVeto && cIso;
	TCut cBgnd = cSel && (!cVeto);

//*******************************************
//N  Date   EndTime FinalPosition
//1   4/10   14:00   MID		1475578800
//2   7/10   18:00   DOWN		1475852400
//3  10/10   13:25   UP		1476095100
//4 12/10   11:45    MID		1476261900
//5  14/10   18:15   DOWN		1476458100
//-- 17/10   12:10 reactor OFF
//6  17/10   12:45   UP		1476697500
//7  19/10   14:30   MID		1476876600
//8  21/10   15:30   DOWN		1477053000
//9  24/10   13:00?  UP		1477303200
//10  26/10   19:30  MID		1477499400
//*********************************************/

	TCut cP0("unixTime < 1475577800");
	TCut cP1("unixTime > 1475578800 && unixTime < 1475851400");
	TCut cP2("unixTime < 1476094100 && unixTime > 1475852400");
	TCut cP3("unixTime > 1476095100 && unixTime < 1476260900");
	TCut cP4("unixTime < 1476457100 && unixTime > 1476261900");
	TCut cP5("unixTime > 1476458100 && unixTime < 1476696500");
	TCut cP6("unixTime < 1476875600 && unixTime > 1476697500");
}

