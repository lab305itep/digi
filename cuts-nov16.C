{
	TCut cX("PositronX[0] < 0 || (PositronX[0] > 2 && PositronX[0] < 94)");
	TCut cY("PositronX[1] < 0 || (PositronX[1] > 2 && PositronX[1] < 94)");
	TCut cZ("PositronX[2] > 3.5 && PositronX[2] < 95.5");
        TCut cGamma("AnnihilationEnergy < 1 && AnnihilationGammas < 5 && AnnihilationGammas > 0");
        TCut cVeto("VetoHits == 0");
        TCut cE("PmtCleanEnergy+SiPmCleanEnergy > 10 && PmtCleanEnergy+SiPmCleanEnergy < 20");
}
