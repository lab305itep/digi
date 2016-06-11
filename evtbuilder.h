#ifndef EVTBUILDER_H
#define EVTBUILDER_H

//	Glaobal
#define GLOBALFREQ	125000000.0		// 125 Mhz

/***********************	Types		****************************/

struct DanssEventStruct {
//		Common parameters
	long long	globalTime;		// time in terms of 125 MHz
	long long	number;			// event number in the file
	int		unixTime;		// linux time, seconds
	float		fineTime;		// fine time of the event (for hit selection)
//		Veto parameters
	int		VetoHits;		// hits in the event record
	float		VetoEnergy;		// Energy Sum of all hits
	int		VetoCleanHits;		// hits above threshold and in time window
	float		VetoCleanEnergy;	// Energy Sum of clean hits
//		PMT parameters
	int		PmtHits;		// the same as above for PMT
	float		PmtEnergy;
	int		PmtCleanHits;
	float		PmtCleanEnergy;
//		SiPM parameters
	int		SiPmHits;		// the same as above for PMT
	float		SiPmEnergy;
	int		SiPmCleanHits;
	float		SiPmCleanEnergy;
	int		SiPmEarlyHits;		// to understand random background
	float		SiPmEarlyEnergy;
//		"positron cluster" parameters
	int		PositronHits;		// hits in the cluster
	float		PositronSiPmEnergy;	// Energy sum of the cluster (SiPM)
	float		MaxHitEnergy;		// Energy of the maximum hit (SiPM)
	float		PositronX[3];		// cluster position
	int		AnnihilationGammas;	// number of possible annihilation gammas
	float		AnnihilationEnergy;	// Energy in annihilation gammas
//		"neutron" parameters
	int		NeutronHits;		// number of hits considered as neutron capture gammas
	float		NeutronSiPmEnergy;	// Energy sum of above (SiPM)
	float		NeutronX[3];		// center of gammas position
	float		NeutronGammaEnergy[5];	// sorted list of the 5 most energetic gammas
	float		NeutronGammaDistance[5];	// distances for the gammas above to the "neutron" center
	float		NeutronRadius;		// average distance between hits and the center
};

struct DanssEventStruct2 {
//		Common parameters
	long long	globalTime;		// time in terms of 125 MHz
	long long	number;			// event number in the file
	int		unixTime;		// linux time, seconds
	float		fineTime;		// fine time of the event (for hit selection)
//		Veto parameters
	int		VetoHits;		// hits in the event record
	float		VetoEnergy;		// Energy Sum of all hits
	int		VetoCleanHits;		// hits above threshold and in time window
	float		VetoCleanEnergy;	// Energy Sum of clean hits
//		PMT parameters
	int		PmtHits;		// the same as above for PMT
	float		PmtEnergy;
	int		PmtCleanHits;
	float		PmtCleanEnergy;
//		SiPM parameters
	int		SiPmHits;		// the same as above for PMT
	float		SiPmEnergy;
	int		SiPmCleanHits;
	float		SiPmCleanEnergy;
	int		SiPmEarlyHits;		// to understand random background
	float		SiPmEarlyEnergy;
//		"positron cluster" parameters
	int		PositronHits;		// hits in the cluster
	float		PositronSiPmEnergy;	// Energy sum of the cluster, corrected (SiPM)
	float		PositronPmtEnergy;	// Energy sum of the cluster, corrected (PMT)
	float		MaxHitEnergy;		// Energy of the maximum hit (SiPM)
	float		PositronX[3];		// cluster position
	int		AnnihilationGammas;	// number of possible annihilation gammas
	float		AnnihilationEnergy;	// Energy in annihilation gammas
//		"neutron" parameters
	float		NeutronX[3];		// center of gammas position
	float		NeutronRadius;		// average distance between hits and the center
};

struct DanssInfoStruct {
	long long	upTime;			// running time in terms of 125 MHz
	int		startTime;		// linux start time, seconds
	int		stopTime;		// linux stop time, seconds
	int		events;			// number of events
};

struct DanssPairStruct {
//		Common parameters
	long long	number[2];		// event numbers in the file
	int		unixTime;		// linux time, seconds
	float		SiPmCleanEnergy[2];	// Full Clean energy SiPm
	float		PmtCleanEnergy[2];	// Full Clean energy Pmt
//		"positron cluster" parameters
	int		PositronHits;		// hits in the cluster
	float		PositronEnergy;		// Energy sum of the cluster (SiPM)
	float		MaxHitEnergy;		// Energy of the maximum hit (SiPM)
	float		PositronX[3];		// cluster position
	int		AnnihilationGammas;	// number of possible annihilation gammas
	float		AnnihilationEnergy;	// Energy in annihilation gammas
//		"neutron" parameters
	int		NeutronHits;		// number of hits considered as neutron capture gammas
	float		NeutronSiPmEnergy;	// Energy sum of above (SiPM)
	float		NeutronX[3];		// center of gammas position
        float           NeutronGammaEnergy[5];  // sorted list of the 5 most energetic gammas
        float           NeutronGammaDistance[5];        // distances for the gammas above to the "neutron" center
	float		NeutronRadius;		// average distance between hits and the center
//		Pair parameters
	float		gtDiff;			// time difference in us (from 125 MHz clock)
	float		Distance;		// distance between neutron and positron, cm
	float		DistanceZ;		// in Z, cm
};

struct DanssPairStruct2 {
//		Common parameters
	long long	number[2];		// event numbers in the file
	int		unixTime;		// linux time, seconds
	float		SiPmCleanEnergy[2];	// Full Clean energy SiPm
	float		PmtCleanEnergy[2];	// Full Clean energy Pmt
//		"positron cluster" parameters
	int		PositronHits;		// hits in the cluster
	float		PositronEnergy;		// Energy sum of the cluster (SiPM)
	float		MaxHitEnergy;		// Energy of the maximum hit (SiPM)
	float		PositronX[3];		// cluster position
	int		AnnihilationGammas;	// number of possible annihilation gammas
	float		AnnihilationEnergy;	// Energy in annihilation gammas
//		"neutron" parameters
	int		NeutronHits;		// number of hits considered as neutron capture gammas
	float		NeutronEnergy;		// Energy sum of above (SiPM)
	float		NeutronX[3];		// center of gammas position
	float		NeutronRadius;		// average distance between hits and the center
//		Pair parameters
	float		gtDiff;			// time difference in us between positron and neutron
	float		Distance;		// distance between neutron and positron, cm
	float		DistanceZ;		// in Z, cm
//		Environment
	float		gtFromPrevious;		// time from the previous hit before positron, us
	float		gtToNext;		// time to the next hit after neutron, counted from positron, us
	int		EventsBetween;		// Events between positron and neutron
};

//		248Cm analysis
struct DanssCmStruct {
	long long	number[10];		// event numbers in the file
	int		unixTime;		// linux time, seconds
	int		N;			// number of neutrons + 1
	float		SiPmCleanEnergy[10];	// Full Clean energy SiPm
	float		PmtCleanEnergy[10];	// Full Clean energy Pmt
//		"neutron" parameters
	int		Hits[10];
	int		NeutronHits[10];	// number of hits considered as neutron capture gammas
	float		NeutronEnergy[10];	// Energy sum of above (SiPM)
	float		NeutronX[10][3];	// center of gammas position
	float		PositronX[10][3];	// center of maximum hit clusters
	float		NeutronGammaEnergy[10][5];	// sorted list of the 5 most energetic gammas
	float		NeutronGammaDistance[10][5];	// distances for the gammas above to the "neutron" center
	float		NeutronRadius[10];		// average distance between hits and the center
//		Pair parameters
	float		gtDiff[10];		// time difference in us (from 125 MHz clock)
	float		Distance[10];		// distance between neutron and positron, cm
	float		DistanceZ[10];		// in Z, cm
};

//		MC truth
struct DanssMcStruct {
	float		Energy;
	float		X[3];
};

#endif /* EVTBUILDER_H */

