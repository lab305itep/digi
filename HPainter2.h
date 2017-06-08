#ifndef HPAINTER2_H
#define HPAINTER2_H

class TChain;

class HPainter2 {
private:
	TChain *tSig;
	TChain *tRand;
	double upTime;
	unsigned int tBegin;
	unsigned int tEnd;
	int Make_file_list(int *list, int size, int mask, int run_from, int run_to);
	TFile *fRes;
	int ClosefRes;
public:
	HPainter2(int mask, int run_from, int run_to, const char *root2dir = "/space/danss_pair/");
	~HPainter2(void);
	inline int IsOpen(void) { return tSig && tRand; };
	inline void SetFile(TFile *f) { fRes = f; };
	void OpenFile(const char *name);
	void Project(TH1 *hist, const char *what, TCut cut);
	inline double GetUpTime(void) { return upTime; };
	inline void SetUpTime(double tm) { upTime = tm; };
	void SetUpTime(unsigned int t0, unsigned int t1);
	inline TChain *GetPairChain(void) { return tSig;};
	inline TChain *GetRandomChain(void) { return tRand;};
};

#endif

