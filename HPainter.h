#ifndef HPAINTER_H
#define HPAINTER_H

class TChain;

class HPainter {
private:
	TChain *tSig;
	TChain *tRand;
	double upTime;
	unsigned int tBegin;
	unsigned int tEnd;
	int Make_file_list(int *list, int size, int mask, int run_from, int run_to);
	void Init(const char *sname, const char *rname);
	TFile *fRes;
	int ClosefRes;
public:
	HPainter(const char *base);
	HPainter(const char *sname, const char *rname);
	HPainter(int mask, int run_from, int run_to);
	~HPainter(void);
	inline int IsOpen(void) { return tSig && tRand; };
	inline void SetFile(TFile *f) { fRes = f; };
	void OpenFile(const char *name);
	void Project(TH1 *hist, const char *what, TCut cut);
	inline double GetUpTime(void) { return upTime; };
	inline void SetUpTime(double tm) { upTime = tm; };
	void SetUpTime(unsigned int t0, unsigned int t1);
};

#endif

