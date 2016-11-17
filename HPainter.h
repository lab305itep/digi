#ifndef HPAINTER_H
#define HPAINTER_H

class TTree;
class TFile;

class HPainter {
private:
	TFile *fSig;
	TFile *fRand;
	TTree *tSig;
	TTree *tRand;
	double upTime;
	unsigned int tBegin;
	unsigned int tEnd;
public:
	HPainter(const char *base);
	HPainter(const char *sname, const char *rname);
	~HPainter(void);
	void Init(const char *sname, const char *rname);
	inline int IsOpen(void) { return tSig && tRand; };
	void Project(TH1 *hist, const char *what, TCut cut);
	inline double GetUpTime(void) { return upTime; };
	inline void SetUpTime(double tm) { upTime = tm; };
	void SetUpTime(unsigned int t0, unsigned int t1);
};

#endif

