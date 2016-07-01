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
	float upTime;
public:
	HPainter(const char *base);
	HPainter(const char *sname, const char *rname);
	~HPainter(void);
	void Init(const char *sname, const char *rname);
	inline int IsOpen(void) { return tSig && tRand; };
	void Project(TH1 *hist, const char *what, TCut cut);
};

#endif

