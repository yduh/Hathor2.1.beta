#!/usr/bin/python
from ROOT import gStyle, TCanvas, TPad, TH1D
import ROOT as r
from math import log

Uncert = 0.07
Nexp = 500
RecoBins = 100
EvtRatio = 1.0050307907


fSM = r.TFile("~yduh/local/Hathor-2.1-beta/Powheg/results/yukawa_2Dreweighing/tt_PowhegP8.root")
fNP = r.TFile("~yduh/local/Hathor-2.1-beta/Powheg/results/yukawa2_2Dreweighing/tt_PowhegP8.root")
hmtt_SM = fSM.Get("YUKAWA_RECO/yukawa_delY")
hmtt_NP = fNP.Get("YUKAWA_RECO/yukawa_delY")
#print int(hmtt_SM.GetXaxis().GetNbins()/RecoBins)
mergebins = 12 #int(hmtt_SM.GetXaxis().GetNBins()/RecoBins)
mttmerge_SM = hmtt_SM.Rebin(mergebins)
mttmerge_NP = hmtt_NP.Rebin(mergebins)


histQ_SM = r.TH1D("LikeQ_SM", "SM like Q dist", 25, -150, -50)
histQ_NP = r.TH1D("LikeQ_NP", "NP like Q dist", 25, -150, -50)
'''
def smear():
	Gaus_i = {} #pdfbin, mttbin
	for mtt in range(300, 1000, 30):
		SM_PDF_i = mttmerge_SM.GetBinContent(mttmerge_SM.FindFixBin(mtt))
		Eval()
		Gaus_i[pdfbin] = SM_PDF_i* Gaus(SM_PDF_i, 0.12)

	element = []
	for mtt in range(300, 1000, 30):
		element.append(Gaus_i[pdfbin].Eval(mtt))
		smear_SM_PDF.SetBinContent(mttmerge_SM.FindFixBin(mtt), SM_PDF_i* Gaus(SM_PDF_i, 0.12))
'''

def Q(Nsig, modelType):
	NsigList = []
	for exp in range(Nexp):
		#NsigList.append(r.gRandom.Poisson(Nsig))
		a = r.gRandom.Poisson(Nsig)
		b = r.gRandom.Gaus(1, Uncert)
		NsigList.append(int(Nsig*b+0.5))

	mttrandom_SM = {}
	mttrandom_NP = {}
	likeQ_SM = []
	likeQ_NP = []
	for exp in range(Nexp):
		mttrandom_SM[exp] = TH1D("SMPseudo%s"%exp, "SM pseudo", RecoBins, -6, 6)
		mttrandom_SM[exp].FillRandom(mttmerge_SM, NsigList[exp])
		mttrandom_NP[exp] = TH1D("NPPseudo%s"%exp, "NP pseudo", RecoBins, -6, 6)
		mttrandom_NP[exp].FillRandom(mttmerge_NP, NsigList[exp])

		logL_SM = 0
		logL_NP = 0
		for mtt in range(-4,4, 1):
			SM_mu_i = mttmerge_SM.GetBinContent(mttmerge_SM.FindFixBin(mtt)) # expected num
			SM_k_i = mttrandom_SM[exp].GetBinContent(mttrandom_SM[exp].FindFixBin(mtt)) # observed num
			NP_mu_i = mttmerge_NP.GetBinContent(mttmerge_NP.FindFixBin(mtt)) # expected num
			NP_k_i = mttrandom_NP[exp].GetBinContent(mttrandom_NP[exp].FindFixBin(mtt)) # observed num

			#logL_SM = SM_k_i* log(SM_mu_i) - SM_mu_i + logL_SM
			#logL_NP = NP_k_i* log(NP_mu_i) - NP_mu_i + logL_NP
			#print logL_SM, logL_NP
			logL_SM = SM_k_i* (log(NP_mu_i) - log(SM_mu_i)) + logL_SM
			logL_NP = NP_k_i* (log(NP_mu_i) - log(SM_mu_i)) + logL_NP
		likeQ_SM.append(-2* logL_SM)
		likeQ_NP.append(-2* logL_NP)

	if modelType == 'SM':
		model = likeQ_SM
	if modelType == 'NP':
		model = likeQ_NP

	#print likeQ
	#return likeQ
	return model

def ensemblefit(hist):
	minf = -1100
	maxf = -600
	gaus = r.TF1("Gaussian", "gaus", minf, maxf)
	hist.Fit(gaus)
	#hist.SetParameters(1, 0, 1)
	#print hist.GetParameter(0), histGetParError(0)
	#print hist.GetParameter(1), histGetParError(1)


for exp in range(Nexp):
	histQ_SM.Fill(Q(30000, 'SM')[exp])
	histQ_NP.Fill(Q(30000* EvtRatio, 'NP')[exp])



c1 = r.TCanvas("c1","comparison", 800,1000)
c1.Divide(1, 2)
gStyle.SetOptStat(0)
gStyle.SetLabelSize(2)



c1.cd(1)
mttmerge_SM.Draw()
mttmerge_NP.Draw("same")
c1.cd(2)
#mttrandom_SM[0].Draw()
#mttrandom_NP[0].Draw("same")
#c1.cd(3)
histQ_SM.SetLineColor(r.kBlue)
histQ_NP.SetLineColor(r.kRed)
histQ_NP.Draw()
histQ_SM.Draw("same")
c1.cd()



c2 = r.TCanvas("c2", "fit", 800, 1000)
c2.Divide(1, 2)
c2.cd(1)
ensemblefit(histQ_SM)
c2.cd(2)
ensemblefit(histQ_NP)
