#!/usr/bin/python
from ROOT import gStyle, TCanvas, TPad, TH1D
import ROOT as r
from math import log

Nexp = 500
RecoBins = 100



fSM = r.TFile("~yduh/local/Hathor-2.1-beta/Powheg/results/yukawa_2Dreweighing/tt_PowhegP8.root")
fNP = r.TFile("~yduh/local/Hathor-2.1-beta/Powheg/results/yukawa2_2Dreweighing/tt_PowhegP8.root")
hmtt_SM = fSM.Get("YUKAWA_RECO/yukawa_Mtt")
hmtt_NP = fNP.Get("YUKAWA_RECO/yukawa_Mtt")
mergebins = 15 #int(hmtt_SM.GetXaxis().GetNbins()/RecoBins)
mttmerge_SM = hmtt_SM.Rebin(mergebins)
mttmerge_NP = hmtt_NP.Rebin(mergebins)


histQ_SM = r.TH1D("LikeQ_SM", "SM like Q dist", 40, -900, -750)
histQ_NP = r.TH1D("LikeQ_NP", "NP like Q dist", 40, -900, -750)

def Q(Nsig, modelType):
	NsigList = []
	for exp in range(Nexp):
		NsigList.append(r.gRandom.Poisson(Nsig))

	mttrandom_SM = {}
	mttrandom_NP = {}
	likeQ_SM = []
	likeQ_NP = []
	for exp in range(Nexp):
		mttrandom_SM[exp] = TH1D("SMPseudo%s"%exp, "SM pseudo", RecoBins, 0, 3000)
		mttrandom_SM[exp].FillRandom(mttmerge_SM, NsigList[exp])
		mttrandom_NP[exp] = TH1D("NPPseudo%s"%exp, "NP pseudo", RecoBins, 0, 3000)
		mttrandom_NP[exp].FillRandom(mttmerge_NP, NsigList[exp])

		logL_SM = 0
		logL_NP = 0
		for mtt in range(300,1000, 30):
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



for exp in range(Nexp):
	histQ_SM.Fill(Q(30000, 'SM')[exp])
	histQ_NP.Fill(Q(30000* 1.0050307907, 'NP')[exp])



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
histQ_SM.Draw()
histQ_NP.Draw("same")
c1.cd()




