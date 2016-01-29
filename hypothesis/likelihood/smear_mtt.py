#!/usr/bin/python
from ROOT import gStyle, TCanvas, TPad, TH1D
import ROOT as r
from math import log
from math import sqrt
from SCALE import *

Uncert = [0] #[0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15] #[0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2] #0.07
#TimesLumi = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
timeslumi = 1
Nexp = 2
RecoBins = 150
EvtRatio = 1.0050307907

OutfileName = "./ensemble_mtt_0.1.root"
fSM_sig = r.TFile("~yduh/local/Hathor-2.1-beta/results/reco/1.0y/tt_PowhegP8.root")
fNP_sig = r.TFile("~yduh/local/Hathor-2.1-beta/results/reco/2.0y/tt_PowhegP8.root")
fDATA = r.TFile("~yduh/local/Hathor-2.1-beta/results/reco/1.0y/DATA.root")
hmtt_SM_sig = fSM_sig.Get("YUKAWA_RECO/yukawa_Mtt")
hmtt_NP_sig = fNP_sig.Get("YUKAWA_RECO/yukawa_Mtt")
hmtt_DATA = fDATA.Get("RECO/all_tt_M")
mergebins = 10 #int(hmtt_SM.GetXaxis().GetNbins()/RecoBins)  40 GeV/bin
mttmerge_SM_sig = hmtt_SM_sig.Rebin(mergebins)
mttmerge_NP_sig = hmtt_NP_sig.Rebin(mergebins)
mttmerge_DATA = hmtt_DATA.Rebin(10)


NumData = mttmerge_DATA.GetEntries()
mttmerge_DATA.Scale(timeslumi)
#mttmerge_SM.Scale(NumData*timeslumi/mttmerge_SM.GetEntries()) #hmtt_SM.Rebin(mergebins)
#mttmerge_NP.Scale(NumData*timeslumi* EvtRatio/mttmerge_NP.GetEntries()) #hmtt_NP.Rebin(mergebins)

sum_others("1.0y")
mttmerge_SM = sum_allcomp(mttmerge_SM_sig, sum_others("1.0y"))
mttmerge_NP = sum_allcomp(mttmerge_NP_sig, sum_others("2.0y"))


def sum_others(subpath):
	path = "~yduh/local/Hathor-2.1-beta/results/reco/%s"%subpath
	print path
	histname = "RECO/all_tt_M"

	fSt = r.TFile("%s/STt.root"%path)
	fTbar = r.TFile("%s/Wtbar.root"%path)
	fT = r.TFile("%s/Wt.root"%path)
	fW = r.TFile("%s/WJets.root"%path)
	fDY = r.TFile("%s/DYJets.root"%path)

	hSt = fSt.Get("%s"%histname)
	hTbar = fTbar.Get("%s"%histname)
	hT = fT.Get("%s"%histname)
	hW = fW.Get("%s"%histname)
	hDY = fDY.Get("%s"%histname)

	hSt.Scale(STtscale)
	hTbar.Scale(WTbarscale)
	hT.Scale(WTscale)
	hW.Scale(Wscale)
	hDY.Scale(DYscale)

	hothers.Clone(hSt)
	hothers.Add(hTbar)
	hothers.Add(hT)
	hothers.Add(hW)
	hothers.Add(hDY)

	hothers.Rebin(10)
	return hothers


def sum_allcomp(hsig, hothers):
	hsig.Scale(ttpowheg)
	hsum = TH1D("%s"%hsum, "hsum", RecoBins, 0, 3000)
	for i in range(hsig.GetNbins()):
		if i< hothers.GetNbins():
			binsum = hsig.GetBinContent(i+1) + hothers.GetBinContent(i+1)
		else:
			binsum = hsig.GetBinContent(i+1)
		hsum.SetBinContent(i+1, binsum)

	return hsum




def EvalQ(Nsig, modelType, sys):
	NsigList = []
	for exp in range(Nexp):
		#NsigList.append(r.gRandom.Poisson(Nsig))
		a = r.gRandom.Poisson(Nsig)
		b = r.gRandom.Gaus(1, sys)
		NsigList.append(int(Nsig*b+0.5))

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
		test = 0
		for mtt in range(310,1010, 20):
			SM_mu_i = mttmerge_SM.GetBinContent(mttmerge_SM.FindFixBin(mtt)) # expected num
			SM_k_i = mttrandom_SM[exp].GetBinContent(mttrandom_SM[exp].FindFixBin(mtt)) # observed num
			NP_mu_i = mttmerge_NP.GetBinContent(mttmerge_NP.FindFixBin(mtt)) # expected num
			NP_k_i = mttrandom_NP[exp].GetBinContent(mttrandom_NP[exp].FindFixBin(mtt)) # observed num

			logL_SM = SM_k_i* (log(NP_mu_i) - log(SM_mu_i)) -(NP_mu_i - SM_mu_i) + logL_SM
			logL_NP = NP_k_i* (log(NP_mu_i) - log(SM_mu_i)) -(NP_mu_i - SM_mu_i) + logL_NP
			#print mtt, NP_mu_i, SM_mu_i, NsigList[exp], " random :", NP_k_i, SM_k_i
		#print 'test =', test
		likeQ_SM.append(-2* logL_SM)
		likeQ_NP.append(-2* logL_NP)
		#print -2*logL_SM, -2*logL_NP


	logL_Data = 0
	for mtt in range(310, 1010, 20):
		SM_mu_i = mttmerge_SM.GetBinContent(mttmerge_SM.FindFixBin(mtt))
		NP_mu_i = mttmerge_NP.GetBinContent(mttmerge_NP.FindFixBin(mtt))
		Data_k_i = mttmerge_DATA.GetBinContent(mttmerge_DATA.FindFixBin(mtt))

		logL_Data = Data_k_i* (log(NP_mu_i) - log(SM_mu_i)) -(NP_mu_i - SM_mu_i) + logL_Data
		#print mtt, Data_k_i
	#print "data = ", -2* logL_Data

	if modelType == 'SM':
		model = likeQ_SM
	if modelType == 'NP':
		model = likeQ_NP
	if modelType == 'DATA':
		model = -2* logL_Data

	return model



dev = [] # for diff uncertainty cases
Edev = []
devsigma = []
Edevsigma = []

def ensemblefit(histSM, histNP, minf, maxf):
	gausSM = r.TF1("Gaussian", "gaus", minf, maxf)
	gausNP = r.TF1("Gaussian", "gaus", minf, maxf)
	histSM.Fit(gausSM)
	histNP.Fit(gausNP)

	gausSM_mean = gausSM.GetParameter(1)
	gausSM_Emean = gausSM.GetParError(1)
	gausSM_width = gausSM.GetParameter(2)
	gausSM_Ewidth = gausSM.GetParError(2)
	gausNP_mean = gausNP.GetParameter(1)
	gausNP_Emean = gausNP.GetParError(1)
	gausNP_width = gausNP.GetParameter(2)
	gausNP_Ewidth = gausNP.GetParError(2)

	dev_temp = abs(gausNP_mean - gausSM_mean)
	dev.append(dev_temp)
	Edev_temp = sqrt(gausSM_Emean**2 + gausNP_Emean**2)
	Edev.append(Edev_temp)
	print "derivation = ", dev, "+/-", Edev

	devsigma.append(dev_temp/gausSM_width)
	Edevsigma.append(sqrt((Edev_temp/dev_temp)**2 + (gausSM_Ewidth/gausSM_width)**2) *(dev_temp/gausSM_width))
	print "sigma = ", devsigma, "+/-", Edevsigma



#################################################################################
histQ_SM = {}
histQ_NP = {}

OutFile = r.TFile(OutfileName, "RECREATE")

for i, uncert in enumerate(Uncert):
	OutFile.mkdir("uncert%s"%str(uncert))
	OutFile.cd("uncert%s"%str(uncert))

	minf = -500*timeslumi #-1300 #Q(30000*timeslumi, 'SM', uncert)[]
	maxf = 500*timeslumi #-400
	Nbins = (maxf - minf)/(20*timeslumi)

	histQ_SM[i] = r.TH1D("LikeQ_SM%s"%str(uncert), "SM like Q dist", Nbins, minf, maxf)
	histQ_NP[i] = r.TH1D("LikeQ_NP%s"%str(uncert), "NP like Q dist", Nbins, minf, maxf)

	for exp in range(Nexp):
		histQ_SM[i].Fill(EvalQ(NumData*timeslumi, 'SM', uncert)[exp])
		histQ_NP[i].Fill(EvalQ(NumData*timeslumi* EvtRatio, 'NP', uncert)[exp])

	ensemblefit(histQ_SM[i], histQ_NP[i], minf, maxf)
	histQ_SM[i].Write()
	histQ_NP[i].Write()

	print Q(NumData*timeslumi, 'DATA', uncert)
	#data_line = r.TLine(Q(NumData*timeslumi, 'DATA', uncert)[exp], 0, Q(NumData*timeslumi, 'DATA', uncert)[exp], sqrt(Nexp))
	#data_line.SetLineColor(r.kRed)
	#data_line.Write()


c1 = r.TCanvas("c1","comparison", 800,1000)
#c1.Divide(1, 2)
gStyle.SetOptStat(0)
gStyle.SetLabelSize(2)


#c1.cd(1)
mttmerge_SM.Draw()
mttmerge_NP.Draw("same")
mttmerge_DATA.Draw("same")
#c1.cd(2)
#histQ_SM.SetLineColor(r.kBlue)
#histQ_NP.SetLineColor(r.kRed)
#histQ_NP.Draw()
#histQ_SM.Draw("same")
#c1.cd()


#c2 = r.TCanvas("c2", "fit", 800, 1000)
#c2.Divide(1, 2)
#c2.cd(1)
#ensemblefit(histQ_SM)
#c2.cd(2)
#ensemblefit(histQ_NP)


OutFile.cd()
SummaryPlot = r.TH1F("summaryplot", "summaryplot", len(Uncert), 0, 0.15)

for i in range(len(Uncert)):
	SummaryPlot.SetBinContent(i+1, devsigma[i])
	SummaryPlot.SetBinError(i+1, Edevsigma[i])
SummaryPlot.Write()


OutFile.Close()
