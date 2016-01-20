#!/usr/bin/python
from ROOT import gStyle, TCanvas, TPad, TH1D
import ROOT as r
from math import log
from math import sqrt

Uncert = [0.1] #[0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2] #0.07
#TimesLumi = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
timeslumi = 1
Nexp = 200
RecoBins = 100
EvtRatio = 1.0050307907

OutfileName = "./ensemble_mtt_0.1.root"

fSM = r.TFile("~yduh/local/Hathor-2.1-beta/results/reco/yukawa_2Dreweighing/tt_PowhegP8.root")
fNP = r.TFile("~yduh/local/Hathor-2.1-beta/results/reco/yukawa2_2Dreweighing/tt_PowhegP8.root")
fDATA = r.TFile("~yduh/local/Hathor-2.1-beta/results/reco/Data/DATA.root")
hmtt_SM = fSM.Get("YUKAWA_RECO/yukawa_Mtt")
hmtt_NP = fNP.Get("YUKAWA_RECO/yukawa_Mtt")
hmtt_DATA = fDATA.Get("RECO/all_tt_M")
mergebins = 15 #int(hmtt_SM.GetXaxis().GetNbins()/RecoBins)
mttmerge_SM = hmtt_SM.Rebin(mergebins)
mttmerge_NP = hmtt_NP.Rebin(mergebins)
mttmerge_DATA = hmtt_DATA.Rebin(15)

#minf = -1100
#maxf = -600
#minf = -4000
#maxf = 1000



dev = [] # for diff uncertainty cases
Edev = []
devsigma = []
Edevsigma = []


def Q(Nsig, modelType, sys):
	NsigList = []
	for exp in range(Nexp):
		#NsigList.append(r.gRandom.Poisson(Nsig))
		a = r.gRandom.Poisson(Nsig)
		b = r.gRandom.Gaus(1, sys)
		#print int(Nsig*b+0.5)
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


	for mtt in range(300, 1000, 30):
		SM_mu_i = mttmerge_SM.GetBinContent(mttmerge_SM.FindFixBin(mtt))
		NP_mu_i = mttmerge_NP.GetBinContent(mttmerge_NP.FindFixBin(mtt))
		Data_k_i = mttmerge_Data.GetBinContent(mttmerge_Data.FindFixBin(mtt))

		logL_Data = Data_k_i* (log(NP_mu_i) - log(SM_mu_i)) + logL_Data


	if modelType == 'SM':
		model = likeQ_SM
	if modelType == 'NP':
		model = likeQ_NP
	if modelType == 'DATA':
		model = -2* logL_Data

	#print likeQ
	#return likeQ
	return model



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
#for timeslumi in TimesLumi:
#	OutFile.mkdir("timeslumi%s"%str(timeslumi))
#	OutFile.cd("timeslumi%s"%str(timeslumi))

	minf = -1300 #Q(30000*timeslumi, 'SM', uncert)[]
	maxf = -400
	Nbins = (maxf - minf)/20

	histQ_SM[i] = r.TH1D("LikeQ_SM%s"%str(uncert), "SM like Q dist", Nbins, minf, maxf)
	histQ_NP[i] = r.TH1D("LikeQ_NP%s"%str(uncert), "NP like Q dist", Nbins, minf, maxf)
#for timeslumi in TimesLumi:
#	histQ_SM = r.TH1D("LikeQ_SM%s"%str(timeslumi), "SM like Q dist", 250, minf, maxf)
#	histQ_NP = r.TH1D("LikeQ_NP%s"%str(timeslumi), "NP like Q dist", 250, minf, maxf)

	for exp in range(Nexp):
		histQ_SM[i].Fill(Q(30000*timeslumi, 'SM', uncert)[exp])
		histQ_NP[i].Fill(Q(30000*timeslumi* EvtRatio, 'NP', uncert)[exp])

	ensemblefit(histQ_SM[i], histQ_NP[i], minf, maxf)
	histQ_SM[i].Write()
	histQ_NP[i].Write()

	data_line = r.TLine(Q(30000*timeslumi, 'DATA', uncert)[exp], 0, Q(30000*timeslumi, 'DATA', uncert)[exp], sqrt(Nexp))
	data_line.SetLineColor(r.kRed)
	data_line.Write()

#c1 = r.TCanvas("c1","comparison", 800,1000)
#c1.Divide(1, 2)
#gStyle.SetOptStat(0)
#gStyle.SetLabelSize(2)


#c1.cd(1)
#mttmerge_SM.Draw()
#mttmerge_NP.Draw("same")
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
SummaryPlot = r.TH1F("summaryplot", "summaryplot", len(Uncert), 0, 0.2)
#SummaryPlot = r.TH1F("summaryplot", "summaryplot", len(TimesLumi), 1, 10)

for i in range(len(Uncert)):
#for i in range(len(TimesLumi)):
	SummaryPlot.SetBinContent(i+1, devsigma[i])
	SummaryPlot.SetBinError(i+1, Edevsigma[i])
SummaryPlot.Write()


OutFile.Close()
