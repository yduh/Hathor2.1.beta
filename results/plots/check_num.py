#!/usr/bin/python
from ROOT import gStyle, TCanvas, TPad, TH1, TAxis
import ROOT as r



#f1 = r.TFile("./../../results/yukawa/mtt/tt_PowhegP8.root")
f1 = r.TFile("./../../results/yukawa2/mtt/tt_PowhegP8.root")
#f1 = r.TFile("./../../results/noEW/tt_PowhegP8.root")

#h1 = f1.Get("TRUTH/right_tt_M")
h1 = f1.Get("TRUTH/right_tt_M")

bin_mttCut = h1.GetXaxis().FindBin(450)
bin_mttMax = h1.GetXaxis().FindBin(2000)
print bin_mttCut, bin_mttMax
totNumless = 0
totNumlarger = 0

for i in range (1, bin_mttCut+1):
	#print i
	totNumless = h1.GetBinContent(i) + totNumless
print totNumless

for j in range(bin_mttCut, bin_mttMax+1):
	#print j
	totNumlarger = h1.GetBinContent(j) + totNumlarger
print totNumlarger


print totNumless/totNumlarger


