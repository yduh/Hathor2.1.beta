#!/usr/bin/python
from ROOT import gStyle, TCanvas, TPad, TH1, TH2, TAxis
import ROOT as r


c1 = r.TCanvas("c1","comparison", 600,500)
gStyle.SetOptStat(0)
gStyle.SetLabelSize(2)

#outfile = "./Mttless450_dely.pdf"

f1X = r.TFile("./../OUTPUT/reweighing_346_2000.root")
f2X = r.TFile("./../OUTPUT/reweighing2_346_2000.root")

h1X = f1X.Get("EW_px")
h2X = f2X.Get("EW_px")



#c1.Print(outfile+"[")
'''
h1X_projY.SetLineWidth(3)
h2X_projY.SetLineWidth(3)
hnoEW_projY.SetLineWidth(3)
h2X_projY.SetMarkerColor(r.kRed)
h1X_projY.SetMarkerColor(r.kBlue)
hnoEW_projY.SetMarkerColor(r.kBlack)
#h1.GetYaxis().SetTitleOffset(1.2)
h2X_projY.Draw()
hnoEW_projY.Draw("same")
h1X_projY.Draw("same")
#h1.GetXaxis().SetTitle("M_{tt} (GeV)")
#h1.GetYaxis().SetTitle("(dW/dM)/(dLO/dM)")

leg = r.TLegend(0.65, 0.75, 0.85, 0.85)
leg.AddEntry(h2X_projY, "2X yukawa const", "l")
leg.AddEntry(h1X_projY, "1X", "l")
leg.AddEntry(hnoEW_projY, "no EW", "l")
leg.Draw()
'''
#c1.Print(outfile)

#for i in range(h2X.GetYaxis().FindBin(-2), h2X.GetYaxis().FindBin(2)):
#	print h2X_projY.GetBinContent(i), hnoEW_projY.GetBinContent(i)

h1X.Clone()
h2X.Divide(h1X)
h2X.Draw("hist")
h2X.GetXaxis().SetRange(175, 300)

#h1X_projY.Clone()
#h2X_projY.Divide(h1X_projY)
#h2X_projY.Draw()
#c1.Print(outfile)
