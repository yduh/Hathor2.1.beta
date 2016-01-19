#!/usr/bin/python
from ROOT import gStyle, TCanvas, TPad
import ROOT as r


c1 = r.TCanvas("c1","comparison", 600,500)
gStyle.SetOptStat(0)
gStyle.SetLabelSize(2)

outfile = "./plots_mtt.pdf"

f1 = r.TFile("./../OUTPUT/reweighing_346_1000.root")
f2 = r.TFile("./../OUTPUT/reweighing2_346_1000.root") # 346 2000 here is just for a the dely plot with a cut, doesn't matter for 2D and overall ratio

#h1 = f1.Get("TRUTH/right_tt_M")
#h2 = f2.Get("TRUTH/right_tt_M")
h1 = f1.Get("test2")
h2 = f2.Get("test2")
#h1 = h1pre.Rebin(5)
#h2 = h2pre.Rebin(5)

c1.Print(outfile+"[")
h1.SetLineWidth(3)
h2.SetLineWidth(3)
h1.SetLineColor(r.kRed)
h2.SetLineColor(r.kBlue)
#h1.GetYaxis().SetTitleOffset(1.2)
h2.GetXaxis().SetRangeUser(-3.8, 3.8)
#h2.SetMaximum(0.03)
#h2.SetMinimum(-0.02)
h2.Draw()
h1.Draw("same")
h1.GetXaxis().SetTitle("M_{tt} (GeV)")
#h1.GetYaxis().SetTitle("(dW/dM)/(dLO/dM)")

leg = r.TLegend(0.65, 0.75, 0.85, 0.85)
leg.AddEntry(h1, "1X yukawa const", "l")
leg.AddEntry(h2, "2X", "l")
leg.Draw()

c1.Print(outfile)


