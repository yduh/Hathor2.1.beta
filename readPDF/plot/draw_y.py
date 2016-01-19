#!/usr/bin/python
from ROOT import gStyle, TCanvas, TPad
import ROOT as r


c1 = r.TCanvas("c1","comparison", 600,500)
gStyle.SetOptStat(0)
gStyle.SetLabelSize(2)

outfile = "./plots_mtt.pdf"

f1 = r.TFile("./../../results/yukawa_2Dreweighing/tt_PowhegP8.root")
f2 = r.TFile("./../../results/yukawa2_2Dreweighing/tt_PowhegP8.root")

h1pre = f1.Get("YUKAWA_GEN/yukawa_Y")
h2pre = f2.Get("YUKAWA_GEN/yukawa_Y")
h1 = h1pre.Rebin(8)
h2 = h2pre.Rebin(8)


c1.Print(outfile+"[")
h1.SetLineWidth(3)
h2.SetLineWidth(3)
h1.SetLineColor(r.kRed)
h2.SetLineColor(r.kBlue)
h2.GetXaxis().SetRangeUser(-2.4, 2.4)
#h1.GetYaxis().SetTitleOffset(1.2)
h2.Draw()
h1.Draw("same")
#h1.GetXaxis().SetTitle("")
#h1.GetYaxis().SetTitle("(dW/dM)/(dLO/dM)")

leg = r.TLegend(0.65, 0.75, 0.85, 0.85)
leg.AddEntry(h1, "1X yukawa const", "l")
leg.AddEntry(h2, "2X", "l")
leg.Draw()

c1.Print(outfile)


h2.Clone()
h2.Divide(h1)
h2.Draw("hist")
h2.SetMarkerSize(3)
c1.Print(outfile)

