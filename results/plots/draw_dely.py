#!/usr/bin/python
from ROOT import gStyle, TCanvas, TPad
import ROOT as r
import sys

c1 = r.TCanvas("c1","comparison", 600,500)
gStyle.SetOptStat(0)
gStyle.SetLabelSize(2)

outfile = "./plots_dely.pdf"

f1 = r.TFile("./../results/yukawa_2Dreweighing/tt_PowhegP8.root")
f2 = r.TFile("./../results/yukawa2_2Dreweighing/tt_PowhegP8.root")

#h1prepre = f1.Get("YUKAWA_GEN/yukawa_delY")
#h2prepre = f2.Get("YUKAWA_GEN/yukawa_delY")
h1prepre = f1.Get("YUKAWA_GEN/yukawa_Mtt_delY")
h2prepre = f2.Get("YUKAWA_GEN/yukawa_Mtt_delY")
Mtt1 = h1prepre.ProjectionX()
h1pre = h1prepre.ProjectionY("test", Mtt1.FindBin(int(sys.argv[1])), Mtt1.FindBin(int(sys.argv[2])))
h2pre = h2prepre.ProjectionY("test2", Mtt1.FindBin(int(sys.argv[1])), Mtt1.FindBin(int(sys.argv[2])))
h1 = h1pre.Rebin(20)
h2 = h2pre.Rebin(20)

c1.Print(outfile+"[")
h1.SetLineWidth(3)
h2.SetLineWidth(3)
h1.SetLineColor(r.kRed)
h2.SetLineColor(r.kBlue)
#h1.GetYaxis().SetTitleOffset(1.2)
h2.GetXaxis().SetRangeUser(float(sys.argv[3]), float(sys.argv[4]))
h2.Draw()
h1.Draw("same")
h1.GetXaxis().SetTitle("M_{tt} (GeV)")
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

