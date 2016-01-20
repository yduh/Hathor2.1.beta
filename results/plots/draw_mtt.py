#!/usr/bin/python
from ROOT import gStyle, TCanvas, TPad
import ROOT as r


c1 = r.TCanvas("c1","comparison", 600,500)
gStyle.SetOptStat(0)
gStyle.SetLabelSize(2)

outfile = "./plots_mtt.pdf"

f1 = r.TFile("./../reco/yukawa_2Dreweighing/tt_PowhegP8.root")
f2 = r.TFile("./../reco/yukawa2_2Dreweighing/tt_PowhegP8.root")
#f3 = r.TFile("./../reco/noEW/tt_PowhegP8.root")
f3 = r.TFile("./../reco/DATA/DATA.root")

#h1 = f1.Get("TRUTH/right_tt_M")
#h2 = f2.Get("TRUTH/right_tt_M")
h1pre = f1.Get("YUKAWA_RECO/yukawa_Mtt")
h2pre = f2.Get("YUKAWA_RECO/yukawa_Mtt")
#h3pre = f3.Get("YUKAWA_RECO/yukawa_Mtt")
h3pre = f3.Get("RECO/all_tt_M")
h1 = h1pre.Rebin(10)
h2 = h2pre.Rebin(10)
h3 = h3pre.Rebin(10)
h1.Scale(h3.GetEntries()/h1.GetEntries())
h2.Scale(h3.GetEntries()/h2.GetEntries())

c1.Print(outfile+"[")
h1.SetLineWidth(3)
h2.SetLineWidth(3)
h1.SetLineColor(r.kRed)
h2.SetLineColor(r.kBlue)
h3.SetLineColor(r.kBlack)
#h1.GetYaxis().SetTitleOffset(1.2)
h2.GetXaxis().SetRangeUser(346, 800)
h3.Draw()
h2.Draw("same")
h1.Draw("same")
h1.GetXaxis().SetTitle("M_{tt} (GeV)")
#h1.GetYaxis().SetTitle("(dW/dM)/(dLO/dM)")

leg = r.TLegend(0.65, 0.75, 0.85, 0.85)
leg.AddEntry(h1, "1X yukawa const", "l")
leg.AddEntry(h2, "2X", "l")
leg.Draw()

c1.Print(outfile)

'''
h2.Clone()
h2.Divide(h1)
h2.Draw("hist")
c1.Print(outfile)
'''
