#!/usr/bin/python
from ROOT import TCanvas, TGraph, TH1F, gStyle
import ROOT as r
from array import array

PDFcomp = ["hgg", "huub", "hddb", "hccb", "hssb", "hbbb"]
hist = {"hgg":"LOXS_gg_px", "huub":"LOXS_qq_px", "hddb":"LOXS_qq_px", "hccb":"LOXS_qq_px", "hssb":"LOXS_qq_px", "hbbb":"LOXS_qq_px"}

inputPDFfile = "partonxsec_0.root"
PDFfile = r.TFile(inputPDFfile, "READ")

#inputCroSecfile = "../DiffCroSec_final/DiffCroSec.root"
inputCroSecfile = "../DiffCroSec.root"
CroSecfile = r.TFile(inputCroSecfile, "READ")

convolution = []
sum_convolution = []
tot=[]

for pdfcomp in PDFcomp:
	PDF = PDFfile.Get(pdfcomp)
	print hist[pdfcomp]
	CroSec = CroSecfile.Get(hist[pdfcomp])
#for comp, pdfcomp in enumerate(PDFcomp):
#	PDF = PDFfile.Get(pdfcomp)
#	CroSec = CroSecfile.Get("LOXS_px")
	pdfList = []
	crosecList = []

	for i in range(10000):
		tot.append(PDF.GetBinContent(i+1))
	#norpdf = sum(tot)

	for i in range(2000):
		pdf = PDF.GetBinContent(i+351)
		pdfList.append(pdf)
		crosec = CroSec.GetBinContent(i+1)
		crosecList.append(crosec)

	smallList = []
	for a, b in zip(crosecList, pdfList):
		smallList.append(a*b)
	convolution.append(smallList)

#print "convolution = ", convolution
#print tot
norpdf = sum(tot)
convert_convolution = zip(*convolution)
#print "convert_convolution = ", convert_convolution


for i in range(2000):
	sum_convolution.append(sum(convert_convolution[i])/norpdf)
	#print convert_convolution[i]
#print "sum_convolution =", sum_convolution



result = r.TH1F("PDF", "", 1651, 350, 2000)
for j in range(2000):
	result.SetBinContent(j, sum_convolution[j])

c = r.TCanvas()
result.Draw()
