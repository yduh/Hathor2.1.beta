#!/usr/bin/python
from ROOT import TCanvas, TGraph, TH1F, gStyle
import ROOT as r
from array import array

Num_ybins = 976
Num_sbins = 2000
mergebins = 10
PDFcomp = ["hgg", "huub", "hddb", "hccb", "hssb", "hbbb"]
hist = {"hgg":"LOXS_gg", "huub":"LOXS_qq", "hddb":"LOXS_qq", "hccb":"LOXS_qq", "hssb":"LOXS_qq", "hbbb":"LOXS_qq"}

inputPDFfile = "partonxsec_0.root"
PDFfile = r.TFile(inputPDFfile, "READ")

inputCroSecfile = "../DiffCroSec.root"
CroSecfile = r.TFile(inputCroSecfile, "READ")

convolution = []
sum_convolution = []
tot=[]


#___________________________________________________________________||
# start to read the 2D histograms and save in a 3D list in the end
# read the 1D partonic histograms too
# multiply them together bin by bin

for pdfcomp in PDFcomp:
	PDF = PDFfile.Get(pdfcomp)
	rebinPDF = PDF.Rebin(mergebins)
	CroSec = CroSecfile.Get(hist[pdfcomp])
#for comp, pdfcomp in enumerate(PDFcomp):
#	PDF = PDFfile.Get(pdfcomp)
#	CroSec = CroSecfile.Get("LOXS_gg")
	pdfList = []
	crosecList = []

	for i in range(10000/mergebins):
		tot.append(rebinPDF.GetBinContent(i+1))
	#norpdf = sum(tot)

	for y in range(Num_ybins):
		crosecList_y = []
		pdfList_y = []
		for i in range(Num_sbins):
			pdf = rebinPDF.GetBinContent(int((i+351)/5 +0.5))
			pdfList_y.append(pdf)
			crosec_y = CroSec.GetBinContent(i+1, y+1) #[[x1, x2, x3, ...]-y1, [x1, x2, x3, ...]-y2, []-y3, ...]
			crosecList_y.append(crosec_y)
		#print crosecList_y[1:10]
		crosecList.append(crosecList_y)
		pdfList.append(pdfList_y)

	#print crosecList
	#print pdfList
	convolution_y = []
	for y in range(Num_ybins):
		smallList_y = []
		for a, b in zip(crosecList[y], pdfList[y]):
			smallList_y.append(a*b)
		convolution_y.append(smallList_y)
		#print convolution_y
	convolution.append(convolution_y)

PDFfile.Close()
CroSecfile.Close()

#print "convolution = ", convolution
#print tot
norpdf = sum(tot)
convert_convolution = zip(*convolution)
#print "convert_convolution = ", convert_convolution


convertdouble_convolution = []
for y in range(Num_ybins):
	convert2_convolution = zip(*convert_convolution[y])
	convertdouble_convolution.append(convert2_convolution)
#print "convertdouble_convolution = ", convertdouble_convolution


#_______________________________________________________________________||
# sum over and normalize to all the PDFs

for y in range(Num_ybins):
	sum_convolution_y = []
	for i in range(Num_sbins):
		sum_convolution_y.append(sum(convertdouble_convolution[y][i])/norpdf)
	sum_convolution.append(sum_convolution_y)
#print "sum_convolution =", sum_convolution


#________________________________________________________________________||
# fill the well disposal 3D list to the 2D histogram

outputreweighfile = r.TFile("reweighing.root", "recreate")
LOresult = r.TH2F("LO", "", 1650, 350, 2000, Num_ybins, -4.88, 4.88)

for y in range(Num_ybins):
	for j in range(Num_sbins):
		LOresult.SetBinContent(j, y, sum_convolution[y][j])

#c = r.TCanvas()
#result.Draw("colz")
#projX = result.ProjectionX()
#projX.Draw()
#projY = result.ProjectionY()
#projY.Draw()

result.Write()
outputreweighfile.Close()


