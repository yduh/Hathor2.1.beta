#!/usr/bin/python
from ROOT import TCanvas, TGraph, TH1F, gStyle
import ROOT as r
from array import array
import sys


s_low = 346
s_high = 3000
y_low = -5.69 #-4.88
y_high = 5.69 #4.88
s_step = 1
y_step = 0.01

# in INPUT .root
#s_step = 1
#y_step = 0.01

Num_ybins = int((y_high-y_low)/y_step +0.5) #976 #[-4.88, 4.88]/0.01
Num_sbins = int((s_high-s_low)/s_step +0.5) #1654 #[346, 2000]/1
mergebins = 10 # should only for the PDF



PDFcomp = ["hgg", "huub", "hddb", "hccb", "hssb", "hbbb"]
hist = {"hgg":"LOXS_gg", "huub":"LOXS_qq", "hddb":"LOXS_qq", "hccb":"LOXS_qq", "hssb":"LOXS_qq", "hbbb":"LOXS_qq"}
hist_EW = {"hgg":"WXS_gg", "huub":"WXSup_qq", "hddb":"WXSdown_qq", "hccb":"WXSup_qq", "hssb":"WXSdown_qq", "hbbb":"WXSdown_qq"}

inputPDFfile = "~yduh/local/Hathor-2.1-beta/readPDF/INPUT_auto/partonxsec_0.root"
PDFfile = r.TFile(inputPDFfile, "READ")

inputCroSecfile = "%s" %sys.argv[1]
CroSecfile = r.TFile(inputCroSecfile, "READ")

outputReweighfile = "%s" %sys.argv[2]

convolution = []
convolution_EW = []
convert_convolution = []
convert_convolution_EW = []
normal_convolution = []
normal_convolution_EW = []
tot=[]


#___________________________________________________________________||
# start to read the 2D histograms and save in a 3D list in the end
# read the 1D partonic histograms too
# multiply them together bin by bin

for pdfcomp in PDFcomp:
	PDF = PDFfile.Get(pdfcomp)
	rebinPDF = PDF.Rebin(mergebins)
	CroSec = CroSecfile.Get(hist[pdfcomp])
	CroSec_EW = CroSecfile.Get(hist_EW[pdfcomp])

	# sum tot as the normalization factor #
	for s in range(PDF.GetNbinsX()/mergebins):
		tot.append(rebinPDF.GetBinContent(s+1))


	# double for loop to save bin contents in list #
	pdfList = []
	crosecList = []
	crosecList_EW = []

	for y in range(Num_ybins):
		crosecList_yi = []
		crosecList_yi_EW = []
		pdfList_yi = []

		for s in range(Num_sbins):
			pdf = rebinPDF.GetBinContent(int(s/mergebins +(1/mergebins)/2) +(s_low/mergebins+1))
			# GetBinContent(36 repeat 10 times, 37 repeat, ..., 200 repeat)
			pdfList_yi.append(pdf) # [P-s1, P-s2, P-s3, ...]-yi,pdfi
			crosec_yi = CroSec.GetBinContent(s+1, y+1) # [s1, s2, s3, ...]-yi,pdfi
			crosec_yi_EW = CroSec_EW.GetBinContent(s+1, y+1)
			crosecList_yi.append(crosec_yi)
			crosecList_yi_EW.append(crosec_yi_EW)
		crosecList.append(crosecList_yi) # [[s1, s2, s3, ...]-y1, [s1, s2, s3, ...]-y2, [s1, s2, s3, ...]-y3, ...]-pdfi
		crosecList_EW.append(crosecList_yi_EW)
		pdfList.append(pdfList_yi) # [[P-s1, P-s2, P-s3, ...]-y1, [P-s1, P-s2, P-s3, ...]-y2, [P-s1, P-s2, P-s3, ...]-y3]-pdfi


	# multiply bin contents list with PDF lsit #
	convolution_yi = []
	convolution_yi_EW = []

	for y in range(Num_ybins):
		smallList_yi = []
		smallList_yi_EW = []

		# the convolution list for LO #
		for a, b in zip(crosecList[y], pdfList[y]):
			smallList_yi.append(a*b) # [s1*P-s1, s2*P-s2, s3*P-s3, ...]-yi,pdfi
		convolution_yi.append(smallList_yi) # in order to make like this [[ ], [ ], [ ], ...]

		# the convolution list for EW #
		for a, b in zip(crosecList_EW[y], pdfList[y]):
			smallList_yi_EW.append(a*b)
		convolution_yi_EW.append(smallList_yi_EW)

	convolution.append(convolution_yi) # [[s1*P-s1, s2*P-s2, s3*P-s3, ...]-y1, [s1*P-s1, s2*P-s2, s3*P-s3, ...]-y2, ...]-pdfi
	convolution_EW.append(convolution_yi_EW)


# After running the 3 layers (pdf, y, s) for loop, the list is [ [[ ], [ ], [ ], ...], [[ ], [ ], [ ], ...], [[ ], [ ], [ ], ...], ...]
# In the form of [ pdf1, pdf2, pdf3, ... ], where pdfi = [[s1, s2, s3, ...]-y1, [s1, s2, s3, ...]-y2, [s1, s2, s3, ...]-y3, ...]
PDFfile.Close()
CroSecfile.Close()



norpdf = sum(tot)
#___________________________________________________________________||
# convert the saved 3D list

# first convert as [y1, y2, y3, ...] where yi = [[s1, s2, s3, ...]-pdf1, [s1, s2, s3, ...]-pdf2, [s1, s2, s3, ...]-pdf3, ...]
convert1_convolution = zip(*convolution)
convert1_convolution_EW = zip(*convolution_EW)

# second convert as [[pdf1, pdf2, pdf3, ...]-s1, [pdf1, pdf2, pdf3, ...]-s2, [pdf1, pdf2, pdf3, ...]-s3, ...] in different yi
for y in range(Num_ybins):
	convert2_convolution = zip(*convert1_convolution[y])
	convert2_convolution_EW = zip(*convert1_convolution_EW[y])
	convert_convolution.append(convert2_convolution) # need it, otherwise will only have the last y
	convert_convolution_EW.append(convert2_convolution_EW)


#_______________________________________________________________________||
# sum and normalize

for y in range(Num_ybins):
	normal_convolution_yi = []
	normal_convolution_yi_EW = []

	for s in range(Num_sbins):
		normal_convolution_yi.append(sum(convert_convolution[y][s])/norpdf)
		normal_convolution_yi_EW.append(sum(convert_convolution_EW[y][s])/norpdf)
	normal_convolution.append(normal_convolution_yi) # save as [[s1', s2', s3', ...]-y1, [s1', s2', s3', ...]-y2, [s1', s2', s3']-y3, ...]
	normal_convolution_EW.append(normal_convolution_yi_EW)



#________________________________________________________________________||
# fill the well disposal 3D list to the 2D histogram

reweighfile = r.TFile(outputReweighfile, "recreate")
LOresult = r.TH2F("LO", "", Num_sbins, s_low, s_high, Num_ybins, y_low, y_high)
EWresult = r.TH2F("EW", "", Num_sbins, s_low, s_high, Num_ybins, y_low, y_high)
Rresult = r.TH2F("EWtoLO", "", Num_sbins, s_low, s_high, Num_ybins, y_low, y_high)
Rresult_s = r.TH1F("EWtoLO_s", "", Num_sbins, s_low, s_high)
Rresult_y = r.TH1F("EWtoLO_y", "", Num_ybins, y_low, y_high)

Rresult_y_cut = r.TH1F("EWtoLO_y_wcut", "", Num_ybins, y_low, y_high)


for y in range(Num_ybins):
	for s in range(Num_sbins):
		LOresult.SetBinContent(s, y, normal_convolution[y][s])
		EWresult.SetBinContent(s, y, normal_convolution_EW[y][s])

		if(normal_convolution[y][s]!=0):
			Rresult.SetBinContent(s, y, normal_convolution_EW[y][s]/normal_convolution[y][s])

c = r.TCanvas()

LOprojX = LOresult.ProjectionX()
LOprojY = LOresult.ProjectionY()
EWprojX = EWresult.ProjectionX()
EWprojY = EWresult.ProjectionY()
LOprojY_cut= LOresult.ProjectionY("test", LOprojX.FindBin(int(sys.argv[3])), LOprojX.FindBin(int(sys.argv[4])))
EWprojY_cut= EWresult.ProjectionY("test2", EWprojX.FindBin(int(sys.argv[3])), EWprojX.FindBin(int(sys.argv[4])))

EWprojX.Divide(LOprojX)
EWprojY.Divide(LOprojY)
EWprojY_cut.Divide(LOprojY_cut)
Rresult_s = EWprojX
Rresult_y = EWprojY
Rresult_y_cut = EWprojY_cut
#Rresult_s.SetName("EWtoLO_s")
#Rresult_y.SetName("EWtoLO_y")
#Rresult_y_cut.SetName("EWtoLO_y_wcut")

LOresult.Write()
EWresult.Write()
Rresult.Write()
Rresult_s.Write()
Rresult_y.Write()
Rresult_y_cut.Write()

reweighfile.Close()


