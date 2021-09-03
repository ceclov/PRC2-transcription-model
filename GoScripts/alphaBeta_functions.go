package main

import (
	"fmt"
	"os"
	"runtime"
)

func Fig3b() {
	numCPU = runtime.NumCPU() - 1

	saveParams := false
	importParams := false
	text := "notFile"

	var alphas, betas, fmaxs, pdems []float64

	pdem = 4e-3
	pex = 1e-3
	fmax = 4e-3

	isAlphaRun := false
	isBetaRun := true

	newAlpha := float64(100)

	if !importParams {

		alphas = getParam(0.25, 10, 50, true)

		betas = getParam(0.25, 10, 50, true)

		fmaxs = []float64{fmax}

		pdems = []float64{pdem}

		alphas, betas, fmaxs, pdems = paramsFour(alphas, betas, fmaxs, pdems)

	}

	if saveParams {
		divideParams(alphas, betas, fmaxs, pdems, 30)
		fmt.Println(len(alphas), len(betas), len(fmaxs))
	}

	if importParams {

		if len(os.Args) < 2 {
			fmt.Println("Missing parameter, provide file name!")
			return
		}
		data := os.Args[1]

		text = string(data)

		alphas, betas, fmaxs, pdems = importJsonParam(text)

	}

	if !saveParams {
		currSim := make(Vector, 3*len(alphas))

		defineParams()

		fmt.Println("alpha", "\t", "beta", "\t", "fmax", "\t", "me3", "\t", "transcription", "\t", "type", "\t", "Pdem", "\t", "Pex", "\t", "nSim", nSim, "fmin", fmin, "kme", kme, text)

		currSim.runAB(alphas, betas, fmaxs, pdems, newAlpha, isAlphaRun, isBetaRun)
	}

}

func (v Vector) runAB(alphas, betas, fmaxs, pdems []float64, newAlpha float64, isAlphaRun, isBetaRun bool) {
	l := len(alphas)
	var n3 nucleosomes
	var ind int
	var currWild, currDelta, currRescue, transWild, transDelta, transRescue []float64
	var resWild, resDelta, resRescue float64

	for i := 0; i < l; i++ {

		currWild = make([]float64, nSim)
		currDelta = make([]float64, nSim)
		currRescue = make([]float64, nSim)

		transWild = make([]float64, nSim)
		transDelta = make([]float64, nSim)
		transRescue = make([]float64, nSim)

		for k := 0; k < nSim; k++ {

			n3 = nucleosomes{}
			n3 = initSystem(N, 3)
			n3.nCellCycle = 0

			if isAlphaRun {
				n3 = AlphaRun(alphas[i], betas[i], fmaxs[i], pdems[i], newAlpha)
			}

			if isBetaRun {

				n3 = BetaRun(alphas[i], betas[i], fmaxs[i], pdems[i])
			}

			transWild[k] = n3.trans[nGen-1]
			transDelta[k] = n3.trans[(2*nGen)-1]
			transRescue[k] = n3.trans[(3*nGen)-1]

			currWild[k] = n3.Pme3end[nGen-1]
			currDelta[k] = n3.Pme3end[(2*nGen)-1]
			currRescue[k] = n3.Pme3end[(3*nGen)-1]

		}
		ind = i * 3
		v[ind] = meanFloat(currWild)
		v[ind+1] = meanFloat(currDelta)
		v[ind+2] = meanFloat(currRescue)

		resWild = meanFloat(transWild)
		resDelta = meanFloat(transDelta)
		resRescue = meanFloat(transRescue)

		printRes(v[ind], v[ind+1], v[ind+2], resWild, resDelta, resRescue, alphas[i], betas[i], fmax, pdems[i], pex)

	}

}

func (v Vector) printAll(alphas, betas, fmaxs []float64) {
	fmt.Println("alpha", "\t", "beta", "\t", "transcription", "\t", "type", nSim)

	var ind int
	for i := 0; i < len(alphas); i++ {
		ind = i * 3
		printRes(v[ind], v[ind+1], v[ind+2], v[ind], v[ind+1], v[ind+2], alphas[i], betas[i], fmax, pdem, pex)

	}
}

func printRes(res1, res2, res3, res4, res5, res6, alpha, beta, fmax, pdem, pex float64) {
	res := make([]float64, 8)

	res[0] = alpha
	res[1] = beta
	res[2] = fmax

	res[6] = pdem
	res[7] = pex

	res[3] = res1
	res[4] = res4
	res[5] = 0.0
	printFloat64(res)

	res[3] = res2
	res[4] = res5
	res[5] = 1.0
	printFloat64(res)

	res[3] = res3
	res[4] = res6
	res[5] = 2.0
	printFloat64(res)

}

func BetaRun(thisAlpha, thisBeta, thisfmax, thisPdem float64) nucleosomes {
	var n3 nucleosomes
	var t float64

	fmax = thisfmax
	pdem = thisPdem
	pex = 0.25 * pdem

	defineParams()

	n3 = nucleosomes{}
	n3 = initSystem(N, 3)

	n3.trans = make([]float64, nGen)

	n3.Pme3end = make([]float64, nGen)

	n3, _ = GStep(n3, tCC, 20, 0.0, 0)

	n3.Rep()

	n3 = ResetSystem(n3)

	n3.trans = make([]float64, (nGen*2 + DelGen))

	n3.Pme3end = make([]float64, (nGen*2 + DelGen))

	alpha = thisAlpha
	beta = thisBeta

	n3, t = GStep(n3, tCC, nGen, 0.0, 0)

	beta = 0.0
	tCC = tCCvar

	n3, t = GStep(n3, tCCvar, nGen+DelGen, t, nGen)

	tCC = tCCnorm

	beta = thisBeta

	n3, t = GStep(n3, tCC, (nGen*2 + DelGen), t, nGen+DelGen)

	return n3
}

func AlphaRun(thisAlpha, thisBeta, thisfmax, thisPdem, newAlpha float64) nucleosomes {
	var n3 nucleosomes
	var t float64

	fmax = thisfmax
	pdem = thisPdem
	pex = 0.25 * pdem

	defineParams()

	n3 = nucleosomes{}
	n3 = initSystem(N, 3)

	alpha = thisAlpha
	beta = thisBeta

	n3.trans = make([]float64, nGen)
	n3.Pme3end = make([]float64, nGen)

	n3, _ = GStep(n3, tCC, nGen, 0.0, 0)

	n3.Rep()
	n3 = ResetSystem(n3)

	n3.trans = make([]float64, nGen*2+DelGen)
	n3.Pme3end = make([]float64, nGen*2+DelGen)

	n3, t = GStep(n3, tCC, nGen, 0.0, 0)

	alpha = newAlpha

	n3, t = GStep(n3, tCC, (nGen + DelGen), t, nGen)

	alpha = thisAlpha

	n3, t = GStep(n3, tCC, (nGen*2 + DelGen), t, (nGen + DelGen))

	return n3
}
