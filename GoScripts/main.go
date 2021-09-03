package main

import (
	"fmt"
	"math/rand"

	"github.com/seehuhn/mt19937"
)

//pre-initalizations of variables. Values might be redefined.
var kme = 8e-6
var pdem = 4e-3
var fmax = 4e-3

var kme01 = 9.0 * kme
var kme12 = 6.0 * kme
var gme01 = kme01 / 20.0
var gme12 = kme12 / 20.0
var gme23 = kme / 20.0

const pme2 = 0.1

var fmin = 1e-04

const PT = 1.0 / 3.0

var gdem = fmin * pdem

var N = 30

var tCC = 22.0
var tCCvar = 22.0
var tCCnorm = 22.0

var alpha = 1.0

var beta = 1.0

var pex = 1e-3

var count0 = false
var count2 = true
var count3 = true
var count1 = false

var usePT = true

var rng = rand.New(mt19937.New())

var nGen = 20
var DelGen = 20
var nSim = 1000
var timeRec = 21.0

var transMax = 1 / 60.0

func main() {

	//Fig3b()
	//Fig3d()
	//FigS7d()
	FigS10d()
}

func defineParams() {
	nSim = 1000
	tCCnorm = 22
	nGen = 20

	count0 = true
	count2 = true
	count3 = true
	count1 = true

	usePT = true

	timeRec = 0.0

	kme = 8e-6

	kme01 = 9.0 * kme
	kme12 = 6.0 * kme
	gme01 = kme01 / 20.0
	gme12 = kme12 / 20.0
	gme23 = kme / 20.0
	gdem = fmin * pdem
}

func Fig3d() {

	pdem = 4e-3
	pex = 1e-3
	fmax = 4e-3

	alpha = 3.5
	beta = 2.1

	/*
		alpha = 2
		beta = 4
	*/

	defineParams()

	n3, a := runManySims()
	outputResTot(n3, a)
	fmt.Println(alpha, kme, beta, a[0], n3.me3[20], "pdem", pdem, "pex", pex, "fmax", fmax, "delGen", DelGen, "nGen", nGen)
}

func FigS10d() {
	//Change to alpharun in runManySims

	delGens := []int{1, 2, 3, 4, 5, 20}
	pdem = 4e-3
	pex = 1e-3
	fmax = 4e-3

	var a []float64
	var n3 averages

	alpha = 3.5
	beta = 2.1

	/*
		alpha = 5.966389493
		beta = 1.58113883
	*/

	printInt(delGens)

	for i := 0; i < len(delGens); i++ {

		defineParams()
		DelGen = delGens[i]

		n3, a = runManySims()
		outputResTot(n3, a)

	}

	fmt.Println("alpha", alpha, "kme", kme, "beta", beta, a[0], n3.me3[20], "pdem", pdem, "pex", pex, "fmax", fmax, "delGen", DelGen, "nGen", nGen, nSim)
}

func FigS7d() {
	deltCC := []float64{10, 11, 22, 44, 55, 88, 110}
	delGens := []int{44, 40, 20, 10, 8, 5, 4}
	pdem = 4e-3
	pex = 1e-3
	fmax = 4e-3

	var a []float64
	var n3 averages

	alpha = 3.5
	beta = 2.1

	/*
		alpha = 5.966389493
		beta = 1.58113883
	*/

	printFloat64(deltCC)

	for i := 0; i < len(deltCC); i++ {

		tCCvar = deltCC[i]
		DelGen = delGens[i]
		defineParams()

		n3, a = runManySims()
		outputResTot(n3, a)

	}

	fmt.Println("alpha", alpha, "kme", kme, "beta", beta, a[0], n3.me3[20], "pdem", pdem, "pex", pex, "fmax", fmax, "delGen", DelGen, "nGen", nGen, nSim)
}

func runManySims() (averages, []float64) {

	var n3 nucleosomes
	maxJ := (nGen*2)*int(tCC) + DelGen*int(tCCvar)

	step := 1.0

	tot := initResTot(nSim)

	for i := 0; i < nSim; i++ {

		n3 = nucleosomes{}
		n3 = initSystem(N, 3)
		n3.nCellCycle = 0

		n3 = BetaRun(alpha, beta, fmax, pdem)
		//n3 = AlphaRun(alpha, beta, fmax, pdem, 100)

		tot.me3[i] = make([]float64, maxJ+1)
		tot.me0[i] = make([]float64, maxJ+1)
		tot.trans[i] = make([]float64, (nGen*2 + DelGen))

		tot.trans[i] = n3.trans

		n3 = WriteTimeAll(n3, step, maxJ)

		tot.me3[i] = n3.me3
		tot.me0[i] = n3.me0

	}
	a := averageResTot(tot, nSim, maxJ+1)

	a = averageResTrans(tot, nSim, (nGen*2 + DelGen), a)

	return a, n3.time
}
