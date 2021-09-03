package main

type nucleosome struct {
	h1  int `validate:"nonzero"`
	h2  int `validate:"number,min=0,max=3"`
	me1 float64
	me2 float64
	de1 float64
	de2 float64
	lab string
}

type reaction struct {
	rate       float64
	nucleosome int
	histone    string
	action     string
}

type nucleosomes struct {
	values        []nucleosome
	transcription float64
	time          []float64
	Btime         []float64
	me3           []float64

	me2 []float64

	me1 []float64

	me0 []float64

	cellcycle  []int
	Pme3end    []float64
	Pme1end    float64
	Pme2end    float64
	nCellCycle int
	trans      []float64
}

type rates interface {
	Rates()
	Rtf() float64
	Props() []reaction
}

type conversion interface {
	n2r()
}

type transcription interface {
	DeM()
	Ex()
}

type Vector []float64

type Matrix [][]float64

type results struct {
	me3   [][]float64
	me0   [][]float64
	trans [][]float64
}

type averages struct {
	me3   []float64
	me1   []float64
	me0   []float64
	trans []float64
}
