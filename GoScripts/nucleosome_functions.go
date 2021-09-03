package main

func Kdelta(x int, y int) float64 {
	switch x {
	case y:
		return 1
	}
	return 0
}

func Eihelp(item int) float64 {
	return pme2*(Kdelta(item, 2)) + (Kdelta(item, 3))
}

func Ei(n []nucleosome, ind int) []float64 {

	var neighbors []int
	Ein := make([]float64, 2)
	sum := float64(0)
	L := len(n)

	switch ind {

	case (L - 1):
		neighbors = append(neighbors, ind-1)

	case 0:
		neighbors = append(neighbors, ind+1)
	default:
		neighbors = append(neighbors, ind-1, ind+1)
	}

	for i := 0; i < len(neighbors); i++ {
		sum = sum + Eihelp(n[neighbors[i]].h1) + Eihelp(n[neighbors[i]].h2)
	}

	Ein[0] = sum + Eihelp(n[ind].h2)
	Ein[1] = sum + Eihelp(n[ind].h1)

	return Ein
}

func (n nucleosomes) Props() []reaction {

	nuc := n.values

	var ind, ind2 int
	var currNuc nucleosome
	numR := 1 + 4*len(n.values)
	d := make([]reaction, numR)
	currRate := 0.0

	for i := 0; i < len(nuc); i++ {
		ind = i
		currEi := Ei(nuc, ind)

		currNuc = nuc[ind]
		nuc[ind].me1 = helpRm(currNuc.h1, currEi[0])
		nuc[ind].me2 = helpRm(currNuc.h2, currEi[1])
		nuc[ind].de1 = helpDm(currNuc.h1)
		nuc[ind].de2 = helpDm(currNuc.h2)

		ind2 = (i * 4)

		currRate = nuc[ind].me1 + currRate
		d[ind2].rate = currRate
		d[ind2].action = "meth"
		d[ind2].nucleosome = i
		d[ind2].histone = "h1"

		currRate = currRate + nuc[ind].me2
		d[ind2+1].rate = currRate
		d[ind2+1].action = "meth"
		d[ind2+1].nucleosome = i
		d[ind2+1].histone = "h2"

		currRate = currRate + nuc[ind].de1
		d[ind2+2].rate = currRate
		d[ind2+2].action = "demeth"
		d[ind2+2].nucleosome = i
		d[ind2+2].histone = "h1"

		currRate = currRate + nuc[ind].de2
		d[ind2+3].rate = currRate
		d[ind2+3].action = "demeth"
		d[ind2+3].nucleosome = i
		d[ind2+3].histone = "h2"

	}

	n.transcription = n.Rtf()
	currRate = currRate + n.transcription
	d[numR-1].rate = currRate
	d[numR-1].action = "transcription"

	return d
}

func helpRm(h int, currEi float64) float64 {

	return beta * (Kdelta(h, 0)*(gme01+kme01*currEi) + Kdelta(h, 1)*(gme12+kme12*currEi) + Kdelta(h, 2)*(gme23+kme*currEi))
}

func helpDm(h int) float64 {

	return gdem * (Kdelta(h, 1) + Kdelta(h, 2) + Kdelta(h, 3))
}

func (n nucleosomes) Rtf() float64 {

	nuc := n.values
	L := len(nuc)
	sum := float64(0)
	var res float64

	for i := 0; i < L; i++ {
		sum = Kdelta(n.values[i].h1, 2) + Kdelta(n.values[i].h2, 2) + Kdelta(n.values[i].h1, 3) + Kdelta(n.values[i].h2, 3) + sum

	}
	sum = sum / float64(L*2)

	if usePT {
		switch {
		case sum < PT:
			res = alpha * (fmax - (sum/PT)*(fmax-fmin))
		case sum > PT:
			res = alpha * fmin
		case sum == PT:
			res = alpha * fmin
		}
	} else {
		res = alpha * (fmax - sum*(fmax-fmin))
	}

	switch {
	case res > transMax:
		res = transMax

	case res == transMax:
		res = transMax
	}

	return res

}

func (n nucleosomes) DeM() {

	for i := 0; i < len(n.values); i++ {
		switch {
		case rng.Float64() < pdem:
			n.values[i].h1 = deMeth(n.values[i].h1)
		}

		switch {
		case rng.Float64() < pdem:
			n.values[i].h2 = deMeth(n.values[i].h2)
		}

	}

}

func (n nucleosomes) Ex() {

	nuc := n.values

	for i := 0; i < len(nuc); i++ {
		switch {
		case rng.Float64() < pex:
			nuc[i].h1 = 0
			nuc[i].h2 = 0

		}

		switch {
		case rng.Float64() < pex:
			nuc[i].h1 = 0
			nuc[i].h2 = 0

		}

	}
}

func (n nucleosomes) Rep() {

	nuc := n.values

	for i := 0; i < len(nuc); i++ {
		switch {
		case rng.Float64() < 0.5:
			nuc[i].h1 = 0
			nuc[i].h2 = 0

		}

	}

}

func r2n(n nucleosomes, r reaction) {
	currAction := r.action

	switch currAction {
	case "meth":
		if r.histone == "h1" {
			n.values[r.nucleosome].h1 = addMeth(n.values[r.nucleosome].h1)
		} else {
			n.values[r.nucleosome].h2 = addMeth(n.values[r.nucleosome].h2)
		}

	case "demeth":
		if r.histone == "h1" {
			n.values[r.nucleosome].h1 = deMeth(n.values[r.nucleosome].h1)
		} else {
			n.values[r.nucleosome].h2 = deMeth(n.values[r.nucleosome].h2)
		}

	case "transcription":
		n.DeM()
		n.Ex()
		n.trans[n.nCellCycle]++
	}

}
