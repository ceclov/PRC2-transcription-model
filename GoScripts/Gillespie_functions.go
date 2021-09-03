package main

import (
	"math"
)

func DirectGillespie(n nucleosomes) float64 {
	var thisTime float64

	d := n.Props()

	r := rng.Float64()

	a0 := d[len(d)-1].rate

	r2 := rng.Float64()

	thisTime = (1 / a0) * math.Log(1/r)

	dr := FindReac(d, r2*a0)

	r2n(n, dr)

	return thisTime
}

func FindReac(d []reaction, item float64) reaction {
	var dr reaction
	found := false

	for i := 0; i < len(d); i++ {
		if d[i].rate > item || d[i].rate == item {
			dr = d[i]
			found = true
			break
		}
	}
	if !found {
		panic("no reaction found")
	}

	return dr
}

func GStep(n nucleosomes, tmax float64, nGen int, currTime float64, genStart int) (nucleosomes, float64) {

	var currCCtime float64
	currCCtime = currTime
	genTime := 0.0
	Btime := 0.0

	var currRate float64

	for i := genStart; i < nGen; i++ {

		for {

			currRate = DirectGillespie(n)

			genTime = genTime + currRate/3600.0

			if genTime > tmax {
				n.Pme3end[i] = meanState(n, 3)

				n.nCellCycle++
				if i < (nGen - 1) {

					n.Rep()
				}

				genTime = 0.0
				currCCtime = currCCtime + tmax

				currTime = currCCtime
				break
			}

			currTime = currTime + (currRate / 3600.0)

			if genTime > timeRec {

				n.time = append(n.time, currTime)

				Btime = genTime - timeRec + float64(i)

				n.Btime = append(n.Btime, Btime)

				n.cellcycle = append(n.cellcycle, i)

				if count0 {

					n.me0 = append(n.me0, meanState(n, 0))

				}

				if count3 {

					n.me3 = append(n.me3, meanState(n, 3))

				}

				if count2 {

					n.me2 = append(n.me2, meanState(n, 2))

				}

				if count1 {

					n.me1 = append(n.me1, meanState(n, 1))

				}

			}

		}

	}

	return n, currTime
}
