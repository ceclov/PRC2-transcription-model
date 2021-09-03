package main

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"math"
	"os"
	"strconv"
	"strings"
)

func initSystem(N int, state int) nucleosomes {

	var n nucleosomes
	n.values = make([]nucleosome, N)

	switch {

	case (N < 0):
		panic("Length of system to small")
	case N == 0:
		panic("System size is zero")

	case state < 0:
		panic("Init state to low")
	case state > 3:
		panic("init state to high")
	}

	for i := 0; i < N; i++ {
		n.values[i].h1 = state
		n.values[i].h2 = state

	}

	n.transcription = 0.0
	if count0 {
		n.me0 = []float64{}

	}

	if count3 {
		n.me3 = []float64{}

	}

	if count2 {
		n.me2 = []float64{}

	}

	if count1 {

		n.me1 = []float64{}

	}

	n.time = []float64{}

	n.Btime = []float64{}

	return n

}

func paramsFour(a, b, c, d []float64) ([]float64, []float64, []float64, []float64) {
	L := len(a) * len(b) * len(c) * len(d)
	resA := make([]float64, L)
	resB := make([]float64, L)
	resC := make([]float64, L)
	resD := make([]float64, L)

	ind := 0
	for i := 0; i < len(a); i++ {

		for j := 0; j < len(b); j++ {

			for u := 0; u < len(c); u++ {
				for k := 0; k < len(d); k++ {
					resA[ind] = a[i]
					resB[ind] = b[j]
					resC[ind] = c[u]
					resD[ind] = d[k]

					ind++
				}
			}
		}
	}
	return resA, resB, resC, resD
}

func divideParams(a, b, c, d []float64, r float64) {
	l := len(a)
	div := math.Trunc(float64(l) / r)

	var t string
	ind := 0

	currMap := make(map[string][]float64)

	for i := 0; i < int(r); i++ {
		t = "Params_div" + strconv.Itoa(i) + ".json"
		currMap = map[string][]float64{}
		for j := 0; j < int(div); j++ {

			currMap["alpha"] = append(currMap["alpha"], a[ind])
			currMap["beta"] = append(currMap["beta"], b[ind])
			currMap["fmaxs"] = append(currMap["fmaxs"], c[ind])
			currMap["pdems"] = append(currMap["pdems"], d[ind])
			ind++
		}
		maptoFile(currMap, t)

	}

	for {

		if ind > (l - 1) {

			maptoFile(currMap, t)
			break
		}
		currMap["alpha"] = append(currMap["alpha"], a[ind])
		currMap["beta"] = append(currMap["beta"], b[ind])
		currMap["fmaxs"] = append(currMap["fmaxs"], c[ind])
		currMap["pdems"] = append(currMap["pdems"], c[ind])

		ind++
	}

}

func maptoFile(m map[string][]float64, name string) {

	jsonData, err := json.Marshal(m)

	if err != nil {
		panic("error in writing to file")
	}
	jsonFile, err := os.Create(name)

	defer jsonFile.Close()

	jsonFile.Write(jsonData)
	jsonFile.Close()
	fmt.Println("JSON data written to ", jsonFile.Name())

}

func ResetSystem(n nucleosomes) nucleosomes {

	if count0 {
		n.me0 = []float64{}

	}

	if count3 {
		n.me3 = []float64{}

	}

	if count2 {
		n.me2 = []float64{}

	}

	if count1 {

		n.me1 = []float64{}

	}

	n.time = []float64{}

	n.Btime = []float64{}
	n.nCellCycle = 0
	n.trans = []float64{}
	n.Pme3end = []float64{}
	return n

}

func sumFloat(v []float64) float64 {
	sum := 0.0
	for i := 0; i < len(v); i++ {
		sum = sum + v[i]
	}
	return sum
}

func meanState(n nucleosomes, state int) float64 {
	N := len(n.values)
	if N == 0 {

		return 0
	}
	sum := 0
	for i := 0; i < N; i++ {
		switch n.values[i].h1 {
		case state:
			sum++
		}

		switch n.values[i].h2 {
		case state:
			sum++
		}

	}

	return float64(sum) / float64(2*N)
}

func outputResTot(n averages, time []float64) {

	printFloat64(time)

	if count3 {

		printFloat64(n.me3)

	}

	if count0 {
		printFloat64(n.me0)
	}

	printFloat64(n.trans)

}

func printFloat64(v []float64) {

	for i := 0; i < len(v); i++ {
		fmt.Print(v[i], "\t")
	}
	fmt.Println()
}

func printInt(v []int) {

	for i := 0; i < len(v); i++ {
		fmt.Print(v[i], "\t")
	}
	fmt.Println()
}

func deMeth(value int) int {

	switch value {
	case 1, 2, 3:
		return value - 1

	}
	return value
}

func addMeth(value int) int {

	switch value {
	case 0, 1, 2:
		return value + 1

	}
	return value
}

func BM(n nucleosomes) float64 {

	sumOFF := 0.0
	sumON := 0.0
	var ON, OFF []float64
	var res float64

	ON = ONOFF(n.me2, n.me3, true)
	OFF = ONOFF(n.me2, n.me3, false)

	sumOFF = timeAverage(n.Btime, OFF)
	sumON = timeAverage(n.Btime, ON)

	res = 4 * sumOFF * sumON

	return res

}

func ONOFF(me2, me3 []float64, ON bool) []float64 {
	l := len(me2)
	res := make([]float64, l)
	var sum float64

	for i := 0; i < l; i++ {
		sum = me2[i] + me3[i]

		res[i] = 0

		if ON {
			if sum < 0.25 {
				res[i] = 1.0
			}
		} else {
			if sum > 0.75 {
				res[i] = 1.0
			}
		}

	}

	return res
}

func timeAverage(time, states []float64) float64 {

	var timeDiff float64
	K := len(states)
	res := 0.0
	nom := time[len(time)-1] - time[0]

	for i := 0; i < (K - 1); i++ {
		timeDiff = time[i+1] - time[i]
		switch {
		case timeDiff < 0:
			fmt.Println(timeDiff, time[i], time[i+1])
			panic("In timeAverage timeDiff is less than 0.0")

		}

		res = (states[i] * (timeDiff) / nom) + res
	}

	return res
}

func getParam(start, end float64, l int, log bool) []float64 {
	res := make([]float64, l)
	step := (end - start) / float64(l)

	if log {

		step = (math.Log(end) - math.Log(start)) / float64(l)
	}
	currRes := start
	res[0] = currRes
	for i := 1; i < l; i++ {
		if log {
			currRes = currRes * math.Exp(step)
		} else {
			currRes = currRes + step
		}
		res[i] = currRes
	}

	return res
}

func meanFloat(v []float64) float64 {

	l := len(v)
	sum := 0.0

	for i := 0; i < l; i++ {
		sum += v[i]
	}

	return sum / float64(l)
}

func importParam(name string) []float64 {

	var input []float64
	raw, err := ioutil.ReadFile(name)
	if err != nil {
		fmt.Println(err.Error())

	}
	json.Unmarshal(raw, &input)

	return input
}

func importJsonParam(name string) ([]float64, []float64, []float64, []float64) {

	var input map[string][]float64
	raw, err := ioutil.ReadFile(name)
	if err != nil {
		fmt.Println(err.Error())

	}
	json.Unmarshal(raw, &input)

	return input["alpha"], input["beta"], input["fmaxs"], input["pdems"]
}
func ImportParamCsv(name string) []float64 {
	f, _ := os.Open(name)
	var input []float64
	r := csv.NewReader(bufio.NewReader(f))

	for {
		record, err := r.Read()

		if err == io.EOF {
			break
		}
		v, err := strconv.ParseFloat(strings.Join(record, ""), 64)
		input = append(input, float64(v))

	}
	return input
}

func ReadCsv(name string) map[string][]float64 {
	f, _ := os.Open(name)
	input := make(map[string][]float64)
	r := csv.NewReader(bufio.NewReader(f))

	var names []string
	i := 0

	for {
		record, err := r.Read()

		if err == io.EOF {
			break
		}
		if i == 0 {
			names = record

			for j := 0; j < len(names); j++ {
				input[names[j]] = []float64{}
			}

		} else {
			for j := 0; j < len(names); j++ {
				v, _ := strconv.ParseFloat(record[j], 64)
				input[names[j]] = append(input[names[j]], float64(v))

			}
		}
		i++
	}
	return input
}

func WriteTimeAll(n nucleosomes, step float64, l int) nucleosomes {
	tout := 0.0
	currT := 0.0
	output := nucleosomes{}

	output.time = make([]float64, l+1)

	output.me3 = make([]float64, l+1)
	output.me1 = make([]float64, l+1)
	output.me0 = make([]float64, l+1)

	L := len(n.time)
	output.time[0] = tout

	output.me3[0] = n.me3[0]
	output.me1[0] = n.me1[0]
	output.me0[0] = n.me0[0]

	j := 0
	for i := 1; i < L; i++ {

		for {

			if tout > n.time[i] {

				break
			} else {

				currT = tout

				output.me3[j] = n.me3[i]
				output.me1[j] = n.me1[i]
				output.me0[j] = n.me0[i]

				output.time[j] = currT

				tout += step

				j++
			}

		}
	}

	maxJ := float64(l)
	if !(currT > maxJ || currT == maxJ) {

		for {

			if tout > maxJ {
				break
			}
			currT = tout
			output.time[j] = currT

			output.me3[j] = n.me3[L-1]
			output.me1[j] = n.me1[L-1]
			output.me0[j] = n.me0[L-1]

			j++
			tout += step
		}
	}

	return output
}

func initResTot(m int) results {

	var r results

	r.me3 = make([][]float64, m)
	r.me0 = make([][]float64, m)
	r.trans = make([][]float64, m)

	return r

}

func getColumn(mat [][]float64, m, n int) []float64 {
	res := make([]float64, m)

	for k := 0; k < m; k++ {

		res[k] = mat[k][n]

	}

	return res
}

func averageResTot(r results, m, n int) averages {

	var a averages

	a.me3 = make([]float64, n)
	a.me0 = make([]float64, n)

	for i := 0; i < n; i++ {

		a.me3[i] = meanFloat(getColumn(r.me3, m, i))
		a.me0[i] = meanFloat(getColumn(r.me0, m, i))

	}
	return a
}

func averageResTrans(r results, m, n int, a averages) averages {
	a.trans = make([]float64, n)

	for i := 0; i < n; i++ {

		a.trans[i] = meanFloat(getColumn(r.trans, m, i))

	}
	return a
}
