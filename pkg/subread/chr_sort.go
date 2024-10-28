package subread

import (
	"fmt"
	"sort"
	"strconv"
	"strings"
)

type ComparableChrCpg struct {
	Chr string
	Cpg int
}

func sortChrCpgs(chrCpgs []string, refOrderMap map[string]int) ([]string, error) {

	var res []string

	var comChrCpgs []ComparableChrCpg

	//先转成map
	for _, chrCpg := range chrCpgs {
		parts := strings.Split(chrCpg, "->")
		if len(parts) != 2 {
			return nil, fmt.Errorf("invalid chrCpg:%s", chrCpg)
		}
		if _, exist := refOrderMap[parts[0]]; !exist {
			return nil, fmt.Errorf("unknown ref chr:%s", chrCpg)
		}
		cpg, err := strconv.Atoi(parts[1])
		if err != nil {
			return nil, err
		}
		comChrCpg := ComparableChrCpg{
			Chr: parts[0],
			Cpg: cpg,
		}
		comChrCpgs = append(comChrCpgs, comChrCpg)
	}

	sort.Slice(comChrCpgs, func(i, j int) bool {

		chrI := comChrCpgs[i].Chr
		chrJ := comChrCpgs[j].Chr
		if refOrderMap[chrI] < refOrderMap[chrJ] {
			return true
		}
		if refOrderMap[chrI] > refOrderMap[chrJ] {
			return false
		}
		if refOrderMap[chrI] == refOrderMap[chrJ] {
			return comChrCpgs[i].Cpg < comChrCpgs[j].Cpg
		}
		return false
	})

	for _, comChrCpg := range comChrCpgs {
		chrCpg := fmt.Sprintf("%s->%d", comChrCpg.Chr, comChrCpg.Cpg)
		res = append(res, chrCpg)
	}
	return res, nil
}

func sortCpgOutput(cpgOutputs []*CpgOutput, refOrderMap map[string]int) {
	sort.Slice(cpgOutputs, func(i, j int) bool {

		refI := cpgOutputs[i].Ref
		refJ := cpgOutputs[j].Ref

		if refOrderMap[refI] < refOrderMap[refJ] {
			return true
		}

		if refOrderMap[refI] > refOrderMap[refJ] {
			return false
		}
		if refOrderMap[refI] == refOrderMap[refJ] {
			return cpgOutputs[i].Cpg < cpgOutputs[j].Cpg
		}
		return false
	})
}
