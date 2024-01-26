package cpgpos

import (
	"testing"
)

func Test_FindOverlappingCpg(t *testing.T) {
	chrcglist := []int{20, 30, 40, 50, 60}
	cpg := FindOverlappingCpg(chrcglist, 0, 10)
	t.Logf("cpg:%v", cpg)

	cpg = FindOverlappingCpg(chrcglist, 10, 35)
	t.Logf("cpg:%v", cpg)

	cpg = FindOverlappingCpg(chrcglist, 10, 40)
	t.Logf("cpg:%v", cpg)

	cpg = FindOverlappingCpg(chrcglist, 10, 60)
	t.Logf("cpg:%v", cpg)

	cpg = FindOverlappingCpg(chrcglist, 10, 70)
	t.Logf("cpg:%v", cpg)

	cpg = FindOverlappingCpg(chrcglist, 25, 35)
	t.Logf("cpg:%v", cpg)

	cpg = FindOverlappingCpg(chrcglist, 25, 50)
	t.Logf("cpg:%v", cpg)

	cpg = FindOverlappingCpg(chrcglist, 25, 60)
	t.Logf("cpg:%v", cpg)

	cpg = FindOverlappingCpg(chrcglist, 25, 70)
	t.Logf("cpg:%v", cpg)

	cpg = FindOverlappingCpg(chrcglist, 60, 70)
	t.Logf("cpg:%v", cpg)

	cpg = FindOverlappingCpg(chrcglist, 70, 90)
	t.Logf("cpg:%v", cpg)
}

func Test_whichCpgMatched(t *testing.T) {

	overlappingCpg := []int{100, 200, 300, 400, 500}

	matchedCpgs, nextIdx := whichCpgMatched(overlappingCpg, 0, 5, 10)
	t.Logf("matched cpgs:%v, nextIdx:%d", matchedCpgs, nextIdx)

	matchedCpgs, nextIdx = whichCpgMatched(overlappingCpg, 1, 5, 10)
	t.Logf("matched cpgs:%v, nextIdx:%d", matchedCpgs, nextIdx)

	matchedCpgs, nextIdx = whichCpgMatched(overlappingCpg, 0, 100, 200)
	t.Logf("matched cpgs:%v, nextIdx:%d", matchedCpgs, nextIdx)

	matchedCpgs, nextIdx = whichCpgMatched(overlappingCpg, 2, 100, 200)
	t.Logf("matched cpgs:%v, nextIdx:%d", matchedCpgs, nextIdx)

	matchedCpgs, nextIdx = whichCpgMatched(overlappingCpg, 0, 500, 600)
	t.Logf("matched cpgs:%v, nextIdx:%d", matchedCpgs, nextIdx)

	matchedCpgs, nextIdx = whichCpgMatched(overlappingCpg, 0, 600, 700)
	t.Logf("matched cpgs:%v, nextIdx:%d", matchedCpgs, nextIdx)
}
