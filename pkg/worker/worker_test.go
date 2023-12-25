package worker

import (
	"testing"
)

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

func Test_getkineticswin(t *testing.T) {
	nums := []uint8{1, 2, 3}
	kwin, err := getkineticswin(UInt8Slice(nums), true)
	if err != nil {
		t.Errorf("err:%+v", err)
		return
	}
	t.Logf("kwin:%+v", kwin)
}

func Test_findOverlappingCpg(t *testing.T) {
	chrcglist := []int{20, 30, 40, 50, 60}
	cpg := findOverlappingCpg(chrcglist, 0, 10)
	t.Logf("cpg:%v", cpg)

	cpg = findOverlappingCpg(chrcglist, 10, 35)
	t.Logf("cpg:%v", cpg)

	cpg = findOverlappingCpg(chrcglist, 10, 40)
	t.Logf("cpg:%v", cpg)

	cpg = findOverlappingCpg(chrcglist, 10, 60)
	t.Logf("cpg:%v", cpg)

	cpg = findOverlappingCpg(chrcglist, 10, 70)
	t.Logf("cpg:%v", cpg)

	cpg = findOverlappingCpg(chrcglist, 25, 35)
	t.Logf("cpg:%v", cpg)

	cpg = findOverlappingCpg(chrcglist, 25, 50)
	t.Logf("cpg:%v", cpg)

	cpg = findOverlappingCpg(chrcglist, 25, 60)
	t.Logf("cpg:%v", cpg)

	cpg = findOverlappingCpg(chrcglist, 25, 70)
	t.Logf("cpg:%v", cpg)

	cpg = findOverlappingCpg(chrcglist, 60, 70)
	t.Logf("cpg:%v", cpg)

	cpg = findOverlappingCpg(chrcglist, 70, 90)
	t.Logf("cpg:%v", cpg)
}

func Test_reverseSlice(t *testing.T) {

	s := []uint8{1, 2, 3}
	rev := reverseSlice(UInt8Slice(s))
	t.Logf("rev:%v", rev)

}

func Test_joinSlice(t *testing.T) {
	var s interface{}
	s = []int{1, 2, 3}
	res := joinSlice(s)
	t.Logf("res:%v", res)

	s = []uint8{1, 2, 3}
	res = joinSlice(s)
	t.Logf("res:%v", res)

	s = []float64{1.22, 2, 3}
	res = joinSlice(s)
	t.Logf("res:%v", res)

}
