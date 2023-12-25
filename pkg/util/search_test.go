package util

import "testing"

func TestFindUpBoundIndex(t *testing.T) {
	nums := []int{1, 3, 5, 7, 9, 11, 13}
	flag, index := FindUpBoundIndex(nums, 15)
	t.Logf("flag:%v, index:%v", flag, index)
	flag, index = FindUpBoundIndex(nums, -1)
	t.Logf("flag:%v, index:%v", flag, index)

	flag, index = FindUpBoundIndex(nums, 6)
	t.Logf("flag:%v, index:%v", flag, index)
	flag, index = FindUpBoundIndex(nums, 7)
	t.Logf("flag:%v, index:%v", flag, index)

	nums = []int{1, 3, 5, 7, 9, 11, 13, 15}
	flag, index = FindUpBoundIndex(nums, 17)
	t.Logf("flag:%v, index:%v", flag, index)
	flag, index = FindUpBoundIndex(nums, -1)
	t.Logf("flag:%v, index:%v", flag, index)

	flag, index = FindUpBoundIndex(nums, 6)
	t.Logf("flag:%v, index:%v", flag, index)
	flag, index = FindUpBoundIndex(nums, 7)
	t.Logf("flag:%v, index:%v", flag, index)
}

func TestFindLowBoundIndex(t *testing.T) {
	nums := []int{1, 3, 5, 7, 9, 11, 13}
	flag, index := FindLowBoundIndex(nums, 15)
	t.Logf("flag:%v, index:%v", flag, index)
	flag, index = FindLowBoundIndex(nums, -1)
	t.Logf("flag:%v, index:%v", flag, index)

	flag, index = FindLowBoundIndex(nums, 6)
	t.Logf("flag:%v, index:%v", flag, index)
	flag, index = FindLowBoundIndex(nums, 7)
	t.Logf("flag:%v, index:%v", flag, index)

	nums = []int{1, 3, 5, 7, 9, 11, 13, 15}
	flag, index = FindLowBoundIndex(nums, 17)
	t.Logf("flag:%v, index:%v", flag, index)
	flag, index = FindLowBoundIndex(nums, -1)
	t.Logf("flag:%v, index:%v", flag, index)

	flag, index = FindLowBoundIndex(nums, 6)
	t.Logf("flag:%v, index:%v", flag, index)
	flag, index = FindLowBoundIndex(nums, 7)
	t.Logf("flag:%v, index:%v", flag, index)
}
