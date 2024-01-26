package util

import "testing"

func Test_ReverseSlice(t *testing.T) {

	s := []float32{1, 2, 3}
	rev := ReverseSlice(s)
	t.Logf("rev:%v", rev)

	s = []float32{1.5, 0.6, 3, 4, 6, 8, 6, 5, 9, 9, 9, 18, 7, 9, 12}
	rev = ReverseSlice(s)
	t.Logf("rev:%v", rev)

	s = []float32{1.5, 0.6, 3, 4, 6, 8, 6, 5, 9, 9, 9, 18, 7, 9}
	rev = ReverseSlice(s)
	t.Logf("rev:%v", rev)

	s = []float32{6, 5, 9, 9, 9, 18, 7, 9, 12, 19, 25, 45, 59, 31, 4, 14, 20, 12, 11, 10, 11}
	rev = ReverseSlice(s)
	t.Logf("rev:%v", rev)

	s = []float32{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
	rev = ReverseSlice(s)
	t.Logf("rev:%v", rev)

	s = []float32{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}
	rev = ReverseSlice(s)
	t.Logf("rev:%v", rev)
}

func Test_FormatSlice(t *testing.T) {
	s1 := []int{1, 2, 3}
	res := FormatSlice(s1)
	t.Logf("res:%v", res)

	s2 := []uint8{1, 2, 3}
	res = FormatSlice(s2)
	t.Logf("res:%v", res)

	s3 := []float64{1.22, 2, 3}
	res = FormatSlice(s3)
	t.Logf("res:%v", res)

}
