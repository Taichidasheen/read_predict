package worker

import (
	"testing"
)

func Test_getkineticswin(t *testing.T) {
	nums := []uint8{1, 2, 3}
	kwin, err := getkineticswin(nums, true)
	if err != nil {
		t.Errorf("err:%+v", err)
		return
	}
	t.Logf("kwin:%+v", kwin)

	nums = []uint8{11, 12, 13, 14}
	kwin, err = getkineticswin(nums, true)
	if err != nil {
		t.Errorf("err:%+v", err)
		return
	}
	t.Logf("kwin:%+v", kwin)

	nums = []uint8{1, 1, 1}
	kwin, err = getkineticswin(nums, true)
	if err != nil {
		t.Errorf("err:%+v", err)
		return
	}
	t.Logf("kwin:%+v", kwin)

}

func Test_getkineticswinQuick(t *testing.T) {
	nums := []uint8{1, 2, 3}
	kwin, err := getkineticswinQuick(nums, true)
	if err != nil {
		t.Errorf("err:%+v", err)
		return
	}
	t.Logf("kwin:%+v", kwin)

	nums = []uint8{11, 12, 13, 14}
	kwin, err = getkineticswinQuick(nums, true)
	if err != nil {
		t.Errorf("err:%+v", err)
		return
	}
	t.Logf("kwin:%+v", kwin)

	nums = []uint8{1, 1, 1}
	kwin, err = getkineticswinQuick(nums, true)
	if err != nil {
		t.Errorf("err:%+v", err)
		return
	}
	t.Logf("kwin:%+v", kwin)

}
