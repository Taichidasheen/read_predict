package worker

import (
	"fmt"
	"strings"
)

/*func reverseSlice[T any](arr []T) []T {
	var res []T
	length := len(arr)
	for i := length - 1; i >= 0; i-- {
		res = append(res, arr[i])
	}
	return res
}*/

func reverseSliceGeneric[T any](arr []T) []T {
	res := make([]T, len(arr))
	length := len(arr)
	for i := length - 1; i >= 0; i-- {
		res[length-i-1] = arr[i]
	}
	return res
}

func reverseSlice(arr []float32) []float32 {

	res := make([]float32, len(arr))
	length := len(arr)
	for i := length - 1; i >= 0; i-- {
		res[length-i-1] = arr[i]
	}
	return res
}

func reverseSliceByte(arr []byte) []byte {

	res := make([]byte, len(arr))
	length := len(arr)

	for i := length - 1; i >= 0; i-- {
		res[length-i-1] = arr[i]
	}
	return res
}

// rbegin和rend表示对于反转后的arr数组，返回从rbegin到rend(不包括rend)区间内的元素
func reverseSliceByteByIndex(arr []byte, rbegin, rend int) []byte {

	scope := rend - rbegin
	res := make([]byte, scope)

	arrLen := len(arr)

	//根据rbegin和rend推算出目标元素在反转数组前的区间范围
	fbegin := arrLen - rend
	fend := arrLen - rbegin //不包括fend

	pos := 0
	for i := fend - 1; i >= fbegin; i-- {
		res[pos] = arr[i]
		pos = pos + 1
	}
	return res
}

func formatSlice[T any](arr []T) string {
	strs := make([]string, len(arr))
	for i, num := range arr {
		strs[i] = fmt.Sprintf("%v", num)
	}
	return strings.Join(strs, ",")
}
