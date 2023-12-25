package util

// 在升序数组中查找target的index
// 如果找到返回true和target在数组中index
// 如果找不到返回false和数组中第一个比target大的值的index
// left有可能返回len(nums)
func FindUpBoundIndex(nums []int, target int) (bool, int) {
	left := 0
	right := len(nums) - 1
	for left <= right {
		mid := (left + right) / 2
		if target > nums[mid] {
			left = mid + 1 //注意这里有可能导致nums[left]比target大
		} else if target < nums[mid] {
			right = mid - 1 //注意这里有可能导致nums[right]比target小
		} else {
			return true, mid
		}
	}
	//for循环结束时，right < left
	return false, left
}

// 在升序数组中查找target的index
// 如果找到返回true和target在数组中index
// 如果找不到返回false和在数组中第一个比target小的值的index
// right 有可能返回 -1
func FindLowBoundIndex(nums []int, target int) (bool, int) {
	left := 0
	right := len(nums) - 1
	for left <= right {
		mid := (left + right) / 2
		if target > nums[mid] {
			left = mid + 1 //注意这里有可能导致nums[left]比target大
		} else if target < nums[mid] {
			right = mid - 1 //注意这里有可能导致nums[right]比target小
		} else {
			return true, mid
		}
	}
	//for循环结束时，right < left
	return false, right
}
