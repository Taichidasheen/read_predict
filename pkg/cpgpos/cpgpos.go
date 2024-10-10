package cpgpos

import (
	"github.com/Taichidasheen/read_predict/pkg/util"
	"github.com/biogo/hts/sam"
)

func FindOverlappingCpg(chrcglist []int, refStart, refEnd int) []int {
	var overlappingCpg []int
	//查找比refStart大的第一个元素
	_, leftIndex := util.FindUpBoundIndex(chrcglist, refStart)
	if leftIndex > len(chrcglist) {
		return nil
	}
	//查找比refEnd小的第一个元素
	_, rightIndex := util.FindLowBoundIndex(chrcglist, refEnd)
	if rightIndex < 0 {
		return nil
	}
	if leftIndex <= rightIndex {
		for i := leftIndex; i <= rightIndex; i++ {
			overlappingCpg = append(overlappingCpg, chrcglist[i])
		}
	}
	return overlappingCpg

}

func LocateCpgPosOnSeq(alnRefStart int, readCigar sam.Cigar, overlappingCpg []int) ([]int, map[int]int) {
	var locatedCpgs []int
	cpgPosOnSeq := make(map[int]int)
	cpgBeginIdx := 0
	seqLeadingFlag := 1
	seqPosBlkStart := 0
	seqPosBlkEnd := 0
	refWalkingPosStart := 0
	refWalkingPosEnd := 0
	for _, cigar := range readCigar {
		op := cigar.Type()
		count := cigar.Len()
		if (op == sam.CigarInsertion || op == sam.CigarSoftClipped || op == sam.CigarHardClipped ||
			op == sam.CigarPadded) && seqLeadingFlag == 1 {
			//1,4 consume reads, 5_hard_clip should also be counted on reads, since the aligner doesn't change ipd/pw signal vector as the SEQ
			if op != sam.CigarPadded {
				seqPosBlkStart = seqPosBlkEnd
				seqPosBlkEnd = seqPosBlkStart + count
			}
		} else {
			//注意：下面这段逻辑有点奇怪，可以简化
			if seqLeadingFlag == 1 {
				seqLeadingFlag = 2
				if (op == sam.CigarMatch || op == sam.CigarDeletion || op == sam.CigarSkipped ||
					op == sam.CigarEqual || op == sam.CigarMismatch) && seqLeadingFlag == 2 {
					refWalkingPosStart = alnRefStart
					refWalkingPosEnd = refWalkingPosStart + count
				}
				if (op == sam.CigarMatch || op == sam.CigarInsertion || op == sam.CigarSoftClipped ||
					op == sam.CigarEqual || op == sam.CigarMismatch) && seqLeadingFlag == 2 {
					seqPosBlkStart = seqPosBlkEnd
					seqPosBlkEnd = seqPosBlkStart + count
				}
			} else {
				if op == sam.CigarMatch || op == sam.CigarDeletion || op == sam.CigarSkipped ||
					op == sam.CigarEqual || op == sam.CigarMismatch {
					refWalkingPosStart = refWalkingPosEnd
					refWalkingPosEnd = refWalkingPosStart + count
				}
				if op == sam.CigarMatch || op == sam.CigarInsertion || op == sam.CigarSoftClipped ||
					op == sam.CigarEqual || op == sam.CigarMismatch {
					seqPosBlkStart = seqPosBlkEnd
					seqPosBlkEnd = seqPosBlkStart + count
				}
			}
			if op != sam.CigarDeletion && op != sam.CigarSkipped {
				//检查当前refStart和refEnd是否包含某个cpg
				matchedCpgs, nextIdx := whichCpgMatched(overlappingCpg, cpgBeginIdx, refWalkingPosStart, refWalkingPosEnd)
				//if alnRefStart == 49196398 {
				//	log.Printf("cpgBeginIdx:%d, refWalkingPosStart:%d, refWalkingPosEnd:%d", cpgBeginIdx, refWalkingPosStart, refWalkingPosEnd)
				//	log.Printf("mathcedCpgs:%v, nextIdx:%d", matchedCpgs, nextIdx)
				//}

				//if alnRefStart == 10770114 {
				//	log.Printf("cpgBeginIdx:%d, refWalkingPosStart:%d, refWalkingPosEnd:%d", cpgBeginIdx, refWalkingPosStart, refWalkingPosEnd)
				//	log.Printf("mathcedCpgs:%v, nextIdx:%d", matchedCpgs, nextIdx)
				//}

				cpgBeginIdx = nextIdx
				for _, cpg := range matchedCpgs {
					lastOpNeeded := cpg - refWalkingPosStart
					if op == sam.CigarMatch || op == sam.CigarInsertion || op == sam.CigarSoftClipped ||
						op == sam.CigarEqual || op == sam.CigarMismatch {
						seqPos := seqPosBlkStart + lastOpNeeded
						cpgPosOnSeq[cpg] = seqPos - 1
						locatedCpgs = append(locatedCpgs, cpg)
					}
				}
				if nextIdx >= len(overlappingCpg) {
					//对比结束
					return locatedCpgs, cpgPosOnSeq
				}

				/*if nextIdx > len(overlappingCpg) {
					//对比结束
					return locatedCpgs, cpgPosOnSeq
				} else {
					cpgBeginIdx = nextIdx
					for _, cpg := range matchedCpgs {
						lastOpNeeded := cpg - refWalkingPosStart
						if op == sam.CigarMatch || op == sam.CigarInsertion || op == sam.CigarSoftClipped ||
							op == sam.CigarEqual || op == sam.CigarMismatch {
							seqPos := seqPosBlkStart + lastOpNeeded
							cpgPosOnSeq[cpg] = seqPos - 1
							locatedCpgs = append(locatedCpgs, cpg)
						}
					}
				}*/
			}
		}
	}
	return locatedCpgs, cpgPosOnSeq
}

// 判断哪些cpg在当前给定的范围内，返回这些cpg的值和在cpglist中的坐标
// cpgBeginIdx表示从overlappingCpg的哪个位置开始比对（cpgBeginIdx之前的位置可以不用比对）
func whichCpgMatched(overlappingCpg []int, cpgBeginIdx int, refStart, refEnd int) ([]int, int) {
	var matchedCpgs []int
	var cpgNextIndex int //下一次查找overlappingCpg数组的开始位置

	for i := cpgBeginIdx; i < len(overlappingCpg); i++ {
		if refStart <= overlappingCpg[i] && overlappingCpg[i] <= refEnd {
			matchedCpgs = append(matchedCpgs, overlappingCpg[i])
			cpgNextIndex = i + 1
		}
		if refStart > overlappingCpg[i] {
			cpgNextIndex = i + 1
		}
		if refEnd < overlappingCpg[i] {
			//不用再比对了，等到refstart和refend更新后再从i位置开始比对
			cpgNextIndex = i
			break
		}
	}
	return matchedCpgs, cpgNextIndex
}
