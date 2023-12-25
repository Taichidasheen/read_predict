package worker

import (
	"fmt"
	"github.com/Taichidasheen/subreads_locate/pkg/util"
	"github.com/biogo/hts/sam"
	"github.com/montanaflynn/stats"
	"log"
	"sort"
	"strings"
)

type LocateWorker struct {
	record     *sam.Record
	cgListMap  map[string][]int
	resultChan chan string
	mappingQ   int
	minSubDep  int
	maxSubDep  int
	radius     int
	scaleFlag  bool
	err        error
}

func NewLocateWorker(record *sam.Record, cgListMap map[string][]int, resultChan chan string,
	mappingQ, minSubDep, maxSubDep, radius int, scaleFlag bool) LocateWorker {
	return LocateWorker{
		record:     record,
		cgListMap:  cgListMap,
		resultChan: resultChan,
		mappingQ:   mappingQ,
		minSubDep:  minSubDep,
		maxSubDep:  maxSubDep,
		radius:     radius,
		scaleFlag:  scaleFlag,
	}
}

func (w *LocateWorker) Task(num int) {
	record := w.record
	radius := w.radius
	//winsize := 2*radius + 1
	scaleFlag := w.scaleFlag
	if !isSecondary(record.Flags) && !isSupplementary(record.Flags) && int(record.MapQ) > w.mappingQ && matchingRatio(record) >= 0.85 {
		recordTag, err := extractRecordTag(record)
		if err != nil {
			log.Printf("extractRecordTag err:%v", err)
			w.err = err
			return
		}
		fn := recordTag.Fn
		rn := recordTag.Rn
		totalSubreadsDep := fn + rn
		if totalSubreadsDep >= int32(w.minSubDep) && totalSubreadsDep <= int32(w.maxSubDep) && fn >= 1 && rn >= 1 {
			alnRefChr := record.Ref.Name()
			alnRefStart := record.Pos
			alnRefEnd := record.End()
			log.Printf("alnRefStart:%d, alnRefEnd:%d", alnRefStart, alnRefEnd)
			cgList := w.cgListMap[alnRefChr]
			overlappingCpg := findOverlappingCpg(cgList, alnRefStart, alnRefEnd)
			log.Printf("overlappingCpg:%+v", overlappingCpg)
			if len(overlappingCpg) >= 1 {
				readFiList := recordTag.Fi
				readFpList := recordTag.Fp
				readRiList := recordTag.Ri
				readRpList := recordTag.Rp
				readSeqList := record.Seq.Expand()
				readIsReverse := isReverse(record.Flags)
				readQueryLength := len(readSeqList)
				readCigar := record.Cigar
				readName := record.Name
				parts := strings.Split(readName, "/")
				zmwname := parts[0] + "_" + parts[1]
				haplotype := recordTag.HP
				haploTypeBlock := recordTag.PS

				locatedCpgs, cpgPosOnRead := locateCpgPosOnRead(alnRefStart, readCigar, overlappingCpg)
				log.Printf("readName:%s, len(overlappingCpg):%d, len(locatedCpgs):%d", readName, len(overlappingCpg), len(locatedCpgs))
				if len(locatedCpgs) >= 1 {
					for _, cpg := range locatedCpgs {
						posOnRead := cpgPosOnRead[cpg]
						//log.Printf("readName:%s, cpg:%d, posOnRead:%d", readName, cpg, posOnRead)
						var fIPDList, fPWList, rIPDList, rPWList []float64
						var seqList []byte

						//this hifi read aligned to Ref_R_strand
						if readIsReverse {
							//(1)
							//For subreads mapped to Ref_R_strand, which carrying information when syn the complementary F_strand, key 'F' #
							//information is stored in fi, fp tags
							refFCOnFList := readQueryLength - posOnRead - 1
							leftRefFCOnFList := refFCOnFList - radius - 1
							rightRefFCOnFList := refFCOnFList + radius
							if leftRefFCOnFList > 0 && rightRefFCOnFList < readQueryLength {
								fIPDList, err = getkineticswin(readFiList[leftRefFCOnFList:rightRefFCOnFList], scaleFlag)
								if err != nil {
									log.Printf("getkineticswin err:%+v", err)
									continue
								}
								fPWList, err = getkineticswin(readFpList[leftRefFCOnFList:rightRefFCOnFList], scaleFlag)
								if err != nil {
									log.Printf("getkineticswin err:%+v", err)
									continue
								}
							}
							//(2)
							//For subreads mapped to Ref_F_strand, which carrying information when syn the complementary R_strand, key 'R' #
							//information is stored in ri, rp tags
							refRCOnRList := posOnRead + 1
							leftRefRCOnRList := refRCOnRList - radius + 1
							rightRefRCOnRList := refRCOnRList + radius
							if leftRefRCOnRList > 0 && rightRefRCOnRList < readQueryLength {
								rIPDList, err = getkineticswin(readRiList[leftRefRCOnRList:rightRefRCOnRList], scaleFlag)
								if err != nil {
									log.Printf("getkineticswin err:%+v", err)
									continue
								}
								rPWList, err = getkineticswin(readRpList[leftRefRCOnRList:rightRefRCOnRList], scaleFlag)
								if err != nil {
									log.Printf("getkineticswin err:%+v", err)
									continue
								}
								seqList = readSeqList[leftRefRCOnRList : rightRefRCOnRList+1]
								rIPDList = reverseSlice(rIPDList)
								rPWList = reverseSlice(rPWList)

							}
						} else {
							//(1)
							//For subreads mapped to Ref_F_strand, which carrying information when syn the complementary R_strand, key 'R' #
							//information is stored in fi, fp tags
							refRCOnFList := posOnRead + 1
							leftRefRCOnFList := refRCOnFList - radius - 1
							rightRefRCOnFList := refRCOnFList + radius
							if leftRefRCOnFList > 0 && rightRefRCOnFList < readQueryLength {
								rIPDList, err = getkineticswin(readFiList[leftRefRCOnFList:rightRefRCOnFList], scaleFlag)
								if err != nil {
									log.Printf("getkineticswin err:%+v", err)
									continue
								}
								rPWList, err = getkineticswin(readFpList[leftRefRCOnFList:rightRefRCOnFList], scaleFlag)
								if err != nil {
									log.Printf("getkineticswin err:%+v", err)
									continue
								}
								seqList = readSeqList[leftRefRCOnFList : rightRefRCOnFList+1]

							}
							//(2)
							//For subreads mapped to Ref_R_strand, which carrying information when syn the complementary F_strand, key 'F' #
							//information is stored in ri, rp tags
							refFCOnRList := readQueryLength - posOnRead - 1
							leftRefFCOnRList := refFCOnRList - radius - 1
							rightRefFCOnRList := refFCOnRList + radius
							if leftRefFCOnRList > 0 && rightRefFCOnRList < readQueryLength {
								fIPDList, err = getkineticswin(UInt8Slice(readRiList[leftRefFCOnRList:rightRefFCOnRList]), scaleFlag)
								if err != nil {
									log.Printf("getkineticswin err:%+v", err)
									continue
								}
								fPWList, err = getkineticswin(UInt8Slice(readRpList[leftRefFCOnRList:rightRefFCOnRList]), scaleFlag)
								if err != nil {
									log.Printf("getkineticswin err:%+v", err)
									continue
								}
								fIPDList = reverseSlice(fIPDList)
								fPWList = reverseSlice(fPWList)
							}
						}
						//输出
						ccsSeqString := string(seqList)
						fIPDPart := joinSlice(fIPDList)
						fPWPart := joinSlice(fPWList)
						rIPDPart := joinSlice(rIPDList)
						rPWPart := joinSlice(rPWList)

						if len(fIPDList) == 0 || len(fPWList) == 0 || len(rIPDList) == 0 || len(rPWList) == 0 {
							log.Printf("ipd or pw length == 0")
							continue
						}
						outputZMWLine := fmt.Sprintf("%s\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
							alnRefChr, cpg, ccsSeqString, zmwname, fn, rn, fIPDPart, fPWPart, rIPDPart, rPWPart, haplotype, haploTypeBlock)
						w.resultChan <- outputZMWLine
					}
				}
			}
		}
	}
}

type NumberSlice interface {
	Len() int
	Get(i int) float64
}

// 实现 []uint8 类型的切片
type UInt8Slice []uint8

func (u UInt8Slice) Len() int          { return len(u) }
func (u UInt8Slice) Get(i int) float64 { return float64(u[i]) }

// 实现 []float64 类型的切片
type FloatSlice []float64

func (s FloatSlice) Len() int          { return len(s) }
func (s FloatSlice) Get(i int) float64 { return s[i] }

type RecordTag struct {
	Fn int32
	Rn int32
	Fi []uint8
	Fp []uint8
	Ri []uint8
	Rp []uint8
	// HP:i:2   PS:i:11457
	HP string
	PS string
}

func reverseSlice(s []float64) []float64 {

	sort.Slice(s, func(i, j int) bool {
		return i > j
	})
	return s
}

func joinSlice(arr interface{}) string {
	var res string
	switch v := arr.(type) {
	case []uint8:
		res = formatSlice(v)
	case []float64:
		res = formatSlice(v)
	default:
		log.Printf("unknown type, arr:%v", arr)
	}
	return res
}

func formatSlice[T any](arr []T) string {
	var strs []string
	for _, num := range arr {
		strs = append(strs, fmt.Sprintf("%v", num))
	}
	return strings.Join(strs, ",")
}

func extractRecordTag(record *sam.Record) (*RecordTag, error) {
	//fn
	fnAux, exist := record.Tag([]byte{'f', 'n'})
	if !exist {
		log.Printf("record fn not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("fn not exist")
	}
	fn, ok := fnAux.Value().(int32)
	if !ok {
		log.Printf("record fn invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("fn invalid")
	}

	//rn
	rnAux, exist := record.Tag([]byte{'r', 'n'})
	if !exist {
		log.Printf("record rn not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("rn not exist")
	}
	rn, ok := rnAux.Value().(int32)
	if !ok {
		log.Printf("record rn invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("rn invalid")
	}
	fmt.Println("record name:", record.Name, "fn:", fn, " rn:", rn)

	//fi
	fiAux, exist := record.Tag([]byte{'f', 'i'})
	if !exist {
		log.Printf("record fi not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("fi not exist")
	}
	fi, ok := fiAux.Value().([]uint8)
	if !ok {
		log.Printf("record fi invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("fi invalid")
	}

	//fp
	fpAux, exist := record.Tag([]byte{'f', 'p'})
	if !exist {
		log.Printf("record fp not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("fp not exist")
	}
	fp, ok := fpAux.Value().([]uint8)
	if !ok {
		log.Printf("record fp invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("fp invalid")
	}

	//ri
	riAux, exist := record.Tag([]byte{'r', 'i'})
	if !exist {
		log.Printf("record ri not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("ri not exist")
	}
	ri, ok := riAux.Value().([]uint8)
	if !ok {
		log.Printf("record ri invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("ri invalid")
	}

	//rp
	rpAux, exist := record.Tag([]byte{'r', 'p'})
	if !exist {
		log.Printf("record rp not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("rp not exist")
	}
	rp, ok := rpAux.Value().([]uint8)
	if !ok {
		log.Printf("record rp invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("rp invalid")
	}

	// HP:i:2   PS:i:11457
	var HP, PS string
	HPAux, exist := record.Tag([]byte{'H', 'P'})
	if !exist {
		log.Printf("record HP not exist, record name:%s", record.Name)
		HP = "X"
	} else {
		HPVal, ok := HPAux.Value().(int32)
		if !ok {
			log.Printf("record HP invalid, record name:%s", record.Name)
			return nil, fmt.Errorf("HP invalid")
		}
		HP = fmt.Sprintf("%d", HPVal)
	}

	PSAux, exist := record.Tag([]byte{'P', 'S'})
	if !exist {
		log.Printf("record PS not exist, record name:%s", record.Name)
		PS = "X"
	} else {
		PSVal, ok := PSAux.Value().(int32)
		if !ok {
			log.Printf("record PS invalid, record name:%s", record.Name)
			return nil, fmt.Errorf("PS invalid")
		}
		PS = fmt.Sprintf("%d", PSVal)
	}

	recordTag := &RecordTag{
		Fn: fn,
		Rn: rn,
		Fi: fi,
		Fp: fp,
		Ri: ri,
		Rp: rp,
		HP: HP,
		PS: PS,
	}
	return recordTag, nil
}

func isReverse(flag sam.Flags) bool {
	return flag&sam.Reverse == sam.Reverse
}

func isSecondary(flag sam.Flags) bool {
	return flag&sam.Secondary == sam.Secondary
}

func isSupplementary(flag sam.Flags) bool {
	return flag&sam.Supplementary == sam.Supplementary
}

func findOverlappingCpg(chrcglist []int, refStart, refEnd int) []int {
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

func locateCpgPosOnRead(alnRefStart int, readCigar sam.Cigar, overlappingCpg []int) ([]int, map[int]int) {
	var locatedCpgs []int
	cpgPosOnRead := make(map[int]int)
	cpgBeginIdx := 0

	readLeadingFlag := 1
	readPosBlkStart := 0
	readPosBlkEnd := 0
	refWalkingPosStart := 0
	refWalkingPosEnd := 0
	for _, cigar := range readCigar {
		op := cigar.Type()
		count := cigar.Len()
		if (op == sam.CigarInsertion || op == sam.CigarSoftClipped || op == sam.CigarHardClipped ||
			op == sam.CigarPadded) && readLeadingFlag == 1 {
			//1,4 consume reads, 5_hard_clip should also be counted on reads, since the aligner doesn't change ipd/pw signal vector as the SEQ
			if op != sam.CigarPadded {
				readPosBlkStart = readPosBlkEnd
				readPosBlkEnd = readPosBlkStart + count
			}
		} else {
			//注意：下面这段逻辑有点奇怪，可以简化
			if readLeadingFlag == 1 {
				readLeadingFlag = 2
				if (op == sam.CigarMatch || op == sam.CigarDeletion || op == sam.CigarSkipped ||
					op == sam.CigarEqual || op == sam.CigarMismatch) && readLeadingFlag == 2 {
					refWalkingPosStart = alnRefStart
					refWalkingPosEnd = refWalkingPosStart + count
				}
				if (op == sam.CigarMatch || op == sam.CigarInsertion || op == sam.CigarSoftClipped ||
					op == sam.CigarEqual || op == sam.CigarMismatch) && readLeadingFlag == 2 {
					readPosBlkStart = readPosBlkEnd
					readPosBlkEnd = readPosBlkStart + count
				}
			} else {
				if op == sam.CigarMatch || op == sam.CigarDeletion || op == sam.CigarSkipped ||
					op == sam.CigarEqual || op == sam.CigarMismatch {
					refWalkingPosStart = refWalkingPosEnd
					refWalkingPosEnd = refWalkingPosStart + count
				}
				if op == sam.CigarMatch || op == sam.CigarInsertion || op == sam.CigarSoftClipped ||
					op == sam.CigarEqual || op == sam.CigarMismatch {
					readPosBlkStart = readPosBlkEnd
					readPosBlkEnd = readPosBlkStart + count
				}
			}
			if op != sam.CigarDeletion && op != sam.CigarSkipped {
				//检查当前refStart和refEnd是否包含某个cpg
				matchedCpgs, nextIdx := whichCpgMatched(overlappingCpg, cpgBeginIdx, refWalkingPosStart, refWalkingPosEnd)
				if nextIdx > len(overlappingCpg) {
					//对比结束
					return locatedCpgs, cpgPosOnRead
				} else {
					cpgBeginIdx = nextIdx
					for _, cpg := range matchedCpgs {
						lastOpNeeded := cpg - refWalkingPosStart
						if op == sam.CigarMatch || op == sam.CigarInsertion || op == sam.CigarSoftClipped ||
							op == sam.CigarEqual || op == sam.CigarMismatch {
							readPos := readPosBlkStart + lastOpNeeded
							cpgPosOnRead[cpg] = readPos - 1
							locatedCpgs = append(locatedCpgs, cpg)
						}
					}
				}
			}
		}
	}
	return locatedCpgs, cpgPosOnRead
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

func matchingRatio(record *sam.Record) float32 {
	cigars := record.Cigar
	totalLen := 0
	matchingLen := 0
	for _, cigar := range cigars {
		op := cigar.Type()
		count := cigar.Len()
		// 0:M, 1:I, 4:S, 7:=, 8:X ||||||||||||| sum of "M/I/S/=/X" is the SEQ length
		if op == sam.CigarMatch || op == sam.CigarInsertion || op == sam.CigarSoftClipped ||
			op == sam.CigarEqual || op == sam.CigarMismatch {
			totalLen += count
		}
		if op == sam.CigarMatch || op == sam.CigarEqual || op == sam.CigarMismatch {
			matchingLen += count
		}
	}
	ratio := float32(matchingLen) / float32(totalLen)
	return ratio
}

func getkineticswin(nums []uint8, scaleFlag bool) ([]float64, error) {
	if scaleFlag == false {
		var fnums []float64
		for _, num := range nums {
			fnums = append(fnums, float64(num))
		}
		return fnums, nil
	}
	//计算平均值和标准差
	if len(nums) == 0 {
		return nil, fmt.Errorf("empty input")
	}
	data := stats.LoadRawData(nums)
	mean, _ := stats.Mean(data)
	std, _ := stats.StandardDeviation(data)
	//fmt.Println("std:", std)
	//计算(x-mean)/std
	var kwin []float64
	for _, num := range nums {
		temp, err := stats.Round((float64(num)-mean)/std, 2)
		if err != nil {
			return nil, err
		}
		kwin = append(kwin, temp)
	}
	return kwin, nil
}
